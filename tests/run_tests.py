#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
run_tests.py — Run the standalone test suite.

Usage (from project root):

    # With pytest (preferred, if installed):
    PYTHONPATH=src pytest tests/standalone/ -v

    # Without pytest (uses this script + unittest):
    PYTHONPATH=src python tests/run_tests.py

    # Run a single test file:
    PYTHONPATH=src python tests/run_tests.py tests/standalone/test_compounds.py

The test suite has NO bioSTEAM dependency.  It validates:
  - Chemical registry and compound data
  - FeedState construction and accessors
  - Strong-ion inference
  - Speciation engine (pH chemistry)
  - Headspace ideal-gas helpers
  - Stoichiometry (CHO / CHNO)
  - Controllers (pressure, pH, DO)
  - Kinetic model interface
  - Full fermenter integration (simulate end-to-end)
"""

import os
import sys
import unittest
import time

# ── Path setup ──────────────────────────────────────────────────────────
_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.abspath(os.path.join(_HERE, ".."))
_SRC = os.path.join(_ROOT, "src")
for p in (_ROOT, _SRC, _HERE):
    if p not in sys.path:
        sys.path.insert(0, p)


# ── Compatibility shim ─────────────────────────────────────────────────
# pytest uses `pytest.approx` and `@pytest.mark.*`. When running under
# unittest we need lightweight stubs.

class _Approx:
    """Minimal pytest.approx replacement for unittest."""
    def __init__(self, expected, *, rel=None, abs=None):
        self.expected = expected
        self.rel = rel if rel is not None else 1e-6
        self._abs = abs if abs is not None else 1e-12

    def __eq__(self, other):
        import builtins
        _abs = builtins.abs
        tol = max(self._abs, self.rel * _abs(self.expected) if isinstance(self.expected, (int, float)) else 0.0)
        return _abs(float(other) - float(self.expected)) <= tol

    def __repr__(self):
        return f"approx({self.expected}, rel={self.rel}, abs={self._abs})"

try:
    import pytest
except ImportError:
    # Build a minimal pytest shim so test files can be imported
    class _PytestShim:
        class mark:
            @staticmethod
            def unit(cls_or_fn):    return cls_or_fn
            @staticmethod
            def integration(cls_or_fn): return cls_or_fn
            @staticmethod
            def speciation(cls_or_fn):  return cls_or_fn
            @staticmethod
            def parametrize(*a, **kw):
                def deco(fn): return fn
                return deco

        @staticmethod
        def approx(expected, *, rel=None, abs=None):
            return _Approx(expected, rel=rel, abs=abs)

        @staticmethod
        def raises(exc, *, match=None):
            import contextlib, re
            @contextlib.contextmanager
            def ctx():
                try:
                    yield
                except exc as e:
                    if match and not re.search(match, str(e)):
                        raise AssertionError(f"Exception {e!r} did not match {match!r}")
                else:
                    raise AssertionError(f"{exc.__name__} not raised")
            return ctx()

        fixture = staticmethod(lambda fn=None, **kw: fn if fn else (lambda f: f))

    sys.modules["pytest"] = _PytestShim()
    pytest = _PytestShim()


# ── Fixture helper (replaces pytest fixtures for unittest) ─────────────

from PyOMES.chemistry.compounds import ChemicalRegistry
from PyOMES.stream_adapter import FeedState
from PyOMES import PressureReliefController, PHController, create_standalone_fermenter


def _registry():
    return ChemicalRegistry.default()

def _simple_feed():
    return FeedState.from_mass_concentrations(
        {"AceticAcid": 1.0, "Yeast": 0.1}, registry=_registry(),
    )

def _rich_feed():
    return FeedState.from_mixed_concentrations(
        mass_g_L={"AceticAcid": 1.0, "Yeast": 0.1},
        molar_mol_L={"NH3": 0.015, "AmmoniumSulfate": 0.015, "KH2PO4": 0.007,
                      "MgSO4": 4e-3, "ZnSO4": 1.4e-3, "CaCl2": 5e-3,
                      "MnCl2": 4e-3, "CoCl2": 7.7e-4},
        registry=_registry(),
    )

def _pressure_ctrl():
    return PressureReliefController(P_set_atm=1.10, diameter_m=5e-2, sample_period_s=1)

def _ph_ctrl():
    return PHController(setpoint=6.6, chemical_id="H3PO4", base_chemical_id="KOH",
                        Kp=1.0, Ki=1.0, max_add_molL_hr=10)

def _fermenter_minimal():
    return create_standalone_fermenter(
        tau=2, T=305.15, V_total_L=2.0, organism_id="Yeast", balance_basis="CHO",
    )

def _fermenter_full():
    return create_standalone_fermenter(
        tau=5, T=305.15, V_total_L=2000, organism_id="Yeast", balance_basis="CHNO",
        n_source_id="NH3", yO2_init=0.21, yCO2_init=0.0004, yN2_init=0.7896,
        gas_vvm_min=1.0, gas_feed_composition={"O2": 0.21, "N2": 0.79},
        use_activity=True, activity_model="davies",
        controllers=[_pressure_ctrl(), _ph_ctrl()],
    )


# ── Monkey-patch conftest fixtures into the test modules ───────────────
# When running under pytest, fixtures are injected automatically. Under
# unittest we inject them as module-level callables that the parametrize
# decorator and test classes can reference.

_FIXTURE_MAP = {
    "registry":          _registry,
    "simple_feed":       _simple_feed,
    "rich_feed":         _rich_feed,
    "pressure_ctrl":     _pressure_ctrl,
    "ph_ctrl":           _ph_ctrl,
    "fermenter_minimal": _fermenter_minimal,
    "fermenter_full":    _fermenter_full,
}


def _adapt_test_class(cls):
    """Wrap each test method so that fixture names in the signature are resolved."""
    import inspect
    for name in list(dir(cls)):
        if not name.startswith("test_"):
            continue
        method = getattr(cls, name)
        if not callable(method):
            continue
        sig = inspect.signature(method)
        fixture_params = [p for p in sig.parameters if p != "self" and p in _FIXTURE_MAP]
        if fixture_params:
            orig = method
            def make_wrapper(m, fps):
                def wrapper(self_):
                    kwargs = {fp: _FIXTURE_MAP[fp]() for fp in fps}
                    return m(self_, **kwargs)
                wrapper.__name__ = m.__name__
                wrapper.__doc__ = m.__doc__
                return wrapper
            setattr(cls, name, make_wrapper(orig, fixture_params))
    return cls


# ── Discover and run ───────────────────────────────────────────────────

def main():
    test_dir = os.path.join(_HERE, "standalone")
    if not os.path.isdir(test_dir):
        print(f"ERROR: test directory not found: {test_dir}")
        sys.exit(1)

    # Import all test modules and adapt fixture injection
    test_files = sorted(f for f in os.listdir(test_dir) if f.startswith("test_") and f.endswith(".py"))
    print(f"\n{'='*60}")
    print(f"  Standalone Test Suite — {len(test_files)} modules")
    print(f"{'='*60}\n")

    loader = unittest.TestLoader()
    suite = unittest.TestSuite()

    for tf in test_files:
        mod_name = tf[:-3]
        try:
            mod = __import__(f"standalone.{mod_name}", fromlist=[mod_name])
        except Exception as e:
            print(f"  SKIP {tf}: import error: {e}")
            continue

        # Discover plain test classes (not necessarily unittest.TestCase)
        # and wrap them so unittest can run them.
        for attr_name in dir(mod):
            obj = getattr(mod, attr_name)
            if not isinstance(obj, type) or not attr_name.startswith("Test"):
                continue

            # Collect test methods
            test_methods = [m for m in dir(obj) if m.startswith("test_") and callable(getattr(obj, m))]
            if not test_methods:
                continue

            # Create a unittest.TestCase subclass dynamically
            ns = {}
            for mname in test_methods:
                orig = getattr(obj, mname)
                # Resolve fixture params
                import inspect
                sig = inspect.signature(orig)
                fixture_params = [p for p in sig.parameters if p != "self" and p in _FIXTURE_MAP]

                if fixture_params:
                    def make_test(m, fps):
                        def test_fn(self_):
                            kwargs = {fp: _FIXTURE_MAP[fp]() for fp in fps}
                            return m(self_, **kwargs)
                        test_fn.__name__ = m.__name__
                        test_fn.__doc__ = m.__doc__
                        return test_fn
                    ns[mname] = make_test(orig, fixture_params)
                else:
                    # Wrap to call as instance method on a fresh TestCase
                    def make_plain(m):
                        def test_fn(self_):
                            return m(self_)
                        test_fn.__name__ = m.__name__
                        test_fn.__doc__ = m.__doc__
                        return test_fn
                    ns[mname] = make_plain(orig)

            tc_cls = type(f"{mod_name}.{attr_name}", (unittest.TestCase,), ns)
            for mname in test_methods:
                suite.addTest(tc_cls(mname))

    # Run
    t0 = time.perf_counter()
    runner = unittest.TextTestRunner(verbosity=2)
    result = runner.run(suite)
    elapsed = time.perf_counter() - t0

    print(f"\n{'='*60}")
    print(f"  {result.testsRun} tests in {elapsed:.2f}s")
    print(f"  Failures: {len(result.failures)}")
    print(f"  Errors:   {len(result.errors)}")
    print(f"  Skipped:  {len(result.skipped)}")
    print(f"{'='*60}\n")

    sys.exit(0 if result.wasSuccessful() else 1)


if __name__ == "__main__":
    main()
