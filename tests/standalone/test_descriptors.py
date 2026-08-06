# -*- coding: utf-8 -*-
"""Unit tests for MutableScalar / MutableDict descriptors (C1).

Validates:
- MutableScalar: storage, type coercion, positive/non_negative validation,
  lifecycle gating, _set_unchecked bypass, class-level _ScalarRef, unset access
- MutableDict: storage, whole-dict set, item-level _set_item_unchecked,
  lifecycle gating, _set_unchecked bypass, class-level _DictRef / _DictItemRef
- Dataclass interaction: MutableDict as dataclass default value, explicit kLa,
  default kLa, class-level and instance-level access all work correctly
"""

import pytest
from dataclasses import dataclass, field
from typing import Dict, Optional


# ══════════════════════════════════════════════════════════════════════
#  Helpers
# ══════════════════════════════════════════════════════════════════════

def _running_context(label: str = "sim"):
    from PyOMES.core.lifecycle import RunContext
    return RunContext(is_running=True, label=label)


def _stopped_context(label: str = "sim"):
    from PyOMES.core.lifecycle import RunContext
    return RunContext(is_running=False, label=label)


# Simple host class — not a dataclass — used for most scalar tests.
class _ScalarHost:
    from PyOMES.control.descriptors import MutableScalar
    T_K = MutableScalar(float, positive=True)
    V_L = MutableScalar(float, positive=True)
    score = MutableScalar(float)
    volume = MutableScalar(float, non_negative=True)
    _context = None

    def __init__(self, T_K: float, V_L: float, score: float = 1.0, volume: float = 0.0):
        self.T_K = T_K
        self.V_L = V_L
        self.score = score
        self.volume = volume


# Simple host class for dict tests.
class _DictHost:
    from PyOMES.control.descriptors import MutableDict
    kLa = MutableDict(value_type=float)
    y = MutableDict(value_type=float)
    _context = None

    def __init__(self, kLa: Optional[Dict] = None, y: Optional[Dict] = None):
        self.kLa = kLa or {}
        self.y = y or {}


# ══════════════════════════════════════════════════════════════════════
#  MutableScalar — storage and coercion
# ══════════════════════════════════════════════════════════════════════

class TestMutableScalarStorage:

    def test_stores_and_retrieves_value(self):
        h = _ScalarHost(T_K=308.15, V_L=1.5)
        assert h.T_K == pytest.approx(308.15)
        assert h.V_L == pytest.approx(1.5)

    def test_set_updates_value(self):
        h = _ScalarHost(T_K=308.15, V_L=1.0)
        h.T_K = 310.0
        assert h.T_K == pytest.approx(310.0)

    def test_type_coercion_int_to_float(self):
        h = _ScalarHost(T_K=300, V_L=1)
        assert isinstance(h.T_K, float)
        assert h.T_K == 300.0

    def test_instances_are_independent(self):
        a = _ScalarHost(T_K=300.0, V_L=1.0)
        b = _ScalarHost(T_K=310.0, V_L=2.0)
        a.T_K = 305.0
        assert b.T_K == pytest.approx(310.0)

    def test_backing_store_uses_underscore_name(self):
        h = _ScalarHost(T_K=300.0, V_L=1.0)
        assert "_T_K" in h.__dict__
        assert h.__dict__["_T_K"] == pytest.approx(300.0)

    def test_unconstrained_scalar_allows_zero(self):
        h = _ScalarHost(T_K=300.0, V_L=1.0, score=0.0)
        assert h.score == 0.0

    def test_unconstrained_scalar_allows_negative(self):
        h = _ScalarHost(T_K=300.0, V_L=1.0, score=-5.0)
        assert h.score == pytest.approx(-5.0)


# ══════════════════════════════════════════════════════════════════════
#  MutableScalar — validation
# ══════════════════════════════════════════════════════════════════════

class TestMutableScalarValidation:

    def test_positive_rejects_zero(self):
        with pytest.raises(ValueError, match="positive"):
            _ScalarHost(T_K=0.0, V_L=1.0)

    def test_positive_rejects_negative(self):
        with pytest.raises(ValueError, match="positive"):
            _ScalarHost(T_K=-1.0, V_L=1.0)

    def test_positive_accepts_small_positive(self):
        h = _ScalarHost(T_K=1e-9, V_L=1.0)
        assert h.T_K > 0

    def test_non_negative_rejects_negative(self):
        with pytest.raises(ValueError, match="non-negative"):
            _ScalarHost(T_K=300.0, V_L=1.0, volume=-0.1)

    def test_non_negative_accepts_zero(self):
        h = _ScalarHost(T_K=300.0, V_L=1.0, volume=0.0)
        assert h.volume == 0.0

    def test_non_negative_accepts_positive(self):
        h = _ScalarHost(T_K=300.0, V_L=1.0, volume=5.0)
        assert h.volume == pytest.approx(5.0)


# ══════════════════════════════════════════════════════════════════════
#  MutableScalar — lifecycle gating
# ══════════════════════════════════════════════════════════════════════

class TestMutableScalarGating:

    def test_set_raises_when_running(self):
        h = _ScalarHost(T_K=300.0, V_L=1.0)
        h._context = _running_context()
        with pytest.raises(RuntimeError, match="Cannot mutate"):
            h.T_K = 310.0

    def test_set_allowed_when_stopped(self):
        h = _ScalarHost(T_K=300.0, V_L=1.0)
        h._context = _stopped_context()
        h.T_K = 310.0  # should not raise
        assert h.T_K == pytest.approx(310.0)

    def test_set_allowed_when_context_is_none(self):
        h = _ScalarHost(T_K=300.0, V_L=1.0)
        h._context = None
        h.T_K = 310.0
        assert h.T_K == pytest.approx(310.0)

    def test_set_unchecked_bypasses_gate(self):
        from PyOMES.control.descriptors import MutableScalar
        h = _ScalarHost(T_K=300.0, V_L=1.0)
        h._context = _running_context()
        descriptor = _ScalarHost.__dict__["T_K"]
        assert isinstance(descriptor, MutableScalar)
        descriptor._set_unchecked(h, 320.0)
        assert h.T_K == pytest.approx(320.0)

    def test_set_unchecked_still_validates(self):
        from PyOMES.control.descriptors import MutableScalar
        h = _ScalarHost(T_K=300.0, V_L=1.0)
        descriptor = _ScalarHost.__dict__["T_K"]
        with pytest.raises(ValueError, match="positive"):
            descriptor._set_unchecked(h, -1.0)

    def test_error_message_contains_sim_label(self):
        h = _ScalarHost(T_K=300.0, V_L=1.0)
        h._context = _running_context(label="fermentor_A")
        with pytest.raises(RuntimeError, match="fermentor_A"):
            h.T_K = 310.0


# ══════════════════════════════════════════════════════════════════════
#  MutableScalar — class-level access
# ══════════════════════════════════════════════════════════════════════

class TestMutableScalarClassLevel:

    def test_class_level_returns_descriptor_itself(self):
        from PyOMES.control.descriptors import MutableScalar
        # Class-level access returns the descriptor (SQLAlchemy pattern).
        ref = _ScalarHost.T_K
        assert isinstance(ref, MutableScalar)

    def test_class_level_knows_owner(self):
        ref = _ScalarHost.T_K
        assert ref._owner is _ScalarHost

    def test_class_level_knows_name(self):
        ref = _ScalarHost.T_K
        assert ref.name == "T_K"

    def test_class_level_different_from_instance_level(self):
        from PyOMES.control.descriptors import MutableScalar
        h = _ScalarHost(T_K=300.0, V_L=1.0)
        assert isinstance(_ScalarHost.T_K, MutableScalar)
        assert isinstance(h.T_K, float)


# ══════════════════════════════════════════════════════════════════════
#  MutableScalar — unset access
# ══════════════════════════════════════════════════════════════════════

class TestMutableScalarUnset:

    def test_unset_access_raises_attribute_error(self):
        from PyOMES.control.descriptors import MutableScalar

        class _Bare:
            x = MutableScalar(float)
            _context = None

        b = _Bare.__new__(_Bare)
        b._context = None
        with pytest.raises(AttributeError, match="has not been initialised"):
            _ = b.x


# ══════════════════════════════════════════════════════════════════════
#  MutableDict — storage
# ══════════════════════════════════════════════════════════════════════

class TestMutableDictStorage:

    def test_stores_and_retrieves_dict(self):
        h = _DictHost(kLa={"O2": 150.0, "CO2": 135.0})
        assert h.kLa["O2"] == pytest.approx(150.0)
        assert h.kLa["CO2"] == pytest.approx(135.0)

    def test_empty_dict_on_default(self):
        h = _DictHost()
        assert h.kLa == {}

    def test_instances_are_independent(self):
        a = _DictHost(kLa={"O2": 100.0})
        b = _DictHost(kLa={"O2": 200.0})
        assert a.kLa["O2"] == pytest.approx(100.0)
        assert b.kLa["O2"] == pytest.approx(200.0)

    def test_backing_store_uses_underscore_name(self):
        h = _DictHost(kLa={"O2": 150.0})
        assert "_kLa" in h.__dict__

    def test_whole_dict_replacement(self):
        h = _DictHost(kLa={"O2": 150.0})
        h.kLa = {"O2": 200.0, "CO2": 180.0}
        assert h.kLa == {"O2": 200.0, "CO2": 180.0}

    def test_returned_dict_is_the_stored_dict(self):
        # Instance-level access returns the actual dict (not a copy),
        # consistent with existing KineticGasLiquidLink.kLa behaviour.
        h = _DictHost(kLa={"O2": 150.0})
        d = h.kLa
        assert d is h.__dict__["_kLa"]


# ══════════════════════════════════════════════════════════════════════
#  MutableDict — lifecycle gating
# ══════════════════════════════════════════════════════════════════════

class TestMutableDictGating:

    def test_whole_set_raises_when_running(self):
        h = _DictHost(kLa={"O2": 150.0})
        h._context = _running_context()
        with pytest.raises(RuntimeError, match="Cannot mutate"):
            h.kLa = {"O2": 200.0}

    def test_whole_set_allowed_when_stopped(self):
        h = _DictHost(kLa={"O2": 150.0})
        h._context = _stopped_context()
        h.kLa = {"O2": 200.0}
        assert h.kLa["O2"] == pytest.approx(200.0)

    def test_set_unchecked_bypasses_gate(self):
        from PyOMES.control.descriptors import MutableDict
        h = _DictHost(kLa={"O2": 150.0})
        h._context = _running_context()
        descriptor = _DictHost.__dict__["kLa"]
        assert isinstance(descriptor, MutableDict)
        descriptor._set_unchecked(h, {"O2": 300.0})
        assert h.kLa["O2"] == pytest.approx(300.0)

    def test_set_item_unchecked_bypasses_gate(self):
        from PyOMES.control.descriptors import MutableDict
        h = _DictHost(kLa={"O2": 150.0})
        h._context = _running_context()
        descriptor = _DictHost.__dict__["kLa"]
        descriptor._set_item_unchecked(h, "O2", 250.0)
        assert h.kLa["O2"] == pytest.approx(250.0)

    def test_set_item_unchecked_adds_new_key(self):
        from PyOMES.control.descriptors import MutableDict
        h = _DictHost(kLa={"O2": 150.0})
        descriptor = _DictHost.__dict__["kLa"]
        descriptor._set_item_unchecked(h, "CO2", 135.0)
        assert h.kLa["CO2"] == pytest.approx(135.0)
        assert h.kLa["O2"] == pytest.approx(150.0)

    def test_set_item_unchecked_coerces_value_type(self):
        from PyOMES.control.descriptors import MutableDict
        h = _DictHost(kLa={})
        descriptor = _DictHost.__dict__["kLa"]
        descriptor._set_item_unchecked(h, "O2", 150)  # int → float
        assert isinstance(h.kLa["O2"], float)


# ══════════════════════════════════════════════════════════════════════
#  MutableDict — class-level access
# ══════════════════════════════════════════════════════════════════════

class TestMutableDictClassLevel:

    def test_class_level_returns_descriptor_itself(self):
        from PyOMES.control.descriptors import MutableDict
        ref = _DictHost.kLa
        assert isinstance(ref, MutableDict)

    def test_class_level_knows_owner(self):
        ref = _DictHost.kLa
        assert ref._owner is _DictHost

    def test_class_level_knows_name(self):
        ref = _DictHost.kLa
        assert ref.name == "kLa"

    def test_class_level_getitem_returns_dict_item_ref(self):
        from PyOMES.control.descriptors import _DictItemRef
        item_ref = _DictHost.kLa["O2"]
        assert isinstance(item_ref, _DictItemRef)
        assert item_ref._key == "O2"

    def test_dict_item_ref_knows_owner_and_descriptor(self):
        item_ref = _DictHost.kLa["O2"]
        assert item_ref._owner is _DictHost
        assert item_ref._descriptor.name == "kLa"

    def test_dict_item_ref_repr(self):
        r = repr(_DictHost.kLa["O2"])
        assert "O2" in r
        assert "kLa" in r

    def test_class_level_different_from_instance_level(self):
        from PyOMES.control.descriptors import MutableDict
        h = _DictHost(kLa={"O2": 150.0})
        assert isinstance(_DictHost.kLa, MutableDict)
        assert isinstance(h.kLa, dict)


# ══════════════════════════════════════════════════════════════════════
#  Dataclass interaction
# ══════════════════════════════════════════════════════════════════════

class TestDataclassInteraction:
    """Verify that MutableDict works as a @dataclass field default value.

    This is the KineticGasLiquidLink.kLa pattern: the descriptor is used
    as the default value in a dataclass annotation so that:
    - The generated __init__ accepts kLa as a keyword argument.
    - self.kLa = arg goes through MutableDict.__set__.
    - Class-level access returns a _DictRef.
    - Instance-level access returns the stored dict.
    """

    def setup_method(self):
        from PyOMES.control.descriptors import MutableDict

        @dataclass
        class _Link:
            label: str
            kLa: Dict[str, float] = MutableDict(value_type=float)
            _context: object = field(default=None, repr=False)

        self._Link = _Link

    def test_explicit_kla_is_stored(self):
        link = self._Link(label="test", kLa={"O2": 150.0, "CO2": 135.0})
        assert link.kLa["O2"] == pytest.approx(150.0)

    def test_default_kla_is_empty_dict(self):
        link = self._Link(label="test")
        assert link.kLa == {}

    def test_class_level_access_returns_descriptor(self):
        from PyOMES.control.descriptors import MutableDict
        ref = self._Link.kLa
        assert isinstance(ref, MutableDict)

    def test_instance_access_returns_dict_not_descriptor(self):
        link = self._Link(label="test", kLa={"O2": 100.0})
        assert isinstance(link.kLa, dict)
        assert link.kLa["O2"] == pytest.approx(100.0)

    def test_instances_are_independent(self):
        a = self._Link(label="a", kLa={"O2": 100.0})
        b = self._Link(label="b", kLa={"O2": 200.0})
        assert a.kLa["O2"] == pytest.approx(100.0)
        assert b.kLa["O2"] == pytest.approx(200.0)

    def test_gating_works_on_dataclass_instance(self):
        from PyOMES.core.lifecycle import RunContext
        link = self._Link(label="test", kLa={"O2": 150.0})
        link._context = RunContext(is_running=True, label="sim")
        with pytest.raises(RuntimeError, match="Cannot mutate"):
            link.kLa = {"O2": 200.0}

    def test_set_item_unchecked_works_on_dataclass_instance(self):
        from PyOMES.control.descriptors import MutableDict
        link = self._Link(label="test", kLa={"O2": 150.0})
        link._context = _running_context()
        descriptor = type(link).__dict__["kLa"]
        descriptor._set_item_unchecked(link, "O2", 300.0)
        assert link.kLa["O2"] == pytest.approx(300.0)

    def test_typed_path_form_descriptor_getitem(self):
        from PyOMES.control.descriptors import _DictItemRef
        path_stub = self._Link.kLa["O2"]
        assert isinstance(path_stub, _DictItemRef)
        assert path_stub._key == "O2"


# ══════════════════════════════════════════════════════════════════════
#  MutableScalar dataclass interaction
# ══════════════════════════════════════════════════════════════════════

class TestMutableScalarDataclassInteraction:

    def setup_method(self):
        from PyOMES.control.descriptors import MutableScalar

        @dataclass
        class _Phase:
            T_K: float = MutableScalar(float, positive=True)
            V_L: float = MutableScalar(float, positive=True)
            _context: object = field(default=None, repr=False)

        self._Phase = _Phase

    def test_explicit_values_stored(self):
        p = self._Phase(T_K=308.15, V_L=1.0)
        assert p.T_K == pytest.approx(308.15)
        assert p.V_L == pytest.approx(1.0)

    def test_class_level_returns_descriptor(self):
        from PyOMES.control.descriptors import MutableScalar
        assert isinstance(self._Phase.T_K, MutableScalar)

    def test_positive_validation_on_init(self):
        with pytest.raises(ValueError, match="positive"):
            self._Phase(T_K=0.0, V_L=1.0)

    def test_gating_on_dataclass_instance(self):
        p = self._Phase(T_K=308.15, V_L=1.0)
        p._context = _running_context()
        with pytest.raises(RuntimeError):
            p.T_K = 310.0

    def test_set_unchecked_on_dataclass_instance(self):
        from PyOMES.control.descriptors import MutableScalar
        p = self._Phase(T_K=308.15, V_L=1.0)
        p._context = _running_context()
        descriptor = type(p).__dict__["T_K"]
        descriptor._set_unchecked(p, 320.0)
        assert p.T_K == pytest.approx(320.0)


# ══════════════════════════════════════════════════════════════════════
#  Import / export
# ══════════════════════════════════════════════════════════════════════

class TestImports:

    def test_descriptors_importable_from_module(self):
        from PyOMES.control.descriptors import (
            MutableScalar, MutableDict, _DictItemRef,
        )
        assert MutableScalar is not None
        assert MutableDict is not None
        assert _DictItemRef is not None
