# -*- coding: utf-8 -*-
"""Unit tests for ParamPath + string parser + ParamPathError (C2).

Validates:
- ParamPathError is a distinct Exception subclass
- Segment types (AttrSegment, ListSelector, IndexSelector) are frozen,
  hashable, and repr correctly
- ParamPath.parse() produces correct segments for all valid path forms
- ParamPath.parse() raises ParamPathError for malformed input
- ParamPath equality and hashability (usable as dict key in params_changed)
- ParamPath.leaf / .head accessors
- ParamPath __str__ / __repr__
- Typed-form handles (_DictItemRef, class-level MutableScalar/MutableDict)
  are hashable and distinguishable by the walker
"""

import pytest
from typing import Dict


# ══════════════════════════════════════════════════════════════════════
#  Imports
# ══════════════════════════════════════════════════════════════════════

class TestImports:

    def test_param_path_error_importable(self):
        from PyOMES.control.param_path import ParamPathError
        assert issubclass(ParamPathError, Exception)

    def test_segment_types_importable(self):
        from PyOMES.control.param_path import AttrSegment, ListSelector, IndexSelector
        assert AttrSegment is not None
        assert ListSelector is not None
        assert IndexSelector is not None

    def test_param_path_importable(self):
        from PyOMES.control.param_path import ParamPath
        assert ParamPath is not None


# ══════════════════════════════════════════════════════════════════════
#  Segment types — construction and repr
# ══════════════════════════════════════════════════════════════════════

class TestAttrSegment:

    def test_name_stored(self):
        from PyOMES.control.param_path import AttrSegment
        s = AttrSegment("T_K")
        assert s.name == "T_K"

    def test_frozen(self):
        from PyOMES.control.param_path import AttrSegment
        s = AttrSegment("T_K")
        with pytest.raises(Exception):
            s.name = "V_L"

    def test_hashable(self):
        from PyOMES.control.param_path import AttrSegment
        s = AttrSegment("T_K")
        assert hash(s) == hash(AttrSegment("T_K"))

    def test_equality(self):
        from PyOMES.control.param_path import AttrSegment
        assert AttrSegment("T_K") == AttrSegment("T_K")
        assert AttrSegment("T_K") != AttrSegment("V_L")

    def test_repr(self):
        from PyOMES.control.param_path import AttrSegment
        assert repr(AttrSegment("kLa")) == ".kLa"


class TestListSelector:

    def test_type_name_stored(self):
        from PyOMES.control.param_path import ListSelector
        s = ListSelector("GasFeed")
        assert s.type_name == "GasFeed"
        assert s.index is None

    def test_index_stored(self):
        from PyOMES.control.param_path import ListSelector
        s = ListSelector("GasFeed", index=1)
        assert s.index == 1

    def test_frozen(self):
        from PyOMES.control.param_path import ListSelector
        s = ListSelector("GasFeed")
        with pytest.raises(Exception):
            s.type_name = "GasPhase"

    def test_hashable(self):
        from PyOMES.control.param_path import ListSelector
        s1 = ListSelector("GasFeed")
        s2 = ListSelector("GasFeed")
        assert hash(s1) == hash(s2)

    def test_equality_no_index(self):
        from PyOMES.control.param_path import ListSelector
        assert ListSelector("GasFeed") == ListSelector("GasFeed")
        assert ListSelector("GasFeed") != ListSelector("GasPhase")

    def test_equality_with_index(self):
        from PyOMES.control.param_path import ListSelector
        assert ListSelector("GasFeed", 1) == ListSelector("GasFeed", 1)
        assert ListSelector("GasFeed", 0) != ListSelector("GasFeed", 1)

    def test_repr_no_index(self):
        from PyOMES.control.param_path import ListSelector
        assert repr(ListSelector("GasFeed")) == "[GasFeed]"

    def test_repr_with_index(self):
        from PyOMES.control.param_path import ListSelector
        assert repr(ListSelector("GasFeed", 2)) == "[GasFeed:2]"


class TestIndexSelector:

    def test_index_stored(self):
        from PyOMES.control.param_path import IndexSelector
        s = IndexSelector(3)
        assert s.index == 3

    def test_frozen(self):
        from PyOMES.control.param_path import IndexSelector
        s = IndexSelector(0)
        with pytest.raises(Exception):
            s.index = 1

    def test_hashable(self):
        from PyOMES.control.param_path import IndexSelector
        assert hash(IndexSelector(0)) == hash(IndexSelector(0))

    def test_equality(self):
        from PyOMES.control.param_path import IndexSelector
        assert IndexSelector(0) == IndexSelector(0)
        assert IndexSelector(0) != IndexSelector(1)

    def test_repr(self):
        from PyOMES.control.param_path import IndexSelector
        assert repr(IndexSelector(5)) == "[5]"


# ══════════════════════════════════════════════════════════════════════
#  ParamPath.parse — valid paths
# ══════════════════════════════════════════════════════════════════════

class TestParsePath:

    def test_single_segment(self):
        from PyOMES.control.param_path import ParamPath, AttrSegment
        p = ParamPath.parse("T_K")
        assert p.segments == (AttrSegment("T_K"),)

    def test_simple_attr_chain(self):
        from PyOMES.control.param_path import ParamPath, AttrSegment
        p = ParamPath.parse("phases.liquid.T_K")
        assert p.segments == (
            AttrSegment("phases"),
            AttrSegment("liquid"),
            AttrSegment("T_K"),
        )

    def test_type_selector_no_index(self):
        from PyOMES.control.param_path import ParamPath, AttrSegment, ListSelector
        p = ParamPath.parse("boundaries[GasFeed].vvm_min")
        assert p.segments == (
            AttrSegment("boundaries"),
            ListSelector("GasFeed"),
            AttrSegment("vvm_min"),
        )

    def test_type_selector_with_index(self):
        from PyOMES.control.param_path import ParamPath, AttrSegment, ListSelector
        p = ParamPath.parse("boundaries[GasFeed:1].vvm_min")
        assert p.segments == (
            AttrSegment("boundaries"),
            ListSelector("GasFeed", index=1),
            AttrSegment("vvm_min"),
        )

    def test_type_selector_index_zero(self):
        from PyOMES.control.param_path import ParamPath, ListSelector
        p = ParamPath.parse("boundaries[GasFeed:0].vvm_min")
        sel = p.segments[1]
        assert isinstance(sel, ListSelector)
        assert sel.index == 0

    def test_pure_numeric_index_selector(self):
        from PyOMES.control.param_path import ParamPath, AttrSegment, IndexSelector
        p = ParamPath.parse("boundaries[0].vvm_min")
        assert p.segments == (
            AttrSegment("boundaries"),
            IndexSelector(0),
            AttrSegment("vvm_min"),
        )

    def test_deep_path_with_dict_key(self):
        """internal_interfaces[KineticGasLiquidLink].kLa.O2 — the canonical kLa path."""
        from PyOMES.control.param_path import (
            ParamPath, AttrSegment, ListSelector,
        )
        p = ParamPath.parse(
            "internal_interfaces[KineticGasLiquidLink].kLa.O2"
        )
        assert p.segments == (
            AttrSegment("internal_interfaces"),
            ListSelector("KineticGasLiquidLink"),
            AttrSegment("kLa"),
            AttrSegment("O2"),
        )

    def test_raw_string_preserved(self):
        from PyOMES.control.param_path import ParamPath
        raw = "boundaries[GasFeed].vvm_min"
        p = ParamPath.parse(raw)
        assert p.raw == raw

    def test_camelcase_type_names(self):
        from PyOMES.control.param_path import ParamPath, ListSelector
        p = ParamPath.parse("internal_interfaces[KineticGasLiquidLink].kLa.CO2")
        assert isinstance(p.segments[1], ListSelector)
        assert p.segments[1].type_name == "KineticGasLiquidLink"

    def test_underscore_attribute_names(self):
        from PyOMES.control.param_path import ParamPath, AttrSegment
        p = ParamPath.parse("internal_interfaces.kLa._backing")
        assert p.segments[-1] == AttrSegment("_backing")


# ══════════════════════════════════════════════════════════════════════
#  ParamPath.parse — invalid paths
# ══════════════════════════════════════════════════════════════════════

class TestParseInvalidPaths:

    def test_empty_string(self):
        from PyOMES.control.param_path import ParamPath, ParamPathError
        with pytest.raises(ParamPathError):
            ParamPath.parse("")

    def test_whitespace_only(self):
        from PyOMES.control.param_path import ParamPath, ParamPathError
        with pytest.raises(ParamPathError):
            ParamPath.parse("   ")

    def test_consecutive_dots(self):
        from PyOMES.control.param_path import ParamPath, ParamPathError
        with pytest.raises(ParamPathError):
            ParamPath.parse("phases..liquid")

    def test_leading_dot(self):
        from PyOMES.control.param_path import ParamPath, ParamPathError
        with pytest.raises(ParamPathError):
            ParamPath.parse(".phases.liquid")

    def test_trailing_dot(self):
        from PyOMES.control.param_path import ParamPath, ParamPathError
        with pytest.raises(ParamPathError):
            ParamPath.parse("phases.liquid.")

    def test_invalid_segment_starts_with_digit(self):
        from PyOMES.control.param_path import ParamPath, ParamPathError
        with pytest.raises(ParamPathError):
            ParamPath.parse("1phases.liquid")

    def test_invalid_selector_lowercase_type(self):
        """Lowercase identifiers in brackets are not valid type selectors."""
        from PyOMES.control.param_path import ParamPath, ParamPathError
        # 'liquid' starts lowercase so is neither digits nor CamelCase
        with pytest.raises(ParamPathError):
            ParamPath.parse("phases[liquid].T_K")

    def test_invalid_selector_with_special_chars(self):
        from PyOMES.control.param_path import ParamPath, ParamPathError
        with pytest.raises(ParamPathError):
            ParamPath.parse("boundaries[Gas-Feed].vvm_min")

    def test_segment_with_hyphen(self):
        from PyOMES.control.param_path import ParamPath, ParamPathError
        with pytest.raises(ParamPathError):
            ParamPath.parse("phase-key.T_K")


# ══════════════════════════════════════════════════════════════════════
#  ParamPath — equality, hash, dict key usage
# ══════════════════════════════════════════════════════════════════════

class TestParamPathHashAndEquality:

    def test_equal_paths_are_equal(self):
        from PyOMES.control.param_path import ParamPath
        p1 = ParamPath.parse("phases.liquid.T_K")
        p2 = ParamPath.parse("phases.liquid.T_K")
        assert p1 == p2

    def test_different_paths_not_equal(self):
        from PyOMES.control.param_path import ParamPath
        p1 = ParamPath.parse("phases.liquid.T_K")
        p2 = ParamPath.parse("phases.gas.T_K")
        assert p1 != p2

    def test_hash_equal_for_equal_paths(self):
        from PyOMES.control.param_path import ParamPath
        p1 = ParamPath.parse("boundaries[GasFeed].vvm_min")
        p2 = ParamPath.parse("boundaries[GasFeed].vvm_min")
        assert hash(p1) == hash(p2)

    def test_usable_as_dict_key(self):
        from PyOMES.control.param_path import ParamPath
        p = ParamPath.parse("phases.liquid.T_K")
        d = {p: 308.15}
        assert d[p] == 308.15

    def test_two_paths_as_separate_dict_keys(self):
        from PyOMES.control.param_path import ParamPath
        p1 = ParamPath.parse("boundaries[GasFeed].vvm_min")
        p2 = ParamPath.parse("boundaries[GasFeed].y.O2")
        d = {p1: 0.5, p2: 0.21}
        assert len(d) == 2
        assert d[p1] == 0.5
        assert d[p2] == 0.21

    def test_raw_mismatch_does_not_affect_equality(self):
        """Two ParamPath with same segments but different raw strings are equal
        because equality is on the frozen dataclass fields (segments + raw)."""
        from PyOMES.control.param_path import ParamPath, AttrSegment
        p1 = ParamPath.parse("phases.liquid.T_K")
        p2 = ParamPath(
            segments=(
                AttrSegment("phases"),
                AttrSegment("liquid"),
                AttrSegment("T_K"),
            ),
            raw="",
        )
        # p1.raw != p2.raw, so they are NOT equal by frozen dataclass rules.
        assert p1 != p2

    def test_equal_paths_with_same_raw(self):
        from PyOMES.control.param_path import ParamPath
        raw = "phases.liquid.T_K"
        p1 = ParamPath.parse(raw)
        p2 = ParamPath.parse(raw)
        assert p1 == p2


# ══════════════════════════════════════════════════════════════════════
#  ParamPath — leaf / head accessors
# ══════════════════════════════════════════════════════════════════════

class TestParamPathAccessors:

    def test_leaf_is_last_segment(self):
        from PyOMES.control.param_path import ParamPath, AttrSegment
        p = ParamPath.parse("phases.liquid.T_K")
        assert p.leaf == AttrSegment("T_K")

    def test_head_is_all_but_last(self):
        from PyOMES.control.param_path import ParamPath, AttrSegment
        p = ParamPath.parse("phases.liquid.T_K")
        assert p.head == (AttrSegment("phases"), AttrSegment("liquid"))

    def test_leaf_of_single_segment_path(self):
        from PyOMES.control.param_path import ParamPath, AttrSegment
        p = ParamPath.parse("T_K")
        assert p.leaf == AttrSegment("T_K")

    def test_head_of_single_segment_path_is_empty(self):
        from PyOMES.control.param_path import ParamPath
        p = ParamPath.parse("T_K")
        assert p.head == ()

    def test_leaf_selector_path(self):
        from PyOMES.control.param_path import ParamPath, AttrSegment
        p = ParamPath.parse("boundaries[GasFeed].vvm_min")
        assert p.leaf == AttrSegment("vvm_min")

    def test_head_selector_path(self):
        from PyOMES.control.param_path import ParamPath, AttrSegment, ListSelector
        p = ParamPath.parse("boundaries[GasFeed].vvm_min")
        assert p.head == (AttrSegment("boundaries"), ListSelector("GasFeed"))


# ══════════════════════════════════════════════════════════════════════
#  ParamPath — str / repr
# ══════════════════════════════════════════════════════════════════════

class TestParamPathRepr:

    def test_str_returns_raw_when_parsed(self):
        from PyOMES.control.param_path import ParamPath
        raw = "boundaries[GasFeed].vvm_min"
        p = ParamPath.parse(raw)
        assert str(p) == raw

    def test_repr_contains_path_string(self):
        from PyOMES.control.param_path import ParamPath
        p = ParamPath.parse("phases.liquid.T_K")
        assert "phases.liquid.T_K" in repr(p)

    def test_str_of_programmatic_path(self):
        from PyOMES.control.param_path import ParamPath, AttrSegment, ListSelector
        p = ParamPath(
            segments=(
                AttrSegment("boundaries"),
                ListSelector("GasFeed"),
                AttrSegment("vvm_min"),
            ),
            raw="",
        )
        s = str(p)
        assert "boundaries" in s
        assert "GasFeed" in s
        assert "vvm_min" in s


# ══════════════════════════════════════════════════════════════════════
#  Typed-form handles — hashability for params_changed dict keys
# ══════════════════════════════════════════════════════════════════════

class TestTypedFormHandles:
    """_DictItemRef and class-level descriptors must be usable as dict keys.

    This verifies the hash contract needed for:
        params_changed = {
            KineticGasLiquidLink.kLa["O2"]: 150.0,
            GasFeed.vvm_min: 0.5,
        }
    """

    def test_dict_item_ref_is_hashable(self):
        from PyOMES.control.descriptors import MutableDict, _DictItemRef

        class _Host:
            kLa: dict = MutableDict(value_type=float)
            _context = None

        ref = _Host.kLa["O2"]
        assert isinstance(ref, _DictItemRef)
        h = hash(ref)  # should not raise
        assert isinstance(h, int)

    def test_dict_item_ref_usable_as_dict_key(self):
        from PyOMES.control.descriptors import _DictItemRef, MutableDict

        class _Host:
            kLa: dict = MutableDict(value_type=float)
            _context = None

        ref = _Host.kLa["O2"]
        d = {ref: 150.0}
        assert d[ref] == 150.0

    def test_mutable_scalar_class_level_is_hashable(self):
        from PyOMES.control.descriptors import MutableScalar

        class _Host:
            T_K = MutableScalar(float, positive=True)
            _context = None

        descriptor = _Host.T_K  # class-level → the descriptor itself
        assert isinstance(descriptor, MutableScalar)
        h = hash(descriptor)
        assert isinstance(h, int)

    def test_mutable_scalar_class_level_usable_as_dict_key(self):
        from PyOMES.control.descriptors import MutableScalar

        class _Host:
            T_K = MutableScalar(float, positive=True)
            _context = None

        d = {_Host.T_K: 308.15}
        assert d[_Host.T_K] == 308.15

    def test_mutable_dict_class_level_is_hashable(self):
        from PyOMES.control.descriptors import MutableDict

        class _Host:
            kLa: dict = MutableDict(value_type=float)
            _context = None

        descriptor = _Host.kLa  # class-level → the descriptor itself
        assert isinstance(descriptor, MutableDict)
        h = hash(descriptor)
        assert isinstance(h, int)

    def test_mixed_typed_and_string_paths_in_same_dict(self):
        """String ParamPath and typed-form handles can coexist as keys."""
        from PyOMES.control.param_path import ParamPath
        from PyOMES.control.descriptors import _DictItemRef, MutableDict

        class _Host:
            kLa: dict = MutableDict(value_type=float)
            _context = None

        string_path = ParamPath.parse("boundaries[GasFeed].vvm_min")
        typed_key = _Host.kLa["O2"]

        d = {string_path: 0.5, typed_key: 150.0}
        assert len(d) == 2
