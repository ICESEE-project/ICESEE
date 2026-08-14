import pytest
from types import MappingProxyType
from pathlib import Path

import numpy as np

from ICESEE.src.utils.icesee_context import (
    matlab_icesee_kwargs,
    normalize_icesee_kwargs,
)


def test_flat_context_keeps_shared_dictionary():
    context = {"nt": 10, "dt": 0.2}

    normalized = normalize_icesee_kwargs(context)

    assert normalized is context


@pytest.mark.parametrize("legacy_key", ["params", "model_kwargs", "kwargs"])
def test_legacy_nested_context_is_flattened(legacy_key):
    context = {legacy_key: {"nt": 10, "dt": 0.2}, "dt": 0.1}

    with pytest.warns(DeprecationWarning):
        normalized = normalize_icesee_kwargs(context)

    assert normalized == {"nt": 10, "dt": 0.1}
    assert legacy_key not in normalized


def test_direct_values_override_legacy_context():
    with pytest.warns(DeprecationWarning):
        normalized = normalize_icesee_kwargs(
            {"params": {"Nens": 20, "nt": 10}}, Nens=40
        )

    assert normalized == {"Nens": 40, "nt": 10}


def test_non_mapping_context_is_rejected():
    with pytest.raises(TypeError, match="icesee_kwargs must be a mapping"):
        normalize_icesee_kwargs([])


def test_read_only_mapping_becomes_mutable_dictionary():
    normalized = normalize_icesee_kwargs(MappingProxyType({"nt": 10}))

    normalized["dt"] = 0.2
    assert normalized == {"nt": 10, "dt": 0.2}


def test_matlab_view_keeps_serializable_values_and_paths():
    context = {
        "count": 3,
        "coordinates": np.array([[0.0, 1.0]]),
        "output": Path("results"),
        "names": ["h", "u"],
    }

    matlab_view = matlab_icesee_kwargs(context)

    assert matlab_view["count"] == 3
    np.testing.assert_array_equal(matlab_view["coordinates"], context["coordinates"])
    assert matlab_view["output"] == "results"
    assert matlab_view["names"] == ["h", "u"]


def test_matlab_view_excludes_runtime_only_values():
    context = {
        "server": object(),
        "callback": lambda: None,
        "optional": None,
        "dt": 0.2,
    }

    assert matlab_icesee_kwargs(context) == {"dt": 0.2}
