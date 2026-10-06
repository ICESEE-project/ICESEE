# ==============================================================================
# @des: Structural/unit test for ISSM basal_friction_variation's state
# layout, using resolve_variable_block_sizes / generate_initial_member_increment
# directly against this application's actual configured
# vec_inputs/scalar_inputs/var_nd values (real MATLAB/ISSM dependencies are
# unavailable in this environment -- see the 6 pre-existing ISSM-data-file
# test failures -- so this cannot be a real integration run).
#
# basal_friction_variation/params.yaml sets scalar_inputs=['IceVolume',
# 'IceVolumeAboveFloatation'], but neither name appears in its own
# vec_inputs=['Thickness','Surface','Vx','Vy','Bed','FrictionCoefficient']
# -- those two scalar_inputs entries are diagnostic quantities read
# elsewhere (_issm_enkf.py), not state-vector blocks. run_da_issm.py
# separately always builds var_nd = {var: (1 if var in scalar_inputs else
# variable_size) for var in vec_inputs} -- since no vec_inputs name is in
# scalar_inputs here, every entry resolves to the same variable_size, i.e.
# this application's layout is (and always was) effectively uniform. These
# tests prove the new heterogeneous-layout-aware code preserves that
# exactly -- a non-regression check for this application, not a new
# feature it needs.
# ==============================================================================
import numpy as np
import pytest

from ICESEE.src.run_model_da._error_generation import (
    generate_initial_member_increment,
    resolve_variable_block_sizes,
)

_VEC_INPUTS = ["Thickness", "Surface", "Vx", "Vy", "Bed", "FrictionCoefficient"]
_SCALAR_INPUTS = ["IceVolume", "IceVolumeAboveFloatation"]  # not in _VEC_INPUTS


def test_scalar_inputs_that_are_not_state_variables_are_a_no_op():
    # Matches run_da_issm.py's own var_nd construction for this app exactly.
    variable_size = 25  # arbitrary "number of mesh vertices" stand-in
    var_nd = {
        var: (1 if var in _SCALAR_INPUTS else variable_size) for var in _VEC_INPUTS
    }
    assert set(var_nd.values()) == {variable_size}  # confirmed uniform, as expected

    vector_size = variable_size * len(_VEC_INPUTS)
    sizes = resolve_variable_block_sizes(
        _VEC_INPUTS, vector_size, scalar_inputs=_SCALAR_INPUTS, var_nd=var_nd
    )
    assert sizes == [variable_size] * len(_VEC_INPUTS)


def test_generate_initial_member_increment_matches_uniform_case_for_this_layout():
    variable_size = 10
    var_nd = {
        var: (1 if var in _SCALAR_INPUTS else variable_size) for var in _VEC_INPUTS
    }
    vector_size = variable_size * len(_VEC_INPUTS)

    kwargs_heterogeneous_aware = {
        "total_state_param_vars": len(_VEC_INPUTS),
        "vec_inputs": _VEC_INPUTS,
        "scalar_inputs": _SCALAR_INPUTS,
        "var_nd": var_nd,
        "sig_Q": [0.02] * len(_VEC_INPUTS),
        "Lx": 10.0,
        "Ly": 10.0,
        "base_seed": 42,
        "random_field_method": "fft",
    }
    kwargs_uniform_only = {
        k: v for k, v in kwargs_heterogeneous_aware.items()
        if k not in ("scalar_inputs", "var_nd")
    }

    inc_a, raw_a = generate_initial_member_increment(
        variable_size, kwargs_heterogeneous_aware, ensemble_id=3, vector_size=vector_size
    )
    inc_b, raw_b = generate_initial_member_increment(
        variable_size, kwargs_uniform_only, ensemble_id=3, vector_size=vector_size
    )
    np.testing.assert_array_equal(inc_a, inc_b)
    np.testing.assert_array_equal(raw_a, raw_b)


def test_would_correctly_handle_a_genuine_scalar_state_variable_if_configured():
    # If a future basal_friction_variation config ever adds a genuine
    # scalar STATE variable (one that IS in vec_inputs), the layout
    # resolves correctly rather than silently staying uniform.
    vec_inputs = _VEC_INPUTS + ["SomeScalarDiagnosticInState"]
    scalar_inputs = _SCALAR_INPUTS + ["SomeScalarDiagnosticInState"]
    variable_size = 10
    vector_size = variable_size * len(_VEC_INPUTS) + 1

    sizes = resolve_variable_block_sizes(
        vec_inputs, vector_size, scalar_inputs=scalar_inputs
    )
    assert sizes == [variable_size] * len(_VEC_INPUTS) + [1]
