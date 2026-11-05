"""Synthetic binding to the shared whole-member native runtime.

One spatial rank per member is supported. Filtering, localization, parameter
inference and publication use the same operator as file-backed Mode 2.
"""
import numpy as np
import firedrake
from . import _icepack_model, _icepack_enkf
from ICESEE.src.parallelization.distributed_whole_member_adapter import (
    WholeMemberNativeAdapter, WholeMemberState,
)
from ICESEE.src.parallelization.distributed_mode3_registry import register_native_application
from ICESEE.src.parallelization.native_file_analysis import FileBackedNativeAnalysis
from ICESEE.src.parallelization._mpi_ensemble_intialization import _member_initialization_context
from ICESEE.src.parallelization._mpi_forecast_functions import advance_ensemble_member
from ICESEE.src.run_model_da._error_generation import generate_initial_member_increment


def _setup(options):
    return options


def _initialize(member_id, *, topology, icesee_kwargs):
    context = dict(icesee_kwargs)
    context["comm"] = firedrake.COMM_SELF
    context = _icepack_model.initialize_model(**context)
    context["solver"] = context["solver_weertman"]
    context["a_nuged"] = firedrake.Function(context["Q"]).interpolate(
        firedrake.Constant(float(context["a_in_p"]))
        + firedrake.Constant(float(context["da_p"]))*context["x"]/context["Lx"])
    context.update(ens_id=member_id, model_module=_icepack_enkf, sub_rank=0)
    with _member_initialization_context(context, member_id) as member_options:
        variables = _icepack_enkf.initialize_ensemble(member_id, **member_options)
    order = list(context["vec_inputs"])
    fields = {name: np.asarray(variables[name]).copy() for name in order}
    hdim = len(fields[order[0]])
    increment, _ = generate_initial_member_increment(hdim, context, member_id, int(context["nd"]))
    offset = 0
    for name in order:
        count = fields[name].size
        fields[name] += increment[offset:offset+count]
        offset += count
    return WholeMemberState(fields, vec_inputs=order, model_context=context,
                            layout_id="synthetic-ice-stream-whole-member-v1")


def _forecast(state, timestep, *, topology, icesee_kwargs):
    context = state.model_context
    # Preserve each member mesh, solver and physical fields.
    for name in ("k", "km", "alpha", "rho", "obs_index", "number_obs_instants",
                 "m_obs", "base_seed", "data_path", "_modelrun_datasets",
                 "process_noise_schedule", "execution_parity_trace"):
        if name in icesee_kwargs:
            context[name] = icesee_kwargs[name]
    packed = state.registry.pack_owned()
    indices = {}
    offset = 0
    for name in state.vec_inputs:
        count = state.registry.field(name).global_size
        indices[name] = np.arange(offset,offset+count)
        offset += count
    updated = advance_ensemble_member(packed,context["ens_id"],context,indices,sub_rank=0)
    return {name: updated[rows] for name,rows in indices.items()}


register_native_application(
    "icepack", application="synthetic_ice_stream",
    adapter=WholeMemberNativeAdapter(initialize_member=_initialize,forecast_member=_forecast),
    setup=_setup, analysis_cycle=FileBackedNativeAnalysis(),
    supported_features=("localization_flag", "inference_plugin_enabled", "inversion_flag"),
)
