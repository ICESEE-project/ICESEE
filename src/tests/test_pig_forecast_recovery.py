"""Exercise the actual PIG callbacks without compiling a finite-element model."""
import ast
import __future__
from pathlib import Path
import sys
from types import SimpleNamespace as NS, ModuleType
import numpy as np
import pytest


def test_pig_retries_physics_without_changing_macro_forcing_time_or_solver_options(monkeypatch):
    path = Path(__file__).resolve().parents[2] / "applications/icepack_model/examples/idealized_pig/_icepack_native.py"
    names = {"forecast_member", "_recover_forecast_member"}
    nodes = [n for n in ast.parse(path.read_text()).body if isinstance(n,ast.FunctionDef) and n.name in names]
    class ConvergenceError(Exception): pass
    exceptions = ModuleType("firedrake.exceptions")
    exceptions.ConvergenceError = ConvergenceError
    monkeypatch.setitem(sys.modules,"firedrake.exceptions",exceptions)
    state = NS(thickness=10.,velocity=3.,surface=1.,basal_melt=None)
    def pack(): return np.array([state.thickness,state.velocity,state.surface])
    def unpack(values,**kw): state.thickness,state.velocity,state.surface = values
    def update(result):
        state.thickness,state.velocity,state.surface = result.thickness,result.velocity,result.surface
    state.registry = NS(pack_owned=pack,unpack_owned=unpack)
    state.update_from_forecast = update
    ctx = NS(forward_solver=object(),smb=0,bed=0,h0=10,Q=None)
    solves,forcing,opts = [],[],[]
    def solve(solver,h,u,smb,melt,bed,dt,h0,kwargs):
        solves.append(dt)
        if dt > .025:
            state.thickness = -999
            raise ConvergenceError("DIVERGED_MAX_IT")
        return h+dt,u,state.surface,None,None
    ns = dict(np=np,Icepack=solve,IcepackForecastFields=lambda **kw:NS(**kw),
        _shared_context=lambda *a:ctx,bind_melt_reference=lambda options,ctx:options,
        _recompute_floating=lambda *a:None,melt_forecast_mode=lambda kw:"prescribed_wrong",
        _icepack_run_kwargs=lambda *a:{},EXPERIMENT_WRONG="wrong",
        _basal_melt_for_step=lambda kw,k,*a,**b:forcing.append((kw["dt"],k)) or 20,
        icepack=NS(models=NS(IceStream=lambda **kw:None),solvers=NS(
            FlowSolver=lambda model,**kw:opts.append(kw) or object())),
        schoofFriction=None,regViscosity=None)
    exec(compile(ast.Module(body=nodes,type_ignores=[]),str(path),"exec",
                 flags=__future__.annotations.compiler_flag),ns)
    comm = NS(rank=0,allgather=lambda v:[v])
    options = dict(dt=.05,forecast_recovery=dict(max_halvings=2,min_dt=.00625))
    result = ns["forecast_member"](state,299,topology=NS(spatial_comm=comm,ensemble_slot=0),icesee_kwargs=options)
    assert solves == [.05,.025,.025]
    assert forcing == [(.05,299.),(.05,299.),(.05,299.5)]
    assert state.thickness == pytest.approx(10.05)
    assert result.thickness == state.thickness
    assert all(o == dict(dirichlet_ids=[1],diagnostic_solver_parameters={"max_iterations":150,"tolerance":1e-6}) for o in opts)
    assert options["dt"] == .05


def test_fresh_solver_receives_a_typed_threshold_with_the_same_value(monkeypatch):
    path = Path(__file__).resolve().parents[2] / "applications/icepack_model/examples/idealized_pig/_icepack_native.py"
    node = next(n for n in ast.parse(path.read_text()).body
                if isinstance(n,ast.FunctionDef) and n.name == "_icepack_run_kwargs")
    firedrake = ModuleType("firedrake")
    firedrake.Constant = lambda value:NS(value=value,ufl_shape=())
    monkeypatch.setitem(sys.modules,"firedrake",firedrake)
    namespace = {}
    exec(compile(ast.Module(body=[node],type_ignores=[]),str(path),"exec",
                 flags=__future__.annotations.compiler_flag),namespace)
    ctx = NS(Q=None,s0=None,beta0=None,A0=None,bed=None)
    options = dict(uThresh=300,GLThresh=41,hThresh=30)
    result = namespace["_icepack_run_kwargs"](ctx,options)
    assert result["uThresh"].value == 300 and options["uThresh"] == 300
    options["uThresh"] = result["uThresh"]
    assert namespace["_icepack_run_kwargs"](ctx,options)["uThresh"] is options["uThresh"]
