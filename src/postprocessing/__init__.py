"""ICESEE post-processing: load a run's data-assimilation outputs and draw the
examples' scientific result views.

The figures are the ones each example's ``read_results.ipynb`` draws; the
notebooks call these functions, and so can any application. See
:mod:`.outputs` (loading), :mod:`.views` (figures and tables),
:mod:`.contract` (the per-example ``results.yaml``) and :mod:`.render`
(saving, and the command-line entry point).
"""
from .contract import ContractError, ResultsContract, load_contract, parse_contract
from .outputs import DAOutputs, ResultsError, find_data_dir, load_da_outputs, load_params

__all__ = [
    "ContractError", "DAOutputs", "ResultsContract", "ResultsError", "find_data_dir",
    "load_contract", "load_da_outputs", "load_params", "parse_contract",
]
