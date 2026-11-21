import numpy as np
import pytest
from ICESEE.src.utils.analysis_bounds import bound_analysis_rows
from ICESEE.src.parallelization._mpi_analysis_functions import apply_analysis_controls_local


def test_bounds_use_global_rows_and_leave_other_variables_unchanged():
    values=np.array([[-5.,130.],[-5.,130.],[10.,20.]])
    result=bound_analysis_rows(values,np.array([4,0,5]),{"h":(0,3),"forcing":(3,6)},{"forcing":[0,120]})
    np.testing.assert_array_equal(result,[[0,120],[-5,130],[10,20]])
    np.testing.assert_array_equal(values,[[-5,130],[-5,130],[10,20]])
    assert bound_analysis_rows(values,np.array([4,0,5]),{},None) is values


def test_bounds_apply_after_inflation_without_touching_other_blocks():
    forecast=np.zeros((4,2));analysis=np.array([[1.,2.],[3.,4.],[-1.,119.],[0.,130.]])
    result=apply_analysis_controls_local(analysis,forecast,np.arange(4),dict(
        vec_inputs=["h","forcing"],nd=4,num_state_vars=2,inflation_factor=2.,
        analysis_bounds={"forcing":[0,120]}))
    assert result[2:].min()>=0 and result[2:].max()<=120
    np.testing.assert_allclose(result[:2],[[.5,2.5],[2.5,4.5]])


@pytest.mark.parametrize("bounds", [{"missing":[0,1]},{"h":[2,1]},{"h":[0,np.inf]},{"h":[1]}])
def test_invalid_bounds_are_rejected(bounds):
    with pytest.raises(ValueError,match="analysis_bounds"):
        bound_analysis_rows(np.zeros((1,2)),np.array([0]),{"h":(0,1)},bounds)
