# ==============================================================================
# @des: This file contains the lorenz model function and the 
#       4th order Runge-Kutta integrator to run the model.
# @date: 2025-01-18
# @author: Brian Kyanjo
# ==============================================================================

import numpy as np

# --- Utility imports ---
from ICESEE.config._utility_imports import icesee_get_index

# --- model initialization ---
def initialize_model(**kwargs):
    """des: initialize the icepack model"""
    pass

# --- run function for the lorenz96 model ---
def Lorenz96(state, **kwargs):
    """des: Lorenz96 model function
        inputs: state - current state of the model
                **kwargs - additional arguments for the model
       outputs: f - the derivative of the state vector
    """
    # Unpack the arguments
    sigma = kwargs.get('sigma', None)
    beta  = kwargs.get('beta', None)
    rho   = kwargs.get('rho', None)

    x,y,z = state # Unpack the state vector
    f = np.zeros(3) # Create an empty vector to store the derivatives
    f[0] = sigma*(y-x)  
    f[1] = x*(rho-z)-y
    f[2] = x*y - beta*z
    return f

# --- 4th order Runge-Kutta integrator --- 
def RK4(rhs, state, **kwargs):
    """des: 4th order Runge-Kutta integrator
        inputs: rhs - function that computes the right-hand side of the ODE
                state - current state of the model
                dt - time step
                *args - additional arguments for the model
        outputs: state - updated state of the model after one time step
    """
    dt = kwargs.get('dt', None)
    k1 = rhs(state, **kwargs)
    k2 = rhs(state + 0.5*dt*k1, **kwargs)
    k3 = rhs(state + 0.5*dt*k2, **kwargs)
    k4 = rhs(state + dt*k3, **kwargs)
    return state + dt/6*(k1 + 2*k2 + 2*k3 + k4)

# --- Run similation for the Lorenz96 model ---
def run_model(ensemble, **kwargs):
    """des: Lorenz96 model function
        inputs: ensemble - current ensemble state of the model
                **kwargs - additional arguments for the model
        outputs: model run
    """

    # call the icesee_get_index function to get the indices of the state variables
    vecs, indx_map, dim_per_proc = icesee_get_index(ensemble, **kwargs)
     
    # Call the RK4 function to push the state forward in time
    state = RK4(Lorenz96, ensemble, **kwargs)

    # return the updated state
    updated_state = {'x' : state[indx_map['x']],
                     'y' : state[indx_map['y']],
                     'z' : state[indx_map['z']]}
    return updated_state