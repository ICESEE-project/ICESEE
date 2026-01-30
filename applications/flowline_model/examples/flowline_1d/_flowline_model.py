# =============================================================================
# Implicit 1D flowline model (Jax version)
# @author: Brian Kyanjo
# @date: 2024-09-24
# @description: This script includes the flowline model using JAX, 
#               - Jacobian computation using JAX
#               - Implicit solver using JAX
# =============================================================================

# Import libraries ====
import jax
import numpy as np
from jax import jacfwd
import jax.numpy as jnp
import matplotlib.pyplot as plt
from scipy.optimize import root

# --- Utility imports ---
from ICESEE.config._utility_imports import icesee_get_index

jax.config.update("jax_enable_x64", True) # Set the precision in JAX to use float64

# Bed topography function --------------------------------------------------------------
def bed(x, **kwargs):
    """
    Bed topography function, which computes the bed shape based on input x and model parameters.
    
    Parameters:
    x (jax.numpy array): Input spatial grid points.
    
    Returns:
    jax.numpy array: The bed topography values at each x location.
    """
    import jax.numpy as jnp
    
    # Ensure parameters are floats
    sillamp    = kwargs['sillamp']
    sillsmooth = kwargs['sillsmooth']
    xsill      = kwargs['xsill']

    # Compute the bed topography
    b = sillamp * (-2 * jnp.arccos((1 - sillsmooth) * jnp.sin(jnp.pi * x / (2 * xsill))) / jnp.pi - 1)
    return b

# Implicit flowline model function (Jax version) --------------------------------------------------------------
def flowline(varin, varin_old, **kwargs):
    # Unpack grid
    NX          = kwargs["NX"]
    N1          = kwargs["N1"]
    dt          = kwargs["dt"] / kwargs["tscale"]
    ds          = kwargs["grid"]["dsigma"]
    sigma       = kwargs["grid"]["sigma"]
    sigma_elem  = kwargs["grid"]["sigma_elem"]

    # Unpack parameters
    tcurrent    = kwargs["tcurrent"]
    xscale      = kwargs["xscale"]
    hscale      = kwargs["hscale"]
    lambd       = kwargs["lambda"]
    m           = kwargs["m"]
    n           = kwargs["n"]
    a           = kwargs["accum"] / kwargs["ascale"]
    eps         = kwargs["eps"]
    transient   = kwargs["transient"]

    # put a guard on mdot, it could be a scalar or an array
    if isinstance(kwargs["facemelt"], (int, float)):
        mdot = kwargs["facemelt"] / kwargs["uscale"]
    else:
        mdot   = kwargs["facemelt"][tcurrent]/kwargs["uscale"]
    # Unpack variables
    h = varin[0:NX]
    u = varin[NX:2*NX]
    xg = varin[2*NX]

    h_old = varin_old[0:NX]
    xg_old = varin_old[2*NX]


    # Calculate bed 
    hf  = -bed(xg * xscale, **kwargs) / (hscale * (1 - lambd))
    hfm = -bed(xg * sigma_elem[-1] * xscale, **kwargs) / (hscale * (1 - lambd))
    b   = -bed(xg * sigma * xscale, **kwargs) / hscale

    # Initialize the residual vector
    F = jnp.zeros(2 * NX + 1, dtype=jnp.float64)

    # Calculate thickness functions        
    F = F.at[0].set(transient * (h[0] - h_old[0]) / dt + (2 * h[0] * u[0]) / (ds[0] * xg)  - a)
    
    F = F.at[1].set(
        transient * (h[1] - h_old[1]) / dt
        - transient * sigma_elem[1] * (xg - xg_old) * (h[2] - h[0]) / (2 * dt * ds[1] * xg)
        + (h[1] * (u[1] + u[0])) / (2 * xg * ds[1]) - a
    )

    F = F.at[2:NX-1].set(
        transient * (h[2:NX-1] - h_old[2:NX-1]) / dt
        - transient * sigma_elem[2:NX-1] * (xg - xg_old) * (h[3:NX] - h[1:NX-2]) / (2 * dt * ds[2:NX-1] * xg)
        + (h[2:NX-1] * (u[2:NX-1] + u[1:NX-2]) - h[1:NX-2] * (u[1:NX-2] + u[0:NX-3])) / (2 * xg * ds[2:NX-1]) - a
    )

    F = F.at[N1-1].set(
        (1 + 0.5 * (1 + (ds[N1-1] / ds[N1-2]))) * h[N1-1]
        - 0.5 * (1 + (ds[N1-1] / ds[N1-2])) * h[N1-2]
        - h[N1]
    )

    F = F.at[NX-1].set(
    transient * (h[NX-1] - h_old[NX-1]) / dt
    - transient * sigma[NX-1] * (xg - xg_old) * (h[NX-1] - h[NX-2]) / (dt * ds[NX-2] * xg)
    + (h[NX-1] * (u[NX-1] + mdot * hf / h[NX-1] + u[NX-2]) - h[NX-2] * (u[NX-2] + u[NX-3])) / (2 * xg * ds[NX-2])
    - a
    )
    
    # Calculate velocity functions
    F = F.at[NX].set(
        ((4 * eps / (xg * ds[0]) ** ((1 / n) + 1)) * (h[1] * (u[1] - u[0]) * abs(u[1] - u[0]) ** ((1 / n) - 1)
            - h[0] * (2 * u[0]) * abs(2 * u[0]) ** ((1 / n) - 1)))
        - u[0] * abs(u[0]) ** (m - 1)
        - 0.5 * (h[0] + h[1]) * (h[1] - b[1] - h[0] + b[0]) / (xg * ds[0])
    )

    F = F.at[NX+1:2*NX-1].set(
        (4 * eps / (xg * ds[1:NX-1]) ** ((1 / n) + 1))
        * (h[2:NX] * (u[2:NX] - u[1:NX-1]) * abs(u[2:NX] - u[1:NX-1]) ** ((1 / n) - 1)
           - h[1:NX-1] * (u[1:NX-1] - u[0:NX-2]) * abs(u[1:NX-1] - u[0:NX-2]) ** ((1 / n) - 1))
        - u[1:NX-1] * abs(u[1:NX-1]) ** (m - 1)
        - 0.5 * (h[1:NX-1] + h[2:NX]) * (h[2:NX] - b[2:NX] - h[1:NX-1] + b[1:NX-1]) / (xg * ds[1:NX-1])
    )

    F = F.at[NX+N1-1].set((u[N1] - u[N1-1]) / ds[N1-1] - (u[N1-1] - u[N1-2]) / ds[N1-2])
    F = F.at[2*NX-1].set(
        (1 / (xg * ds[NX-2]) ** (1 / n)) * (abs(u[NX-1] - u[NX-2]) ** ((1 / n) - 1)) * (u[NX-1] - u[NX-2])
        - lambd * hf / (8 * eps)
    )

    # Calculate grounding line functions
    F = F.at[2*NX].set(3 * h[NX-1] - h[NX-2] - 2 * hf)

    return F

# Calculate the Jacobian of the flowline model function --------------------------------------------------------------
def Jac_calc(huxg_old, **kwargs):
    """
    Use automatic differentiation to calculate Jacobian for nonlinear solver.
    """

    def f(varin):
        # Call the flowline function with current arguments
        return flowline(varin, huxg_old, **kwargs)

    # Create a function that calculates the Jacobian using JAX
    def Jf(varin):
        # Jacobian of f with respect to varin
        return jax.jacfwd(f)(varin)

    return Jf

# initialize_model function --------------------------------------------------------------
def initialize_model(**kwargs):
    """
    Initialize the flowline model by solving for the initial state.
    """
    xg = 300e3/kwargs['xscale']
    hf = (-bed(xg * kwargs['xscale'], **kwargs)/kwargs["hscale"])/(1 -kwargs["lambda"])
    h = 1 - (1-hf)*kwargs['grid']['sigma']
    u = 1*(kwargs['grid']["sigma_elem"]**(1/3)) + 1e-3
    huxg_old = np.concatenate((h, u, [xg]))

    # print(f"Initializing the flowline model h: {h}, u: {u}, xg: {xg}, hf: {hf}")

    Jf = Jac_calc(huxg_old, **kwargs)
    flowline_func = lambda varin: flowline(varin, huxg_old, **kwargs)
    solve_result = root(
        flowline_func,
        huxg_old,
        jac=Jf,
        method='hybr',  # Hybr is a commonly used solver like nlsolve
        options={'maxiter': 100}
    )
    huxg_out0 = solve_result.x
    return huxg_out0

# Function that runs the flowline model function --------------------------------------------------------------
def run_model(ensemble, **kwargs):
    """
    Run the implicit flowline model over nt time steps.
    """

    varin = ensemble.copy()
    last_entry = varin[-1]
    
    # Jacobian calculation
    Jf = Jac_calc(ensemble, **kwargs)

    # Solve the system of nonlinear equations
    solve_result = root(
        lambda varin: flowline(varin, ensemble, **kwargs), 
        ensemble, 
        jac=Jf, 
        method='hybr',  # Hybr is a commonly used solver like nlsolve
        options={'maxiter': 100}
    )
    
    # Update the old solution
    ensemble = solve_result.x
    vecs, indx_map, dim_per_proc = icesee_get_index(**kwargs)
    # update with last entry
    ensemble[-1] = last_entry
    updated_state = {}
    for key in kwargs['vec_inputs']:
        updated_state[key] = ensemble[indx_map[key]]
    return updated_state










