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
from scipy.stats import norm, multivariate_normal

jax.config.update("jax_enable_x64", True) # Set the precision in JAX to use float64

# Implicit flowline model function (Jax version) --------------------------------------------------------------
def flowline(varin, varin_old, params, grid, bedfun):
    # Unpack grid
    NX          = params["NX"]
    N1          = params["N1"]
    dt          = params["dt"] / params["tscale"]
    ds          = grid["dsigma"]
    sigma       = grid["sigma"]
    sigma_elem  = grid["sigma_elem"]

    # Unpack parameters
    tcurrent    = params["tcurrent"]
    xscale      = params["xscale"]
    hscale      = params["hscale"]
    lambd       = params["lambda"]
    m           = params["m"]
    n           = params["n"]
    a           = params["accum"] / params["ascale"]
    eps         = params["eps"]
    transient   = params["transient"]

    # put a guard on mdot, it could be a scalar or an array
    if isinstance(params["facemelt"], (int, float)):
        mdot = params["facemelt"] / params["uscale"]
    else:
        mdot   = params["facemelt"][tcurrent]/params["uscale"]

    # Unpack variables
    h = varin[0:NX]
    u = varin[NX:2*NX]
    xg = varin[2*NX]

    h_old = varin_old[0:NX]
    xg_old = varin_old[2*NX]


    # Calculate bed 
    hf  = -bedfun(xg * xscale) / (hscale * (1 - lambd))
    hfm = -bedfun(xg * sigma_elem[-1] * xscale) / (hscale * (1 - lambd))
    b   = -bedfun(xg * sigma * xscale) / hscale

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
def Jac_calc(huxg_old, params, grid, bedfun, flowlinefun):
    """
    Use automatic differentiation to calculate Jacobian for nonlinear solver.
    """

    def f(varin):
        # Initialize F as an array of zeros with size 2*NX + 1
        # F = jnp.zeros(2 * params["NX"] + 1, dtype=jnp.float64)
        # Call the flowline function with current arguments
        return flowlinefun(varin, huxg_old, params, grid, bedfun)
        # print(F)
        # return F

    # Create a function that calculates the Jacobian using JAX
    def Jf(varin):
        # Jacobian of f with respect to varin
        return jax.jacfwd(f)(varin)

    return Jf

# Function that runs the flowline model function --------------------------------------------------------------
def flowline_run(varin, params, grid, bedfun, flowlinefun):
    nt = params["NT"]
    huxg_old = varin
    huxg_all = np.zeros((huxg_old.shape[0], nt))

    for i in range(nt):
        if not params["assim"]:
            params["tcurrent"] = i + 1  # Adjusting for 1-based indexing in Julia
        
        # Jacobian calculation
        Jf = Jac_calc(huxg_old, params, grid, bedfun, flowlinefun)
        
        # Solve the system of nonlinear equations
        solve_result = root(
            lambda varin: flowlinefun(varin, huxg_old, params, grid, bedfun), 
            huxg_old, 
            jac=Jf, 
            method='hybr',  # Hybr is a commonly used solver like nlsolve
            options={'maxiter': 100}
        )
        
        # Update the old solution
        huxg_old = solve_result.x
        
        # Store the result for this time step
        huxg_all[:, i] = huxg_old

        if not params["assim"]:
            print(f"Step {i + 1}\n")
    
    return huxg_all