

def bed(self, x):
    """
    Bed topography function, which computes the bed shape based on input x and model parameters.
    
    Parameters:
    x (jax.numpy array): Input spatial grid points.
    
    Returns:
    jax.numpy array: The bed topography values at each x location.
    """
    import jax.numpy as jnp
    
    # Ensure parameters are floats
    params     = self.params
    sillamp    = float(params['sillamp'])
    sillsmooth = float(params['sillsmooth'])
    xsill      = float(params['xsill'])

    # Compute the bed topography
    b = sillamp * (-2 * jnp.arccos((1 - sillsmooth) * jnp.sin(jnp.pi * x / (2 * xsill))) / jnp.pi - 1)
    return b