# run_jax.py

# imports
import jax
import jax.numpy as jnp
import numpy as np
import matplotlib.pyplot as plt
import optax
from utils import constants

# python -m jax.run_jax

@jax.jit
def planck_function(nu, T):
    """
    B_nu(T): the blackbody function
    """
    # clipping T to avoid any division by zero or negative temps
    T = jnp.maximum(T, 1e-5)
    x = constants.h * nu / (constants.k_b * T)
    # using jnp.expm1 for numerical stability
    return (2 * constants.h * nu**3 / constants.c**2) / jnp.expm1(x)

@jax.jit
def modified_blackbody(params, nu, distance_cm, kappa_0, nu_0, beta):
    """
    Generates model flux S_nu given parameters
    params: [log10_mass, temperature]
    """
    log_M, T = params
    M_dust = jnp.power(10, log_M) * 1.989e33 # converting log(M_sun) to grams

    # opacity power law
    # kappa = kappa_0 * (nu / nu_0)^beta
    kappa = kappa_0 * (nu / nu_0)**beta

    # RTE (optically thin)
    flux = (M_dust * kappa * planck_function(nu, T)) / distance_cm**2
    return flux

@jax.jit
def loss_function(params, nu_obs, flux_obs, distance_cm, kappa_0, nu_0, beta):
    flux_model = modified_blackbody(params, nu_obs, distance_cm, kappa_0, nu_0, beta)
    # log-spaced error
    diff = jnp.log10(flux_model) - jnp.log10(flux_obs)
    return jnp.mean(diff**2)

