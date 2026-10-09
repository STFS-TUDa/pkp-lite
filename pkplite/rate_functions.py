"""Rate functions for the detailed (Cantera) and simplified pyrolysis models.

Any rate function passed to `reactor` must accept (T, Y, ...) and return
dY/dt, and expose `.n_coefficients` / `.coefficient_names` so the fitting
code knows how to build the coefficients dict it returns.

Unit note: the Arrhenius term below uses `exp(-Ea / (T * 8314))`, i.e. R is
taken as 8314 instead of 8.314 J/(mol*K). This means Ea as stored in
coefficients.yml is 1000x the SI value in J/mol (divide by 1e6 to get kJ/mol).
"""

import cantera as ct
import numpy as np

def rate_function_detailed(T, Y, mechanism, pressure=ct.one_atm):
    """Species production rates (mass basis) from a Cantera solid-phase mechanism."""

    mechanism.TPY = T, pressure, Y

    return ( mechanism.net_production_rates *
                mechanism.molecular_weights / mechanism.density )

def C2SM(T, Y, coefficients):
    """Competing 2-step model: two parallel first-order reactions, each with
    its own ultimate volatile yield (alpha), pre-exponential factor (A) and
    activation energy (Ea), competing for the same raw fuel.

    Y = [volatiles, char, raw].
    """

    alpha1, A1, Ea1, alpha2, A2, Ea2 = coefficients

    rate1 = ( A1 * np.exp(-Ea1/(T * 8314)) )
    rate2 = ( A2 * np.exp(-Ea2/(T * 8314)) )

    devol = Y[-1] * ( alpha1 * rate1 + alpha2 * rate2 )
    char = Y[-1] * ( (1 - alpha1) * rate1 + (1 - alpha2) * rate2 )

    return np.array([devol,char,-devol-char])

C2SM.n_coefficients = 6
C2SM.coefficient_names = ["alpha1", "A1", "Ea1", "alpha2", "A2", "Ea2"]

def SFOR(T, Y, coefficients):
    """Single First-Order Reaction model: dV/dt = A * exp(-Ea/RT) * alpha * raw.

    Y = [volatiles, char, raw]. alpha is the ultimate volatile yield (0-1).
    """
    alpha, A, Ea = coefficients

    rate = A * np.exp(-Ea / (T * 8314))

    devol = Y[-1] * alpha * rate
    char  = Y[-1] * (1 - alpha) * rate

    return np.array([devol, char, -devol - char])

SFOR.n_coefficients = 3
SFOR.coefficient_names = ["alpha", "A", "Ea"]
