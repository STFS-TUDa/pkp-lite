"""Fit a simplified model's coefficients to detailed-model results using
scipy's TRF least-squares solver (`least_squares_fit`), or the same in
log10-space for parameters spanning many orders of magnitude
(`least_squares_fit_log`). Generally faster than the evolutionary approach
in `evolution_fit.py` for well-behaved problems, at the cost of being more
sensitive to the initial guess / local minima.
"""

import numpy as np
import pandas as pd
from scipy.optimize import least_squares


def _residuals(coefficients, model_reactor, operating_points, results):
    """Concatenated (detailed - model) `volatiles` residuals across all runs,
    normalized the same way as Vascellari et al. (2013)'s objective function:
    each run's residuals are divided by its own (max-min) yield range and by
    sqrt(number of timesteps * number of runs), so that scipy's internal
    sum-of-squares cost gives every operating point (heating rate) equal
    weight regardless of its absolute yield scale or sample count."""
    model_reactor.rate_function_arguments = {"coefficients": coefficients}
    n_runs = len(operating_points)
    res_list = []
    for op, res in zip(operating_points, results):
        model_res = model_reactor.run(op)
        volatiles = res['volatiles'].to_numpy()
        y_range = volatiles.max() - volatiles.min()
        residual = volatiles - model_res[:, 1]
        res_list.append(residual / (y_range * np.sqrt(len(volatiles) * n_runs)))
    return np.concatenate(res_list)


def least_squares_fit(fit_settings, model_reactor, operating_points, results):
    """Run scipy TRF least-squares and return (fitted model results per run,
    coefficients dict). fit_settings: the `fit` block from input.yml (bounds,
    optional parameters_init, verbose)."""
    bounds = {
        key: np.array([float(i) for i in fit_settings['bounds'][key]])
        for key in fit_settings['bounds']
    }

    if 'parameters_init' in fit_settings:
        x0 = np.array(fit_settings['parameters_init'], dtype=float)
    else:
        x0 = (bounds['min'] + bounds['max']) / 2

    call_count = [0]

    def _residuals_counted(coefficients, model_reactor, operating_points, results):
        call_count[0] += 1
        r = _residuals(coefficients, model_reactor, operating_points, results)
        cost = 0.5 * (r ** 2).sum()
        print(f"\r  Fitting... iter {call_count[0]:4d}  cost={cost:.4e}", end="", flush=True)
        return r

    result = least_squares(
        _residuals_counted,
        x0,
        bounds=(bounds['min'], bounds['max']),
        args=(model_reactor, operating_points, results),
        method='trf',
        verbose=fit_settings.get('verbose', 0),
    )
    print()

    coefficients = result.x
    model_reactor.rate_function_arguments = {"coefficients": coefficients}
    results_model = [model_reactor.run(op) for op in operating_points]

    coeff_names = model_reactor.rate_function.coefficient_names
    return (
        [pd.DataFrame(res, columns=['t', 'volatiles', 'char', 'raw', 'T']) for res in results_model],
        {key: val for key, val in zip(coeff_names, coefficients)} | {'cost': result.cost},
    )


def least_squares_fit_log(fit_settings, model_reactor, operating_points, results):
    """Like least_squares_fit, but fits A and Ea parameters in log10-space.

    Parameters with min > 0 and max/min > 100 are automatically log-transformed,
    which makes TRF converge much faster when bounds span many orders of magnitude.
    """
    bounds_raw = {
        key: np.array([float(i) for i in fit_settings['bounds'][key]])
        for key in fit_settings['bounds']
    }

    # Auto-detect which parameters benefit from log-transform
    log_mask = (bounds_raw['min'] > 0) & (bounds_raw['max'] / bounds_raw['min'] > 100)

    lb = np.where(log_mask, np.log10(bounds_raw['min']), bounds_raw['min'])
    ub = np.where(log_mask, np.log10(bounds_raw['max']), bounds_raw['max'])

    # Default init: geometric center in log-space, arithmetic center in linear space
    x0_default = np.where(log_mask,
                          10 ** ((np.log10(np.maximum(bounds_raw['min'], 1e-300)) + np.log10(bounds_raw['max'])) / 2),
                          (bounds_raw['min'] + bounds_raw['max']) / 2)
    if 'parameters_init' in fit_settings:
        x0_raw = np.array(fit_settings['parameters_init'], dtype=float)
        # Clamp to strictly inside bounds so we never start on a boundary
        x0_raw = np.clip(x0_raw, bounds_raw['min'] * 1.001 + 1e-300, bounds_raw['max'] * 0.999)
    else:
        x0_raw = x0_default
    x0 = np.where(log_mask, np.log10(x0_raw), x0_raw)

    def to_physical(x):
        return np.where(log_mask, 10 ** x, x)

    call_count = [0]

    def _residuals_log(x_log, model_reactor, operating_points, results):
        call_count[0] += 1
        r = _residuals(to_physical(x_log), model_reactor, operating_points, results)
        cost = 0.5 * (r ** 2).sum()
        print(f"\r  Fitting (log)... iter {call_count[0]:4d}  cost={cost:.4e}", end="", flush=True)
        return r

    result = least_squares(
        _residuals_log,
        x0,
        bounds=(lb, ub),
        args=(model_reactor, operating_points, results),
        method='trf',
        verbose=fit_settings.get('verbose', 0),
    )
    print()

    coefficients = to_physical(result.x)
    model_reactor.rate_function_arguments = {"coefficients": coefficients}
    results_model = [model_reactor.run(op) for op in operating_points]

    coeff_names = model_reactor.rate_function.coefficient_names
    return (
        [pd.DataFrame(res, columns=['t', 'volatiles', 'char', 'raw', 'T']) for res in results_model],
        {key: val for key, val in zip(coeff_names, coefficients)} | {'cost': result.cost},
    )
