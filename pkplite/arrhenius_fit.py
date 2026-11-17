"""SFOR parameters by linear regression of the normalized mass-loss rate,

    ln( (-dTG/dt) / (TG - TG_min) ) = ln(A) - Ea / (R T),

inside a temperature window of the detailed run. alpha = 1 - TG_min, with
TG_min the TG at the end of the detailed run, so the end of the T programme
sets up to which temperature the yields apply. This is the single-reaction
fit used for the PATO solid model.

Window settings (fit.window in the input YAML):
    T_min: auto | value in K   auto: T where TG is closest to TG_start
    T_max: auto | value in K   auto: T at the maximum of -dTG/dt
    TG_start: 0.95             used by T_min: auto
"""

import numpy as np
import pandas as pd

R = 8314.0  # J/(kmol K); Ea in pkp-lite is in J/kmol


def fit_window(T, TG, dTG, window):
    """Return (T_min, T_max) of the regression window for one run."""
    TG_start = window.get('TG_start', 0.95)
    T_min = window.get('T_min', 'auto')
    T_max = window.get('T_max', 'auto')
    if T_min == 'auto':
        T_min = T[np.argmin(np.abs(TG - TG_start))]
    if T_max == 'auto':
        T_max = T[np.argmax(dTG)]
    return float(T_min), float(T_max)


def arrhenius_fit(fit_settings, model_reactor, operating_points, results):
    """Fit SFOR [alpha, A, Ea] to the detailed results (all runs pooled)."""
    window = fit_settings.get('window', {})
    x, y, alphas, info = [], [], [], {}
    for i, res in enumerate(results):
        t, T = res['t'].to_numpy(), res['T'].to_numpy()
        TG = res['solid'].to_numpy() / (res['solid'].iloc[0] + res['volatiles'].iloc[0])
        dTG = -np.gradient(TG, t)                     # positive for mass loss
        T_min, T_max = fit_window(T, TG, dTG, window)
        mask = (dTG > 0) & (T >= T_min) & (T <= T_max) & (TG > TG.min())
        x.append(1.0 / T[mask])
        y.append(np.log(dTG[mask] / (TG[mask] - TG.min())))
        alphas.append(1.0 - TG.min())
        info[f'T_min_run{i}'], info[f'T_max_run{i}'] = T_min, T_max
        print(f"  run{i}: fit window {T_min:.1f}-{T_max:.1f} K, {mask.sum()} points, TG_min = {TG.min():.4f}")
    x, y = np.concatenate(x), np.concatenate(y)

    slope, intercept = np.polyfit(x, y, 1)
    coefficients = np.array([np.mean(alphas), np.exp(intercept), -slope * R])

    model_reactor.rate_function_arguments = {"coefficients": coefficients}
    results_model = [model_reactor.run(op) for op in operating_points]

    coeff_names = model_reactor.rate_function.coefficient_names
    return (
        [pd.DataFrame(res, columns=['t', 'volatiles', 'char', 'raw', 'T']) for res in results_model],
        {key: val for key, val in zip(coeff_names, coefficients)}
        | {'R2': np.corrcoef(x, y)[0, 1] ** 2} | info,
    )
