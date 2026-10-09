"""Fit a simplified model's coefficients to detailed-model results using a
DEAP mu+lambda evolutionary algorithm.

Coefficients are normalized to [0, 1] internally (each individual is a list
of floats in that range) and mapped back to physical units via `bounds`
before being passed to the rate function; this keeps mutation/crossover
scale-independent across parameters that span very different orders of
magnitude (e.g. alpha in [0,1] vs. A in [1e0, 1e12]).
"""

import numpy as np
import random
import pandas as pd
import multiprocessing as mp

from deap import base
from deap import creator
from deap import tools
from deap import algorithms

# create only once at import
if not hasattr(creator, "FitnessMin"):
    creator.create("FitnessMin", base.Fitness, weights=(-1.0,))
if not hasattr(creator, "Individual"):
    creator.create("Individual", list, fitness=creator.FitnessMin)

def evaluate_error(bounds, model_reactor, operating_points, results, final_error_weight, individual):
    """DEAP fitness function, following the objective function of Vascellari
    et al. (2013): for each operating point (heating rate), the mean squared
    error of the model's `volatiles` trace against the detailed-model
    `results` is normalized by that run's own (max-min) yield range, then
    averaged with equal weight across operating points - so a run with a
    larger absolute volatile-yield swing doesn't dominate the fit just
    because of its scale. The optional final-time error penalty (weighted by
    `final_error_weight`) is normalized the same way."""
    coefficients = np.array(individual) * (bounds['max'] - bounds['min']) + bounds['min']
    model_reactor.rate_function_arguments = {"coefficients": coefficients}
    error = 0

    for op, res in zip(operating_points, results):
        model_res = model_reactor.run(op)

        volatiles = res['volatiles'].to_numpy()
        y_range = volatiles.max() - volatiles.min()

        mse = np.mean((volatiles - model_res[:,1])**2) / y_range**2
        final_err = final_error_weight * np.abs(volatiles[-1] - model_res[-1,1]) / y_range
        error += (mse + final_err) / len(operating_points)
    return error,

def check_bounds(min, max):
    """
    Check if the individual exists in the given boundaries.

    Decorator which fix the limit of the function arguments to min and
    max values.

    Parameters
    ----------
    min: minimum value of the parameter
    max: minimum value of the parameter

    """
    def decorator(func):
        def wrappper(*args, **kargs):
            offspring = func(*args, **kargs)
            for child in offspring:
                for i in range(len(child)):
                    if child[i] > max:
                        child[i] = max
                    elif child[i] < min:
                        child[i] = min
            return offspring
        return wrappper
    return decorator

def coefficients_to_individual(coeffs, bounds):
    """Normalize physical-unit coefficients (e.g. A1, Ea1, ...) to [0,1] so
    they can seed the population as a DEAP Individual."""
    normalized = (coeffs - bounds['min']) / (bounds['max'] - bounds['min'])
    return creator.Individual(normalized.tolist())

def _best_of_generation_stats(bounds, coeff_names):
    """DEAP Statistics object that reports, once per generation, the best
    individual's fitness and its decoded (physical-unit) coefficients -
    instead of printing every single individual evaluated."""
    def best(individuals):
        return min(individuals, key=lambda ind: ind.fitness.values[0])

    stats = tools.Statistics(key=lambda ind: ind)
    stats.register("fitness", lambda pop: best(pop).fitness.values[0])
    for i, name in enumerate(coeff_names):
        stats.register(name, lambda pop, i=i:
            best(pop)[i] * (bounds['max'][i] - bounds['min'][i]) + bounds['min'][i])
    return stats

def evolution_fit(fit_settings,model_reactor,operating_points,results):
    """Run the EA and return (fitted model results per run, coefficients dict).

    fit_settings: the `fit` block from input.yml (bounds, npop, ngen, mu,
    lambda_, cxpb, mutpb, final_error_weight, n_jobs, optional
    parameters_init used to warm-start the first individual).
    """

    # creator.create("FitnessMin", base.Fitness, weights=(-1.0,))
    # creator.create("Individual", list, fitness=creator.FitnessMin)

    toolbox = base.Toolbox()
    toolbox.register("attr_float", random.random)
    # individual uses n chromosomes (as many model parameters)
    # and attr_float
    toolbox.register("individual", tools.initRepeat,
                        creator.Individual, toolbox.attr_float,
                        n=model_reactor.rate_function.n_coefficients)
    toolbox.register("population", tools.initRepeat, list, toolbox.individual)
    bounds = {key: np.array([float(i) for i in fit_settings['bounds'][key]]) for key in fit_settings['bounds']}
    toolbox.register("evaluate", evaluate_error, bounds, model_reactor, operating_points, results, fit_settings.get('final_error_weight',0))

    # Define the algorithms 
    toolbox.register('mate', tools.cxBlend, alpha=0.25)
    toolbox.register('mutate', tools.mutGaussian, mu=0, sigma=1,
                        indpb=0.2)
    toolbox.register('select', tools.selTournament, tournsize=3)
    toolbox.decorate("mate", check_bounds(0, 1))
    toolbox.decorate("mutate", check_bounds(0, 1))

    npop = fit_settings.get("npop",20)
    pop = toolbox.population(n=npop)

    if 'parameters_init' in fit_settings:
        init_coeffs = np.array(fit_settings['parameters_init'], dtype=float)
        pop[0] = coefficients_to_individual(init_coeffs, bounds)

    hof = tools.HallOfFame(1)
    number_processors = fit_settings.get("n_jobs", mp.cpu_count())
    coeff_names = model_reactor.rate_function.coefficient_names
    stats = _best_of_generation_stats(bounds, coeff_names)

    with mp.Pool(processes=number_processors) as pool:
        toolbox.register("map", pool.map)

        pop, log = algorithms.eaMuPlusLambda(
            pop,
            toolbox,
            mu=fit_settings.get("mu", npop),
            lambda_=fit_settings.get("lambda_", npop),
            cxpb=fit_settings.get("cxpb", 0.5),
            mutpb=fit_settings.get("mutpb", 0.5),
            ngen=fit_settings.get("ngen", npop),
            halloffame=hof,
            stats=stats,
            verbose=True,
        )

    fitnesses = np.array([p.fitness.values for p in pop])
    best = pop[fitnesses.argmin()]

    coefficients = best * (bounds['max'] - bounds['min']) + bounds['min']

    print("Best:", dict(zip(coeff_names, coefficients)), " error:", np.min(fitnesses))
    model_reactor.rate_function_arguments = {"coefficients": coefficients}

    results_model = [model_reactor.run(op) for op in operating_points]

    return [pd.DataFrame(res,columns=['t','volatiles','char','raw','T']) for res in results_model], {key:val for key,val in zip(coeff_names,coefficients)}