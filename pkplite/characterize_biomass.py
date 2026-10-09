"""Converts elemental and/or biochemical biomass analysis into mass fractions
of CRECK reference species (CELL, XYHW/GMSW/XYGR, LIGC, LIGH, LIGO, TANN,
TGL, ...), for use as the `fuel.reference_species` input to a mechanism."""

import itertools
from collections import OrderedDict

import numpy as np
from scipy.optimize import minimize
from scipy.spatial import Delaunay

# Atomic weights of C, H, O as used in the original CRECK characterization sheet.
_ATOMIC_WEIGHTS = np.array([12.0, 1.0, 16.0])

# Organic reference species and their (C, H, O) atom counts.
_REF_SPECIES = ['CELL', 'HECELL', 'LIGO', 'LIGH', 'LIGC', 'TANN', 'TGL']
_FORMULAS = np.array([
    [6, 10, 5],     # CELL   C6H10O5
    [5, 8, 4],      # HECELL C5H8O4
    [20, 22, 10],   # LIGO   C20H22O10
    [22, 28, 9],    # LIGH   C22H28O9
    [15, 14, 4],    # LIGC   C15H14O4
    [15, 12, 7],    # TANN   C15H12O7
    [57, 100, 7],   # TGL    C57H100O7
])
_MW = _FORMULAS @ _ATOMIC_WEIGHTS

# Protein species share the nitrogen-derived protein mass in fixed proportions.
_PROTEIN_SPECIES = ['PROTC', 'PROTH', 'PROTO']
_PROTEIN_SPLIT = np.array([0.34, 0.33, 0.33])
_PROTEIN_N_MASS_FRACTION = 0.13008701898456757  # N content of the protein mixture

_SPECIES = _REF_SPECIES + _PROTEIN_SPECIES + ['MOIST', 'ASH']

# CRECK hemicellulose species depend on the biomass type.
_HEMICELLULOSE_NAME = {'hardwood': 'XYHW', 'softwood': 'GMSW', 'grass': 'XYGR'}

_SPLITTING_PARAMS = ['alpha', 'beta', 'gamma', 'delta', 'epsilon']
_DEFAULT_SPLITTING_PARAMS = [0.6, 0.8, 0.8, 0.8, 0.8]

# How each biochemical component is distributed over the organic + protein species.
_BIOCHEM_COMPONENTS = ['cellulose', 'hemicellulose', 'lignin', 'extractives', 'proteins']
_REQUIRED_BIOCHEM_COMPONENTS = {'cellulose', 'hemicellulose', 'lignin'}
_BIOCHEM_SPLIT = np.array([
    # CELL HECELL LIGO LIGH LIGC  TANN TGL  PROTC PROTH PROTO
    [1,    0,     0,   0,   0,    0,   0,   0,    0,    0  ],  # cellulose
    [0,    1,     0,   0,   0,    0,   0,   0,    0,    0  ],  # hemicellulose
    [0,    0,   1/3, 1/3, 1/3,    0,   0,   0,    0,    0  ],  # lignin
    [0,    0,     0,   0,   0,  1/2, 1/2,   0,    0,    0  ],  # extractives
    [0,    0,     0,   0,   0,    0,   0, 1/3,  1/3,  1/3  ],  # proteins
])


# -- this function is by Leon Loni Berkel --#
def calculate_reference_species(composition, ash, moist, characterization_type, splitting_params=None, biochem=None, biomass_type=None):
    """Calculate reference species mass fractions (as received basis) based on the characterization type.

        Possible characterization types:
            elemental_analysis: Calculate reference species based on elemental composition - CRECK method. elemental composition and splitting parameters are required.
            biochemical_analysis: Calculate reference species based on biochemical composition - lignin is split equally over LIGO/LIGH/LIGC, extractives over TANN/TGL, proteins over PROTC/PROTH/PROTO.
            elem_and_bio_optimized: Calculate reference species based on both elemental and biochemical composition - CRECK method with optimization of splits, to satisfy biochem.

        Args:
            composition (dict): Dictionary with elemental composition values (C, H, O, optional N) - DAF.
            ash and moist (float): Ash and moisture content - As Received (AR).
            characterization_type (str): One of the types above.
            splitting_params (dict): alpha, beta, gamma, delta, epsilon of the CRECK reference mixtures. Missing ones use defaults (elemental_analysis) or are optimized (elem_and_bio_optimized).
            biochem (dict): cellulose, hemicellulose, lignin and optionally extractives, proteins.
            biomass_type (str): Biomass type - hardwood, softwood or grass. Used for hemicellulose species name change.

        Returns:
            OrderedDict: Mass fractions of CELL, <hemicellulose>, LIGO, LIGH, LIGC, TANN, TGL, PROTC, PROTH, PROTO, MOIST, ASH.
    """
    if characterization_type == 'elemental_analysis':
        splitting_params = splitting_params or {}
        params = [splitting_params.get(name, default) for name, default in zip(_SPLITTING_PARAMS, _DEFAULT_SPLITTING_PARAMS)]
        fractions = _elemental_to_reference_species(composition, ash, moist, params)

    elif characterization_type == 'biochemical_analysis':
        fractions = _biochemical_to_reference_species(biochem, ash, moist)

    elif characterization_type == 'elem_and_bio_optimized':
        params = _optimize_splitting_params(composition, ash, moist, biochem, splitting_params or {})
        fractions = _elemental_to_reference_species(composition, ash, moist, params)

    else:
        raise ValueError(f"Error: Unsupported characterization type {characterization_type}.")

    if biomass_type not in _HEMICELLULOSE_NAME:
        raise ValueError(f"Error: Unsupported biomass type {biomass_type}.")
    names = [_HEMICELLULOSE_NAME[biomass_type] if s == 'HECELL' else s for s in _SPECIES]
    return OrderedDict(zip(names, fractions.tolist()))


def _reference_mixtures(alpha, beta, gamma, delta, epsilon):
    """Molar composition of the three CRECK reference mixtures (rows) in terms of the organic species (columns)."""
    return np.array([
        # CELL  HECELL     LIGO             LIGH          LIGC                   TANN         TGL
        [alpha, 1 - alpha, 0,               0,            0,                     0,           0        ],  # RM1: carbohydrates
        [0,     0,         0,               delta * beta, delta * (1 - beta),    0,           1 - delta],  # RM2: H-rich lignins + TGL
        [0,     0,         epsilon * gamma, 0,            epsilon * (1 - gamma), 1 - epsilon, 0        ],  # RM3: O-rich lignins + tannins
    ])


def _read_elemental(composition):
    """C/H/O mass fractions and N content of the DAF elemental composition. O is taken as 1 - C - H - N
    (the O given in the input is not used)."""
    try:
        carbon, hydrogen = float(composition['C']), float(composition['H'])
        nitrogen = float(composition.get('N', 0))
    except (KeyError, TypeError, ValueError):
        raise ValueError("Error: Elemental composition not available.")
    return np.array([carbon, hydrogen, 1 - carbon - hydrogen - nitrogen]), nitrogen


def _as_received_to_daf(composition, moist, ash):
    """Convert an as received (AR) elemental composition, whose H includes the hydrogen of the moisture,
    to dry ash-free (DAF). The moisture's O drops out since O is taken by difference (see _read_elemental)."""
    try:
        carbon, hydrogen = float(composition['C']), float(composition['H'])
        nitrogen = float(composition.get('N', 0))
    except (KeyError, TypeError, ValueError):
        raise ValueError("Error: Elemental composition not available.")
    h_in_water = 2 * _ATOMIC_WEIGHTS[1] / (2 * _ATOMIC_WEIGHTS[1] + _ATOMIC_WEIGHTS[2])
    daf = 1 - moist - ash
    return {'C': carbon / daf, 'H': (hydrogen - h_in_water * moist) / daf, 'N': nitrogen / daf}


def _elemental_to_reference_species(composition, ash, moist, splitting_params):
    """CRECK method: find the mix of the 3 reference mixtures matching the C/H/O content of the fuel,
    add proteins based on the N content, and return mass fractions in the order of _SPECIES."""
    cho, nitrogen = _read_elemental(composition)

    w_cho = cho / cho.sum()                       # C/H/O mass fractions, N-free
    w_n = nitrogen / (cho.sum() + nitrogen)

    rm_species = _reference_mixtures(*splitting_params)           # (3 RMs, 7 species) mole fractions
    rm_atoms = rm_species @ _FORMULAS                             # (3 RMs, C/H/O) atoms per mole
    rm_mw = rm_atoms @ _ATOMIC_WEIGHTS
    rm_w = rm_atoms * _ATOMIC_WEIGHTS / rm_mw[:, None]            # element mass fractions of each RM

    # RM mass fractions reproducing the fuel's C/H/O (the O balance enforces that they sum to 1)
    try:
        rm_mass = np.linalg.solve(rm_w.T, w_cho)
    except np.linalg.LinAlgError:
        raise ValueError("Error: Reference mixtures are degenerate for these splitting parameters.")
    rm_moles = rm_mass / rm_mw
    rm_moles /= rm_moles.sum()

    x_daf = rm_moles @ rm_species                 # species mole fractions, DAF
    w_daf = x_daf * _MW
    w_daf /= w_daf.sum()                          # species mass fractions, DAF

    daf = 1 - moist - ash
    w_prot = w_n / _PROTEIN_N_MASS_FRACTION * daf if nitrogen > 0 else 0.0
    fractions = np.concatenate([
        w_daf * daf * (1 - w_prot),
        w_prot * _PROTEIN_SPLIT,
        [moist, ash],
    ])
    fractions /= fractions.sum()

    for name, value in zip(_SPECIES, fractions):
        if value < 0 or value > 1:
            # Possible improvement, characterization zone can be plotted
            raise ValueError(f"Error: Value for {name} is out of bounds: {value}. Your elemental composition and splitting parameters probably take you outside of the valid characterization zone.")
    return fractions


def _van_krevelen(cho_moles):
    """(O/C, H/C) molar ratios of C/H/O mole amounts (last axis)."""
    cho_moles = np.asarray(cho_moles, dtype=float)
    return np.stack([cho_moles[..., 2] / cho_moles[..., 0], cho_moles[..., 1] / cho_moles[..., 0]], axis=-1)


def _read_biochem(biochem):
    """Biochemical composition as an array in the order of _BIOCHEM_COMPONENTS (optional components default to 0)."""
    values = []
    for component in _BIOCHEM_COMPONENTS:
        try:
            values.append(float(biochem[component]))
        except (KeyError, TypeError, ValueError):
            if component in _REQUIRED_BIOCHEM_COMPONENTS:
                raise ValueError("Error: Biochemical composition not available.")
            values.append(0.0)
    return np.array(values)


def _biochemical_to_reference_species(biochem, ash, moist):
    """Distribute the normalized biochemical composition over the reference species (order of _SPECIES)."""
    components = _read_biochem(biochem)
    total = components.sum()
    if total > 1.01 or total < 0.99:
        print(f"Characterization warning: Total biochemical composition is not 100% ({total:.2f})")
    species = components / total @ _BIOCHEM_SPLIT * (1 - moist - ash)
    return np.concatenate([species, [moist, ash]])


def _optimize_splitting_params(composition, ash, moist, biochem, fixed):
    """Splitting parameters for which the CRECK elemental method best reproduces the biochemical composition.
    Parameters given in `fixed` are kept at their value, the others are optimized."""
    target = _read_biochem(biochem)
    target = target # / target.sum() * (1 - moist - ash)       # DAF -> as received, like the fitted fractions
    species_to_component = (_BIOCHEM_SPLIT > 0).T            # (10 species, 5 components)

    def objective(params):
        try:
            fractions = _elemental_to_reference_species(composition, ash, moist, params)
        except ValueError:
            return np.inf
        components = fractions[:len(species_to_component)] @ species_to_component
        return np.sum((target - components) ** 2)

    # Same as the original characterization: gradient-based search from the default guess
    bounds = [(fixed[name], fixed[name]) if name in fixed else (0, 1) for name in _SPLITTING_PARAMS]
    initial_guess = [fixed.get(name, default) for name, default in zip(_SPLITTING_PARAMS, _DEFAULT_SPLITTING_PARAMS)]
    result = minimize(objective, initial_guess, bounds=bounds)
    if not np.isfinite(result.fun):
        result = _feasible_search(objective, composition, fixed, bounds)
    print(f"Optimized splitting parameters: {[f'{x:.4f}' for x in result.x]}")
    return result.x


def _feasible_search(objective, composition, fixed, bounds):
    """Fallback when the default guess gives negative species (objective is inf there, so no gradient):
    check that the sample can be characterized at all, then start a gradient-free search from the
    best feasible points of a coarse grid."""
    # Van Krevelen check: the sample must lie inside the region spanned by the reference species
    cho_moles = _read_elemental(composition)[0] / _ATOMIC_WEIGHTS
    if Delaunay(_van_krevelen(_FORMULAS)).find_simplex(_van_krevelen(cho_moles)) < 0:
        raise ValueError("Error: Elemental composition lies outside the Van Krevelen region of the reference species, "
                         "it cannot be characterized with any splitting parameters.")

    grid = [[fixed[name]] if name in fixed else np.linspace(0, 1, 6) for name in _SPLITTING_PARAMS]
    starts = [start for start in sorted(itertools.product(*grid), key=objective)[:5] if np.isfinite(objective(start))]
    if not starts:
        raise ValueError("Error: No feasible splitting parameters found on the coarse grid.")

    # Nelder-Mead does not need gradients, so it can handle the inf outside the feasible region
    return min((minimize(objective, start, bounds=bounds, method='Nelder-Mead',
                         options={'xatol': 1e-8, 'fatol': 1e-14, 'maxiter': 20000}) for start in starts),
               key=lambda result: result.fun)


def characterize_biomass_composition(settings):
    """Entry point used by runPKP: picks the characterization type available
    in the `fuel` settings block (elemental, biochemical, or both) and
    returns the resulting reference species mass fractions."""

    elemental_analysis = settings.get('elemental_analysis',{})
    biochemical_analysis = settings.get('biochemical_analysis',{})

    if elemental_analysis and biochemical_analysis:
        characterization_type = "elem_and_bio_optimized"
    elif elemental_analysis:
        characterization_type = "elemental_analysis"
    elif biochemical_analysis:
        characterization_type = "biochemical_analysis"
    else: 
        raise ValueError("No elemental_analysis or chemical_analysis data available.")

    composition = elemental_analysis.get('composition', None)
    if elemental_analysis:
        basis = (composition or {}).get('basis')
        if basis == 'AR':
            composition = _as_received_to_daf(composition, settings['moisture'], settings['ash'])
        elif basis != 'DAF':
            raise ValueError(f"elemental_analysis.composition.basis must be 'DAF' (dry ash-free) or 'AR' (as received), got {basis}.")

    reference_species = calculate_reference_species(
        composition,
        settings['ash'],
        settings['moisture'],
        characterization_type,
        splitting_params=elemental_analysis.get('splitting_params'),
        biochem=biochemical_analysis,
        biomass_type=elemental_analysis.get('biomass_type', biochemical_analysis.get('biomass_type', None)))
    if (reference_species["PROTC"] == 0.0 and reference_species["PROTH"] == 0.0 and reference_species["PROTO"] == 0.0):
        ##Delete prot entries from the dict
        del reference_species["PROTC"]
        del reference_species["PROTH"]
        del reference_species["PROTO"]
    return reference_species