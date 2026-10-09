"""Converts coal elemental analysis into mass fractions of the POLIMI
reference coal species (COAL1, COAL2, COAL3, CHAR), for use as the
`fuel.reference_species` input to a mechanism.

Ported from the POLIMI coal devolatilization model (Sommariva et al. 2010):
COAL1/COAL2/COAL3/CHAR are fixed reference species with known molecular
formulas, placed on a Van Krevelen (O/C, H/C) diagram. A target coal's own
(O/C, H/C) point is located inside one of the three triangles they form, and
its mass fractions are the barycentric weights of that triangle's vertices.
"""

import numpy as np

# (O/C, H/C) molar ratios of the POLIMI reference coal species, computed
# directly from their molecular formulas: COAL1 = C12H11, COAL2 = C14H10O1,
# COAL3 = C12H12O5, CHAR = C.
_REFERENCE_COALS = {
    'CHAR':  (0 / 1,  0 / 1),
    'COAL1': (0 / 12, 11 / 12),
    'COAL2': (1 / 14, 10 / 14),
    'COAL3': (5 / 12, 12 / 12),
}

_TRIANGLES = [
    ('CHAR', 'COAL1', 'COAL2'),
    ('CHAR', 'COAL2', 'COAL3'),
    ('COAL1', 'COAL2', 'COAL3'),
]

MW_C, MW_H, MW_O = 12, 1, 16


def characterize_coal_composition(settings):
    """Return reference species (COAL1/COAL2/COAL3/CHAR + ASH/MOIST) mass
    fractions for a coal defined by `settings['elemental_analysis']['composition']`
    (dict with C, H, O in mass %, or any consistent unit - only ratios matter;
    same key/shape as biomass's elemental_analysis) and top-level
    `settings['ash']`/`settings['moisture']` (mass fractions, 0-1)."""

    ua = settings['elemental_analysis']['composition']
    mol_C, mol_H, mol_O = ua['C'] / MW_C, ua['H'] / MW_H, ua['O'] / MW_O
    point = np.array([mol_O / mol_C, mol_H / mol_C])

    for names in _TRIANGLES:
        x0, x1, x2 = (np.array(_REFERENCE_COALS[n]) for n in names)
        a, b = np.linalg.solve(np.column_stack([x1 - x0, x2 - x0]), point - x0)
        w = np.array([1 - a - b, a, b])
        if np.all(w >= -1e-9):
            weights = dict(zip(names, np.clip(w, 0, None)))
            break
    else:
        raise ValueError(
            f"Coal ultimate analysis (O/C={point[0]:.3f}, H/C={point[1]:.3f}) "
            "is outside the COAL1/COAL2/COAL3/CHAR triangulation range.")

    ash = settings.get('ash', 0)
    moist = settings.get('moisture', 0)
    factor = 1 - ash - moist

    reference_species = {name: w * factor for name, w in weights.items()}
    reference_species['ASH'] = ash
    reference_species['MOIST'] = moist
    return reference_species
