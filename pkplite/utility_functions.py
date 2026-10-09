import pandas as pd
import os
import numpy as np
import yaml

_DEFAULT_SPECIES_CLASSIFICATION = os.path.join(
    os.path.dirname(__file__), 'species_classification.yml')

def save_data(results, basename, add_info=''):
    """Write each run's result DataFrame to results/<basename>/run<i><add_info>.csv."""

    for i,df in enumerate(results):
        df.to_csv(os.path.join("results",basename,f"run{str(i)}{add_info}.csv"))


def calculate_volatiles_and_solids_from_species(results, species_classification_path=None):
    """Aggregate per-species mechanism output into lumped volatiles/solid
    components (raw, metaplast, char, light_gas, tar). Species lists come
    from species_classification.yml (bundled default, or a path passed in
    via `species_classification_path` / `fuel.species_classification` in
    input.yml) and are mechanism-specific - see that file's comments. Coal
    and biomass lists are merged; a species absent from a run's columns is
    simply skipped."""

    path = species_classification_path or _DEFAULT_SPECIES_CLASSIFICATION
    with open(path) as f:
        classification = yaml.safe_load(f)

    components = {}
    for key in ('raw', 'metaplast', 'char', 'light_gas', 'tar'):
        merged = classification['coal'][key] + classification['biomass'][key]
        components[key] = list(np.unique(np.array(merged)))

    for res in results:
        for v in components:
            res[v] = res[[spec for spec in components[v] if spec in res.columns]].sum(axis=1)
        res['solid'] = res[['metaplast', 'char', 'raw']].sum(axis=1)
        res['volatiles'] = res[['tar', 'light_gas']].sum(axis=1)

    return results