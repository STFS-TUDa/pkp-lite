# pkp-lite

pkp-lite fits simplified pyrolysis kinetic models (SFOR, C2SM) to detailed
[Cantera](https://cantera.org/) reaction mechanism simulations. Give it a
fuel definition and a set of time-temperature profiles, and it will:

1. run a detailed, ODE-based reactor simulation using your Cantera mechanism, and
2. fit the coefficients of a simplified model (SFOR or C2SM) that best match the detailed results.

This is useful for producing lightweight kinetic models (e.g. for CFD) that
approximate a detailed pyrolysis mechanism without needing the full
mechanism at simulation time.

## Installation

Requires Python 3.10+.

```bash
python3.10 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

## Getting a mechanism

pkp-lite needs a Cantera-format solid-phase reaction mechanism (referenced
via `fuel.mechanism` in the input YAML). No mechanism is bundled with this
repo, but the [CRECK Kinetic-Mechanisms](https://github.com/CRECKMODELING/Kinetic-Mechanisms)
repo publishes biomass and coal solid-phase pyrolysis mechanisms in CHEMKIN
format, with species naming (`CELL`/`XYHW`/.../`COAL1`/`COAL2`/`COAL3`/`CHAR`)
that matches what `utility_functions.py` and `characterize_coal.py`
already expect:

- Biomass: [`Solid-Phase/Biomass`](https://github.com/CRECKMODELING/Kinetic-Mechanisms/tree/main/Solid-Phase/Biomass)
  — use the `Solid_only/` subfolder (solid pyrolysis reactions only) unless
  you also need the much larger secondary gas-phase mechanism in `Solid_gas/`.
  If you use it, cite [[1]–[5]](#references).
- Coal: [`Solid-Phase/Coal`](https://github.com/CRECKMODELING/Kinetic-Mechanisms/tree/main/Solid-Phase/Coal).
  If you use it, cite [[6]–[10]](#references).

### Converting to Cantera format

CRECK's files are CHEMKIN format (a kinetics file + a thermo file), but they
can't be fed to Cantera's converter as-is: they use CRECK's own
`MATERIAL <name>` / `SOLID <species...>` header instead of a standard
CHEMKIN `ELEMENTS`/`SPECIES` block, and that header doesn't list species that
only appear as reaction *products* (e.g. `CO`, `CO2`, `H2O`). The scripts
below build a corrected copy of the file and then run Cantera's `ck2yaml`
converter (installed with `cantera`, e.g. at `.venv/bin/ck2yaml`) on it.
Run them from the folder you downloaded the CRECK files into.

**Biomass** — needs `BIO.solid` (kinetics) and `BIO.CKT` (thermo) from `Solid_only/`:

```bash
# Collect every species referenced in the REACTIONS block. Species that
# never react (e.g. ASH) don't appear there, so add them by hand.
awk '/^REACTIONS/{f=1} f' BIO.solid | sed 's/!.*//' \
  | grep -oE '[A-Za-z][A-Za-z0-9_]*' \
  | grep -vE '^(REACTIONS|END|E|DUPLICATE|FORD|LOW|TROE|REV)$' \
  | sort -u > species.txt
echo ASH >> species.txt

# Write a valid CHEMKIN file: a real ELEMENTS/SPECIES header, followed by
# the original REACTIONS block unchanged.
{
  echo "ELEMENTS"; echo "C H O N K"; echo "END"
  echo "SPECIES"; cat species.txt; echo "END"
  awk '/^REACTIONS/{f=1} f' BIO.solid
} > BIO_fixed.inp

ck2yaml --input=BIO_fixed.inp --thermo=BIO.CKT --output=BIO_solid.yaml
```

**Coal** — needs `COAL_2003.solid` (kinetics) and `COAL_2003.tdc` (thermo):

```bash
awk '/^REACTIONS/{f=1} f' COAL_2003.solid | sed 's/!.*//' \
  | grep -oE '[A-Za-z][A-Za-z0-9_]*' \
  | grep -vE '^(REACTIONS|END|E|DUPLICATE|FORD|LOW|TROE|REV)$' \
  | sort -u > species.txt
echo ASH >> species.txt

{
  echo "ELEMENTS"; echo "C H O N S"; echo "END"
  echo "SPECIES"; cat species.txt; echo "END"
  awk '/^REACTIONS/{f=1} f' COAL_2003.solid
} > COAL_fixed.inp

# --transport is omitted: pkp-lite's 0-D reactor doesn't need transport
# data, and the CRECK transport file doesn't cover the solid species anyway.
# --permissive is needed because COAL_2003.tdc has a few duplicate/malformed
# thermo entries for unrelated gas-phase species.
ck2yaml --input=COAL_fixed.inp --thermo=COAL_2003.tdc --output=COAL_2003.yaml --permissive
```

Both recipes were run end to end against the real CRECK files while writing
this README, and the resulting `.yaml` files load cleanly in Cantera.

Point `fuel.mechanism` in your `input.yml` at the resulting `.yaml` file.

## Usage

```bash
python scripts/runPKP input.yml [-n N]
```

- `-n N` — number of parallel workers. Only matters for `fit_method: evolution`
  (passed to `multiprocessing.Pool`); ignored by `least_squares`/`least_squares_log`,
  which run single-threaded.

Output always goes to `results/<input-basename>_<YYYY-MM-DD>_<NN>/`
(auto-incrementing `NN` per day):

| File | Contents |
|---|---|
| `<input>.yml` | copy of the input file used for this run |
| `run{i}.csv` | detailed Cantera simulation for operating condition `i`: time, every mechanism species' mass fraction, temperature, plus aggregated `raw`/`metaplast`/`char`/`light_gas`/`tar`/`solid`/`volatiles` columns |
| `run{i}_model.csv` | the fitted simplified model's prediction for the same run (`t`, `volatiles`, `char`, `raw`, `T`) — compare directly against `run{i}.csv` |
| `coefficients.yml` | `{model: SFOR\|C2SM, coefficients: {...}}` — the fitted result |

### Species classification

The `raw`/`metaplast`/`char`/`light_gas`/`tar` columns in `run{i}.csv` (and
the `volatiles` = `tar + light_gas` / `solid` = `raw + metaplast + char`
totals fit against) come from lumping individual mechanism species together,
per [`pkplite/species_classification.yml`](pkplite/species_classification.yml).
These lists are tied to exactly which mechanism version you're using — the
bundled defaults resemble the current POLIMI convention (CRECK `COAL_2003`
for coal, the `Solid_only` biomass mechanism for biomass; see
[Getting a mechanism](#getting-a-mechanism)). If you use a different or
future mechanism version with different species names, copy that file,
adjust it, and point `fuel.species_classification` at your copy.

pkp-lite only ever fits SFOR/C2SM against the combined `volatiles` trace
(`tar + light_gas`) — there's no separate tar vs. light-gas release
kinetics. This mirrors the original PKP method, which also assumed
volatiles release as light gases and tar (heavy hydrocarbons, liquid or
solid at room temperature) in a constant ratio throughout devolatilization
— a simplification the original authors flagged as inaccurate in general
(tar release tends to dominate at lower temperatures, light gas at higher
ones [[12]](#references)) but reasonable for fast devolatilization regimes.

### Fuel composition: three ways to specify it

1. **Direct**: `fuel.reference_species: {SPECIES_NAME: mass_fraction, ...}` — used as-is, no characterization step.

2. **Elemental analysis** — same key, same underlying method, for both biomass and coal:

   `elemental_analysis: {composition: {C, H, O}}` + top-level `ash`, `moisture` (0–1 fractions).

   ("Elemental analysis" and "ultimate analysis" are the same thing — C/H/O/N/S composition by mass — just biomass- vs. coal-community jargon for it. pkp-lite uses one key, `elemental_analysis`, for both.)

   **What this actually does**: your C/H/O composition is converted into a point on a Van Krevelen diagram (O/C vs. H/C molar ratio). A fixed set of reference species with known molecular formulas sit at fixed points on that same diagram — `COAL1`/`COAL2`/`COAL3`/`CHAR` for coal, `CELL`/`XYHW`/`LIGC`/`LIGH`/`LIGO`/`TANN`/`TGL` for biomass — forming a mesh of triangles between neighboring vertices (the CRECK/POLIMI characterization method: [[10]](#references) for coal, [[11]](#references) for biomass). Your fuel's point is located inside one of those triangles, and its reference-species mass fractions are set so that a mix of that triangle's vertices reproduces your O/C and H/C ratios exactly (barycentric interpolation).

   Both fuel types raise a `ValueError` rather than silently extrapolating if your composition falls outside the characterizable range, just via different checks: coal (`characterize_coal.py`) explicitly checks all 3 candidate triangles for containment up front. Biomass (`characterize_biomass.py`) solves the equivalent system directly via 5 splitting parameters (`alpha`/`beta`/`gamma`/`delta`/`epsilon`, with sensible defaults) and instead validates afterward that every resulting reference-species fraction is within [0, 1].

   Biomass additionally requires `basis: DAF|AR` inside `elemental_analysis.composition`: whether C/H/O/N are dry ash-free (`DAF`) or as received (`AR`, H includes the hydrogen of the moisture). `AR` is converted to DAF using `moisture` and `ash`, which are always as received. For biomass, O is always taken as 1 − C − H − N on DAF basis, so a given `O` is not used.

   Biomass also takes `biomass_type: hardwood|softwood|grass` inside `elemental_analysis` (changes which hemicellulose-like species — `XYHW`/`GMSW`/`XYGR` — is used), and an optional `N` in `composition` — if present and nonzero, part of the fuel mass is assigned to `PROTC`/`PROTH`/`PROTO` (protein) reference species instead, proportional to the nitrogen content; omitted entirely from the output if there's no nitrogen.

3. **Biomass-only, alternative**: `biochemical_analysis: {cellulose, hemicellulose, lignin, extractives, proteins}` — simpler equal-split method, no triangulation. Can be combined with `elemental_analysis`: pkp-lite then optimizes the elemental method's splitting parameters to also match the biochemical numbers.

### Multiple heating rates / operating conditions

`operating_conditions.runs` is a list of runs, each its own piecewise-linear
`[[t0,T0],[t1,T1],...]` time-temperature profile — e.g. several different
heating rates or peak temperatures. **All runs are fit together against one
shared set of coefficients** — you get a single SFOR/C2SM parameter set that
best reproduces every given condition at once (the standard multi-heating-rate
kinetic-fitting approach), not a separate fit per run.

### Fitting options

| `fit.model` | coefficients |
|---|---|
| `SFOR` | `alpha, A, Ea` |
| `C2SM` | `alpha1, A1, Ea1, alpha2, A2, Ea2` |

| `fit.fit_method` | behavior |
|---|---|
| `evolution` (default) | DEAP mu+lambda EA. Settings: `npop`, `ngen`, `mu`, `cxpb`, `mutpb` (defaults: `npop=20`, others fall back to `npop`/0.5), `final_error_weight` (default 0 — extra penalty on the final-time error, on top of the mean-squared error over the whole trace). |
| `least_squares` | scipy TRF. Setting: `verbose` (0/1/2, scipy convention). |
| `least_squares_log` | same as `least_squares`, but any bound pair with `max/min > 100` is fit in log10-space — use this instead of `least_squares` when `A`/`Ea` bounds span many orders of magnitude; converges much more reliably. |

Both methods take `fit.bounds.min`/`fit.bounds.max` (in the model's
coefficient order above) and an optional `fit.parameters_init` warm start.

Following Vascellari et al. [[12]](#references)'s original PKP objective function, each
operating condition's error is normalized by that run's own volatile-yield
range and averaged with equal weight across all `operating_conditions.runs`
— so a heating-rate condition with a larger absolute yield swing doesn't
dominate the fit just because of its scale.

`evolution` (a genetic algorithm) was the original PKP method's calibration
approach [[12]](#references), chosen specifically because the objective function can have
several local minima when multiple heating rates are fit at once — gradient-
based methods aren't reliable there. `least_squares` is pkp-lite's own
addition: faster, but more sensitive to the initial guess and to local
minima, so it's worth checking convergence (e.g. trying a couple of
different `parameters_init`, or cross-checking against an `evolution` run)
rather than trusting a single run blindly.

**Watch out for:** a result where a coefficient lands exactly on one of its
`bounds` (e.g. `alpha` at its max), or an `Ea`/`A` that's physically
implausible for the fuel, is a sign of a bad local minimum, not a good fit —
`least_squares` is the likelier culprit, but `evolution` can land on one too
(a small `npop`/`ngen`, or just an unlucky run). Don't take such a result at
face value; it's worth simply rerunning (a GA restart explores a different
random population; a gradient fit needs a different `parameters_init`) to
see if the result changes.

[`analyze_fit.ipynb`](analyze_fit.ipynb) (needs `matplotlib`) plots the
detailed vs. fitted-model `volatiles` trace for each run in a `results/<run>/`
folder side by side, and reports the per-run fit error using the same
normalization as the fitting objective above — point it at a results folder
to sanity-check a fit visually rather than from `coefficients.yml` alone.

### `input.yml` structure

```yaml
fuel:
  name: ...
  type: biomass  # or coal
  reference_species: {...}  # or elemental_analysis: {composition: {C, H, O}, ...} + top-level
                            # ash/moisture; biomass can also use biochemical_analysis
  HHV: ...
  rho_dry: ...
  mechanism: path/to/mechanism.yaml
  species_classification: path/to/species_classification.yml  # optional, overrides the
                            # bundled pkplite/species_classification.yml (see below)

fit:
  model: SFOR          # or C2SM
  fit_method: least_squares   # or evolution / least_squares_log
  bounds:
    min: [0.0, 1e0, 1e4]     # SFOR: [alpha, A, Ea]
    max: [1.0, 1e12, 1e9]
  parameters_init: [0.6, 1e5, 1e8]

operating_conditions:
  runs:
    - [[t0, T0], [t1, T1], ...]  # piecewise T(t) profile per run
```

## References

[1]–[10] are citations for the CRECK solid-phase mechanisms (see
[Getting a mechanism](#getting-a-mechanism)), kept in sync with each folder's
own `README.md` on the [CRECK repo](https://github.com/CRECKMODELING/Kinetic-Mechanisms),
which is the source of truth if this list drifts out of date. [11] is the
biomass elemental-analysis characterization method (see
[Fuel composition](#fuel-composition-three-ways-to-specify-it)). [12] is the
original PKP fitting method pkp-lite is based on, and [13] the DEAP library
its `evolution` fit method uses (see [Fitting options](#fitting-options) and
[Acknowledgments](#acknowledgments)).

**Biomass**

[1] Afessa, Million et al. "Pyrolysis of large biomass particles: model validation and application to Coffee Husks valorization." *Journal of Analytical and Applied Pyrolysis* (2025): 107028. [DOI](https://doi.org/10.1016/j.jaap.2025.107028)

[2] Debiagi, Paulo et al. "Cellulose pyrolysis kinetic model: Detailed description of volatile species." *Proceedings of the Combustion Institute* 40 (2024): 105651. [DOI](https://doi.org/10.1016/j.proci.2024.105651)

[3] Zou, Xun et al. "Kinetic insights into the high-temperature pyrolysis of biomass: Experimental and modeling study." *Journal of Analytical and Applied Pyrolysis* 179 (2024): 106463. [DOI](https://doi.org/10.1016/j.jaap.2024.106463)

[4] Ranzi, Eliseo et al. "Mathematical Modeling of Fast Biomass Pyrolysis and Bio-Oil Formation. Note I: Kinetic Mechanism of Biomass Pyrolysis." *ACS Sustainable Chemistry & Engineering* 5 (2017): 2867–2881. [DOI](https://doi.org/10.1021/acssuschemeng.6b03096)

[5] Ranzi, Eliseo et al. "Mathematical Modeling of Fast Biomass Pyrolysis and Bio-Oil Formation. Note II: Secondary Gas-Phase Reactions and Bio-Oil Formation." *ACS Sustainable Chemistry & Engineering* 5 (2017): 2882–2896. [DOI](https://doi.org/10.1021/acssuschemeng.6b03098)

**Coal**

[6] Debiagi, Paulo E. A. et al. "Calibration and validation of a comprehensive kinetic model of coal conversion in inert, air and oxy-fuel conditions using data from multiple test rigs." *Fuel* 290 (2021): 119682. [DOI](https://doi.org/10.1016/j.fuel.2020.119682)

[7] Maffei, Tiziano et al. "Experimental and modeling study of single coal particle combustion in O2/N2 and Oxy-fuel (O2/CO2) atmospheres." *Combustion and Flame* 160 (2013): 2559–2572. [DOI](https://doi.org/10.1016/j.combustflame.2013.06.002)

[8] Maffei, Tiziano et al. "Predictive one step kinetic model of coal pyrolysis for CFD applications." *Proceedings of the Combustion Institute* 34 (2013): 2401–2410. [DOI](https://doi.org/10.1016/j.proci.2012.08.006)

[9] Maffei, Tiziano et al. "A predictive kinetic model of sulfur release from coal." *Fuel* 91 (2012): 213–223. [DOI](https://doi.org/10.1016/j.fuel.2011.08.017)

[10] Sommariva, Samuele et al. "A predictive multi-step kinetic model of coal devolatilization." *Fuel* 89 (2010): 318–328. [DOI](https://doi.org/10.1016/j.fuel.2009.07.023)

**Characterization method** (referenced from [Fuel composition](#fuel-composition-three-ways-to-specify-it) above, not from a CRECK repo README)

[11] Debiagi, Paulo E. A. et al. "Extractives Extend the Applicability of Multistep Kinetic Scheme of Biomass Pyrolysis." *Energy & Fuels* 29 (2015): 6544–6555. [DOI](https://doi.org/10.1021/acs.energyfuels.5b01753)

**Fitting method**

[12] Vascellari, Michele, Arora, R., Pollack, M., and Hasse, C. "Simulation of Entrained Flow Gasification with Advanced Coal Conversion Submodels. Part 1: Pyrolysis." *Fuel* 113 (2013): 654–669. [DOI](https://doi.org/10.1016/j.fuel.2013.06.014)

[13] Fortin, Félix-Antoine, De Rainville, François-Michel, Gardner, Marc-André, Parizeau, Marc, and Gagné, Christian. "DEAP: Evolutionary Algorithms Made Easy." *Journal of Machine Learning Research* 13 (2012): 2171–2175. [Paper](https://www.jmlr.org/papers/v13/fortin12a.html)

## Authors

Leon Loni Berkel, Pascal Steffens, Hendrik Nicolai, Christian Hasse

If you use pkp-lite itself in your research, please cite it — see
[CITATION.cff](CITATION.cff).

## Acknowledgments

pkp-lite's reactor/fitting pipeline is based on the PKP code originally
developed by Michele Vascellari et al. [[12]](#references). The `evolution`
fit method is built on [DEAP](https://github.com/DEAP/deap)
(LGPL-3.0) [[13]](#references).

## Contact

pkp-lite is maintained by a small group at STFS, TU Darmstadt, in our spare
time alongside other work — response times may vary. Found a bug, have a
question, or want to contribute? Open an issue, or reach out directly:

- Leon Loni Berkel — berkel@stfs.tu-darmstadt.de
- Pascal Steffens — steffens@stfs.tu-darmstadt.de

## License

GPL-3.0-or-later — see [LICENSE.txt](LICENSE.txt).
