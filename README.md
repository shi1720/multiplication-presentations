# Multiplication presentations and focal-face collapse

Research artifact accompanying **Multiplication Presentations and Focal-Face Collapse in Algebraic Distance Certificates**, Shivam Gupta, 17 October 2026.

[Read the paper](paper/main.pdf) · [LaTeX source](paper/main.tex) · [Literature and novelty boundary](docs/literature_review.md)

For the affine moment curve `(t,t²,…,t^d)`, every recursive multiplication presentation has the same `d` variables and `d−1` quadratic equations. The paper classifies how the factor choices change the dimension of a certified exposed Voronoi face. It also proves a near-face volume-loss bound and gives counterexamples to universal dominance of balanced products.

The dependence of SDP exactness on generators, the cubic example, and exactness after full quadratic closure are established results. They are credited as such. The new claim is the stated classification and its quantitative consequences for the specified family, subject to the literature boundary in the paper. This is a research manuscript, not an accepted publication.

## Reproduce

Python 3.9+; the recorded environment used Python 3.9.6 on macOS arm64. A modern environment may require updating dependency pins together; use the recorded pins for the reported results.

```sh
python3 -m venv .venv
. .venv/bin/activate
python -m pip install -r requirements.txt
python tests/verify_exact.py
python tests/verify_weighted.py
python experiments/run.py
python experiments/near_face.py
python experiments/figures.py
```

Build the manuscript with Tectonic from its source directory:

```sh
cd paper
tectonic main.tex
```

The mathematical tests and experiments require no network, external dataset, API, GPU, or commercial solver. Tectonic may download TeX packages during the first build.

## Results and guarantees

- All 3,655 multiplication presentations in dimensions 3–10: symbolic distance identity, exact rational PSD checks, negative principal-minor witnesses for every forced coordinate, and agreement with the dimension enumerator.
- 400 exact weighted-metric identity/PSD checks.
- Exact rational grid: 201 of 315 inputs belong to the closed cell. Chain certifies 148; balanced 136; union 191. Both sparse systems miss valid inputs. No false positives.
- 144 numerical full-closure SDP comparisons: maximum boundary difference `5.15e-9`; all statuses `optimal`. Tiny negative eigenvalues remain numerical residuals, not exact certificates.
- 30,000 seeded full-dimensional near-face inputs: balanced certifies 30,000; chain 16,034, for the deliberately defined sampling geometry.
- Synthetic four-channel consistency example: all 205 inputs have the same unique nearest state; balanced certifies 205 and chain 20. This is not measured industrial data.

`results/*.csv` contain raw rows; JSON files record aggregates and the environment. `figures/*.pdf` are vector figures; PNGs are previews. `experiments/run.py` uses an exact rational root-parity oracle separately from floating polynomial roots and the CVXPY SDP model.

## Code entry points

`src/presentations.py` implements recursive presentations, centered matrices, exact rational PSD testing, exact univariate membership, scalar boundary calculations, and the full quadratic-closure boundary SDP.

```python
import sys
sys.path.insert(0, 'src')
from presentations import exact_matrix, rational_psd, splits
import sympy as sp
u = [0, sp.Rational(1, 2), 0, -2]
assert not rational_psd(exact_matrix(u, splits(4, 'chain')))
assert rational_psd(exact_matrix(u, splits(4, 'balanced')))
```

## Scope

The face theorem concerns the origin of the complete affine moment curve and recursive factorizations, with a positive-diagonal-metric extension. Full closure exactness is classical and is not a new contribution. The sparse presentations are incomparable. This artifact makes no runtime superiority, patentability, industrial validation, or acceptance claim.

## License

Code: MIT (see LICENSE). Manuscript and original figures: CC BY 4.0, copyright 2026 Shivam Gupta. Third-party works are referenced and retain their original licenses.
