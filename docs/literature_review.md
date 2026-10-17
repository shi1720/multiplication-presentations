# Literature review and claim boundary

Review date: 17 October 2026. This is a focused, documented prior-art search, not a proof that no equivalent result exists anywhere.

## Closest precedents

| Source | Material checked | Established result and effect on this paper |
| --- | --- | --- |
| Cifuentes, Harris, Sturmfels, *The geometry of SDP-exactness in quadratic optimization*, Math. Programming 182 (2020), [DOI](https://doi.org/10.1007/s10107-019-01399-8), [author preprint](https://arxiv.org/abs/1804.01796) | Full preprint searched; Examples 6.1–6.2 and Remark 6.3 read closely | Generator dependence, incomplete cubic description, repair by a third quadric, and minimal-degree explanation are already known. None is claimed as new. |
| Cifuentes, Ranestad, Sturmfels, Weinstein, *Voronoi cells of varieties*, JSC 109 (2022), [DOI](https://doi.org/10.1016/j.jsc.2020.07.009), [preprint](https://arxiv.org/abs/1811.08395) | Full preprint searched; Section 4 closely read | Spectrahedral Voronoi approximations and the cubic example are prior work. Problem 4.4 motivates only the restricted contact calculation here. |
| Lindberg, Rodriguez, *Invariants of SDP exactness in quadratic programming*, JSC 122 (2024), [DOI](https://doi.org/10.1016/j.jsc.2023.102258), [preprint](https://arxiv.org/abs/2211.05645) | Full preprint searched; Section 2 closely read | Same-ideal generator dependence and invariance hypotheses are established. The homogeneous complete-intersection hypothesis fails for the chain family. Our affine complete intersections do not contradict their theorem. |
| Blekherman, Smith, Velasco, *Sums of squares and varieties of minimal degree*, JAMS 29 (2016), [DOI](https://doi.org/10.1090/jams/847), [preprint](https://arxiv.org/abs/1308.0751) | Main theorem and bibliographic record | Full quadratic closure on the rational normal curve is classical. The manuscript supplies a direct univariate specialization, labeled accordingly. |
| Mevissen, Kojima, *SDP relaxations for quadratic optimization problems derived from polynomial optimization problems*, APJOR 27 (2010), [DOI](https://doi.org/10.1142/S0217595910002533), [institutional report index](https://www.is.c.titech.ac.jp/report-b/) | Author/institutional record and original preprint abstract | Monomial transformation heuristics and relaxation quality are not new research directions. Our distinction is the exact exposed-face dimension at fixed size. |
| Karia, Adjiman, Chachuat, *Assessment of a two-step approach for global optimization of mixed-integer polynomial programs using quadratic reformulation*, C&CE 165 (2022), [DOI](https://doi.org/10.1016/j.compchemeng.2022.107909) | Publisher text, Sections 1–3; Crossref metadata | Equivalent quadratic reformulations and redundant constraints are established. We do not claim an industrial solver or general reformulation superiority. |
| Henrion, *Geometry of exactness of moment-SOS relaxations for polynomial optimization*, Optimization Letters 19 (2025), [DOI](https://doi.org/10.1007/s11590-025-02202-8), [author manuscript](https://homepages.laas.fr/henrion/papers/exactness.pdf) | Theorem and author publication record | Broader exactness-cone descriptions are known. This paper works out a particular factorization classification. |
| Cai, Zardini, *Optimizing Lyapunov Certificates via Stability-Preserving Quadratization for Polynomial Systems*, [arXiv:2609.18871](https://arxiv.org/abs/2609.18871) | Current abstract | Recent representation/gauge optimization for dynamical Lyapunov certificates; a different object from nearest-point Voronoi faces. No claim that representation-aware certificate optimization in general is new. |

## Search families

Searches included combinations of: Euclidean distance optimization / quadrics / normal bundle / semidefinite exactness; Voronoi cells / master spectrahedron / redundant generators; twisted cubic / Voronoi / SDP; moment curve / Shor / quadratic optimization; rational normal curve / nonnegative quadratics / sums of squares; quadratic reformulation / monomial decomposition / semidefinite; multiplication / SDP-exactness; focal / Voronoi / semidefinite; invariants of SDP exactness; and recent stability-preserving quadratization.

Primary sources and author/institutional manuscripts were preferred. Bibliographic details were checked against publisher, arXiv, author, and Crossref records. Paywall/access failures were not treated as evidence of absence.

## What survived the search

The proposed contribution is the exact dimension formula for *every recursive one-product-per-coordinate presentation* of the full affine moment curve, its generating polynomial, the full-dimensional normal-slice volume bound, and the explicit anisotropic contact/unbounded-gap consequences. The reviewed sources contain the cubic precedent and general tools, but no matching statement of this entire classification was located. This remains a provisional literature judgment; a specialist may know an equivalent result or a stronger existing formulation.

The proof itself establishes the stated mathematics independently of that judgment. The paper deliberately includes the known cubic case only as context and does not sell full closure as a new solution.

## Practical significance and limits

At identical problem size, one can lose an arbitrarily large number of exposed-face dimensions. This gives a useful analytical stress test for polynomial-model compilers. It does not solve general polynomial optimization, demonstrate a new LLM architecture, establish patentability, or imply a high acceptance probability. The intended presentation is a focused mathematical contribution suitable for discussion at an algebraic distance workshop.
