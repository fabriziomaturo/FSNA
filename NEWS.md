# FSNA 0.2.2

First public release.

* Fuzzy degree, closeness and betweenness of Porreca, Maturo and Ventre (2025): FTOWA of the ties, best paths selected by a ranking score, betweenness on best paths with tied paths sharing the count (Freeman fraction). Every function accepts the ranking score (centre of gravity, core, Yager, graded mean, Liou-Wang). `variant = "definition"` follows the definitions of the paper; `variant = "alternative"` changes the path strength, the treatment of equal scores and the normalization, and each choice can be set individually (`path_tfn`, `tie_break`, `normalize`, `direct_rule`, `ties`).
* Exact fuzzy centralities of Porreca and Maturo (2026, manuscript submitted for publication; preprint forthcoming on arXiv) for every index that is continuous and non-decreasing in the ties (strength, OWA, t-norm closeness, Katz, spectral radius, user-defined indices): alpha-cuts from two crisp evaluations per level, with triangular and Monte Carlo approximations for comparison.
* Comparison of actors: possibility and necessity degrees with suprema located by interpolation for affine indices and by bisection on alpha otherwise, returned with numerical brackets; certified orderings; rank intervals; fuzzy set of the most central actor; nine ranking methods.
* Unknown ties for pairs without information, with interval cores.
* Data: a simulated survey network (`fsn_simulated()`), simulators of fuzzy networks and of questionnaire ratings, a simulation study, readers of the raw mouse-tracking export. The questionnaire data of the two departments are not distributed; `fsna_questionnaire()` and `fsn_department()` read them from a local path.
* Data collection and reporting: web form with three-point elicitation and an optional Google Sheets backend (`make_collector()`), importers of CSV and Excel answers, `fsna_analyze()` for figures, tables and an HTML report.
