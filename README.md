# FSNA: Fuzzy Social Network Analysis

R package for directed social networks whose ties are triangular fuzzy numbers.

* **Indices of Porreca, Maturo and Ventre (2025)**, *International Journal of
  Approximate Reasoning* 176, 109319: fuzzy degree (FTOWA), best paths, fuzzy
  closeness and betweenness on best paths, computed by applying the crisp
  formula to the parameters of the ties and ranking fuzzy numbers by their
  centre of gravity. The ranking score (centre of gravity, core, Yager, graded
  mean, Liou-Wang) can be chosen in every function.
* **Exact fuzzy centralities** (Porreca and Maturo, 2026, manuscript
  submitted for publication): for every index that is continuous and non-decreasing in the
  ties (strength, OWA, t-norm closeness, Katz, spectral radius, user-defined),
  the alpha-cuts of the fuzzy centrality defined by the extension principle are
  obtained exactly from two crisp evaluations per level; triangular and Monte
  Carlo approximations for comparison.
* **Comparison of actors**: possibility and necessity degrees (suprema by
  interpolation for strength and by bisection on alpha for the other indices),
  certified orderings, rank intervals, fuzzy set of the most central actor,
  nine ranking methods.
* **Data**: a simulated survey network (`fsn_simulated()`), simulators of
  networks and questionnaires, a simulation study with the design of the paper,
  and readers of the raw mouse-tracking export. The questionnaire data of the
  two university departments analysed in the papers are not distributed with
  the package; they are available from the corresponding author upon reasonable
  request (`?fsna_questionnaire`).
* **Data collection and reporting**: `make_collector()` writes a web form
  (three-point elicitation with sliders, optional Google Sheets backend);
  `read_ratings()` and `fuzzify_responses()` import CSV/Excel answers;
  `fsna_analyze()` runs all analyses and writes figures, tables and an HTML
  report; `fsna_demo()` runs everything on a small simulated network.

```r
install.packages("FSNA")
library(FSNA)
x <- fsn_simulated()
fc <- fuzzy_centrality(x, "in_closeness_prod")
summary(fc)
plot(fc)
possibility_matrix(fc)
fsna_analyze(x, output_dir = "analysis", open = TRUE)
```

Vignettes: `vignette("FSNA")`, `vignette("exact-centrality")`,
`vignette("ijar-indices")`, `vignette("data-collection")`.

References: `citation("FSNA")`.
