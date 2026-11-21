## Submission

This is a new submission.

## Test environments

* local: Windows, R-release
* win-builder: R-devel and R-release
* Ubuntu Linux, R-release

## R CMD check results

0 errors | 0 warnings | 1 note

* checking CRAN incoming feasibility ... NOTE
  New submission.
  Possibly misspelled words in DESCRIPTION: Maturo, Porreca, Ventre
  (author names in the references), betweenness (a centrality index),
  suprema (plural of supremum).

## Notes for the reviewers

* The questionnaire data analysed in the cited papers are not distributed,
  because the structure of the networks could allow the identification of
  respondents. The functions that read them (fsna_questionnaire(),
  fsn_department(), read_fuzzy_questionnaire()) take the path of a local
  file; their examples and tests run only when the file is available and are
  skipped otherwise. All other examples, tests and vignettes use simulated
  data.
* Examples that take more than a few seconds are wrapped in \donttest{}.
  Examples and vignettes write only to tempdir() and do not open a browser
  unless the session is interactive.
