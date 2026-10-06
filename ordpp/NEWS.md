# ordpp 0.1.0

* Seeded calls restore the caller's random number state on exit, and `seed`
  defaults to `NULL` in `heldout_test_ordpp()`, `permutation_ordpp()` and
  `stability_ordpp()`.

* Default `maxit` raised to 500: item-specific fits often need more than 150.
* `simulate_ordpp()` gains `strength` for weak-dependence power studies.
* Characteristic functions are evaluated on distinct values (much faster).

* New `index = "shape"` option compares the variance-standardized projection
  with the independence reference, removing the pairwise-covariance component.
* New `decompose_ordpp()` reports Var(Z)/Var(Z0) with the full and shape indices.
* New `heldout_test_ordpp()`: exact permutation test along a frozen map.
* `permutation_ordpp()` accepts `fit =` so the null fits reuse its exact settings.
* New `stability_ordpp()` for (cluster) bootstrap stability of loadings and scores.
* New `as_ordpp_map()` wraps external scores/loadings (e.g. optimal-scaling PCA)
  for evaluation with the same held-out tools.
* `simulate_ordpp()` gains `scenario =` ("independent", "rare", "two_factor",
  "response_style").
* Fitted objects store transformed data, weights, settings, the call and the
  number of gap parameters on their bounds.

# ordpp 0.0.1

* Initial prototype.
