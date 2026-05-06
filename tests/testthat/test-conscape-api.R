test_that("conscape_api_slots exposes the ConScape crosswalk", {
  slots <- conscape_api_slots()

  expect_s3_class(slots, "data.frame")
  expect_true(all(c("slot", "run_conscape", "conscape_prep", "julia") %in% names(slots)))
  expect_true(all(c("target_qualities", "source_qualities", "affinities") %in% slots$slot))
})

test_that("conscape_sensitivity validates and names sensitivity requests", {
  spec <- conscape_sensitivity(wrt = c("Q", "A&C=f(A)"), unitless = TRUE)

  expect_s3_class(spec, "ConScapeSensitivitySpec")
  expect_equal(spec$wrt, c("Q", "A&C=f(A)"))
  expect_equal(
    getFromNamespace("conscape_sensitivity_id", "ConScapeRtools")("A&C=f(A)", TRUE),
    "elasticity_affinity_cost_linked"
  )
  expect_error(conscape_sensitivity(wrt = "bad"), "Unsupported")
})

test_that("ConScape API normalizers preserve compatibility aliases", {
  normalize_metrics <- getFromNamespace("normalize_conscape_metrics", "ConScapeRtools")
  normalize_cost <- getFromNamespace("normalize_conscape_cost_function", "ConScapeRtools")
  normalize_conn <- getFromNamespace("normalize_conscape_connectivity_function", "ConScapeRtools")
  output_specs <- getFromNamespace("conscape_output_specs", "ConScapeRtools")

  expect_equal(normalize_metrics(c("btwn", "fcon")), c("betweenness_kweighted", "connected_habitat"))
  expect_equal(normalize_cost("x -> -log(x)"), "minuslog")
  expect_equal(normalize_conn("survival"), "survival_probability")

  specs <- output_specs(
    metrics = "connected_habitat",
    sensitivity = conscape_sensitivity(wrt = "Q")
  )
  expect_true(all(c("connected_habitat", "elasticity_quality") %in% names(specs)))
  expect_equal(specs$connected_habitat$dir, "fcon")
})
