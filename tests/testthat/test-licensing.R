test_that("installed metadata and notices credit the adapted Julia routines", {
  metadata <- read.dcf(system.file("DESCRIPTION", package = "ConScapeRtools"))
  expect_identical(unname(metadata[1, "License"]), "GPL-2")
  expect_match(metadata[1, "Authors@R"], "ConScape.jl contributors", fixed = TRUE)

  notice_file <- system.file("COPYRIGHTS", package = "ConScapeRtools")
  expect_true(file.exists(notice_file))
  notice <- paste(readLines(notice_file, warn = FALSE), collapse = "\n")
  expect_match(notice, "src/randomizedshortestpath.jl", fixed = TRUE)
  expect_match(notice, "1534b7983fcf3f52d1e84e2d148625342842d614", fixed = TRUE)
  expect_match(notice, "b0d387e05e5135e7375e33718cdcb466117c7bee", fixed = TRUE)
  expect_match(notice, "ConScapeRtools adaptations, 2026-10-08", fixed = TRUE)
  expect_match(notice, "Copyright (c) 2024 ConScapeRtools authors", fixed = TRUE)
  expect_match(notice, "Permission is hereby granted, free of charge", fixed = TRUE)
  expect_match(notice, "THE SOFTWARE IS PROVIDED", fixed = TRUE)

  helper <- readLines(system.file("extdata", "conscape.jl", package = "ConScapeRtools"),
                       warn = FALSE)
  expect_match(helper[1], "SPDX-License-Identifier: GPL-2.0-only", fixed = TRUE)
})
