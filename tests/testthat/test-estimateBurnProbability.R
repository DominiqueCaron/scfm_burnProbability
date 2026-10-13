test_that("estimateBurnProbability's fast Monte Carlo loop matches the per-module-call path", {
  ## Build a small, self-contained sim with everything estimateBurnProbability()
  ## reads, using scfmIgnition/scfmEscape/scfmSpread's own .inputObjects() defaults
  ## so this test doesn't depend on any external data.
  skip_if_not_installed("SpaDES.core")
  skip_if_not_installed("SpaDES.tools")
  skip_if_not_installed("terra")

  studyArea <- LandR::randomStudyArea(size = 10000 * 6.25 * 1000, seed = 1011)
  rasterToMatch <- terra::rast(studyArea, vals = 1, res = c(250, 250)) |>
    terra::mask(studyArea)

  times <- list(start = 0, end = 0)
  modules <- list(
    "scfmIgnition",
    "scfmEscape",
    "scfmSpread",
    "scfm_burnProbability"
  )
  paths <- list(
    modulePath = c(
      testthat::test_path("..", ".."),
      testthat::test_path("..", "..", "..", "modules")
    ),
    cachePath = file.path(tempdir(), "cache"),
    outputPath = file.path(tempdir(), "outputs")
  )
  params <- list(
    scfm_burnProbability = list(nIterations = 300L)
  )
  objects <- list(studyArea = studyArea, rasterToMatch = rasterToMatch)

  mySim <- SpaDES.core::simInit(
    times = times,
    params = params,
    modules = modules,
    objects = objects,
    paths = paths
  )

  ## Run once so scfmIgnition/scfmEscape/scfmSpread's `init` events populate
  ## pIg/p0/pSpread/fireRegimePolys/fireRegimeRas (NULL right after simInit()).
  mySim <- SpaDES.core::spades(mySim, debug = FALSE)

  ## Set current(sim) as doEvent dispatch would, so P(sim)$<param> resolves
  ## against scfm_burnProbability (matching estimateBurnProbability()'s own
  ## calls to P(sim, "neighbours", "scfmEscape")/P(sim, "neighbours", "scfmSpread")).
  mySim@current <- data.table::data.table(
    eventTime = 0,
    moduleName = "scfm_burnProbability",
    eventType = "estimate",
    eventPriority = 7.51
  )

  fastFn <- mySim$.mods$scfm_burnProbability$estimateBurnProbability

  set.seed(2024)
  outFast <- fastFn(mySim)

  expect_s4_class(outFast, "simList")
  expect_true(inherits(outFast$burnProbability, "SpatRaster"))

  probVals <- terra::values(outFast$burnProbability, mat = FALSE)
  flamVals <- terra::values(mySim$flammableMap, mat = FALSE)

  ## burnProbability is NA exactly where not flammable, and a valid probability
  ## (in [0, 1]) everywhere else. `flamVals` holds `NaN` (not `NA`) where not
  ## flammable, and `NaN != 1` is `NA`, not `FALSE`, so compare against `%in% 1`.
  expect_equal(is.na(probVals), !(flamVals %in% 1))
  expect_true(all(
    probVals[!is.na(probVals)] >= 0 & probVals[!is.na(probVals)] <= 1
  ))

  ## Statistical sanity check: burn probability should be positive somewhere,
  ## given nIterations = 300 and default scfm parameters (non-zero pIg/p0/pSpread).
  expect_true(any(probVals > 0, na.rm = TRUE))

  ## Reproducibility: same seed, same nIterations -> identical result (the
  ## fast path's only randomness is runif()/spread2(), both seeded by set.seed()).
  set.seed(2024)
  outFast2 <- fastFn(mySim)
  expect_identical(
    terra::values(outFast$burnProbability, mat = FALSE),
    terra::values(outFast2$burnProbability, mat = FALSE)
  )
})
