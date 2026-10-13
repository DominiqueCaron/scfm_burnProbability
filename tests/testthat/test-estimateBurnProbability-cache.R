test_that("estimateBurnProbability's Cache() correctly hits/misses on pIg/p0/pSpread changes", {
  skip_if_not_installed("SpaDES.core")
  skip_if_not_installed("SpaDES.tools")
  skip_if_not_installed("reproducible")
  skip_if_not_installed("terra")

  studyArea <- LandR::randomStudyArea(size = 10000 * 6.25 * 1000, seed = 3031)
  rasterToMatch <- terra::rast(studyArea, vals = 1, res = c(250, 250)) |>
    terra::mask(studyArea)

  times <- list(start = 0, end = 0)
  modules <- list(
    "scfmIgnition",
    "scfmEscape",
    "scfmSpread",
    "scfm_burnProbability"
  )
  cachePath <- file.path(
    tempdir(),
    paste0("cache-estimateBurnProbability-", sample.int(1e6, 1))
  )
  paths <- list(
    modulePath = c(
      testthat::test_path("..", ".."),
      testthat::test_path("..", "..", "..", "modules")
    ),
    cachePath = cachePath,
    outputPath = file.path(tempdir(), "outputs")
  )
  params <- list(
    scfm_burnProbability = list(nIterations = 100L, .useCache = TRUE)
  )
  objects <- list(studyArea = studyArea, rasterToMatch = rasterToMatch)

  mySim <- SpaDES.core::simInit(
    times = times,
    params = params,
    modules = modules,
    objects = objects,
    paths = paths
  )
  mySim <- SpaDES.core::spades(mySim, debug = FALSE)
  mySim@current <- data.table::data.table(
    eventTime = 0,
    moduleName = "scfm_burnProbability",
    eventType = "estimate",
    eventPriority = 7.51
  )
  fastFn <- mySim$.mods$scfm_burnProbability$estimateBurnProbability

  cacheRdsFiles <- function() {
    list.files(file.path(cachePath, "cacheOutputs"), pattern = "\\.rds$")
  }

  ## First call: cache miss, must compute and save a cache entry.
  mySim$pIg <- 0.01
  t1 <- system.time(out1 <- fastFn(mySim))
  nCacheFiles1 <- length(cacheRdsFiles())
  expect_true(nCacheFiles1 > 0)

  ## Second call with identical inputs: cache hit, no new cache entry, much faster.
  t2 <- system.time(out2 <- fastFn(mySim))
  nCacheFiles2 <- length(cacheRdsFiles())
  expect_equal(nCacheFiles2, nCacheFiles1)
  expect_lt(t2[["elapsed"]], t1[["elapsed"]])
  expect_identical(
    terra::values(out1$burnProbability, mat = FALSE),
    terra::values(out2$burnProbability, mat = FALSE)
  )

  ## Changing pIg: cache miss, new cache entry.
  mySim$pIg <- 0.5
  out3 <- fastFn(mySim)
  nCacheFiles3 <- length(cacheRdsFiles())
  expect_true(nCacheFiles3 > nCacheFiles2)
})
