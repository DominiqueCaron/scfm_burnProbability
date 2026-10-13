# Everything in this file gets sourced during simInit, and all functions and objects
# are put into the simList. To use objects and functions, use sim$xxx.
defineModule(sim, list(
  name = "scfm_burnProbability",
  description = paste(
    "Estimates each pixel's probability of burning in a given year by Monte Carlo:",
    "repeatedly re-draws ignitions, escapes, and fire spread."),
  keywords = c("fire", "probability", "monte carlo", "scfm"),
  authors = c(
    person("Dominique", "Caron", email = "", role = c("aut", "cre"))
  ),
  childModules = character(),
  version = list(scfm_burnProbability = "0.1.0", SpaDES.core = "0.1.0"),
  timeframe = as.POSIXlt(c(NA, NA)),
  timeunit = "year",
  citation = list("citation.bib"),
  documentation = list("README.md", "scfm_burnProbability.Rmd"),
  reqdPkgs = list("data.table", "terra", "SpaDES.core (>=3.0.4)", "reproducible"),
  loadOrder = list(after = c("scfmIgnition", "scfmEscape", "scfmSpread")),
  parameters = bindrows(
    defineParameter("nIterations", "integer", 10^6, 1L, NA,
                    desc = paste("Number of Monte Carlo iterations (independent ignition/escape/spread",
                                 "draws) used to estimate burn probability at each `.runInterval`.")),
    defineParameter(".runInitialTime", "numeric", default = start(sim), NA, NA,
                    desc = "Time of the first burn-probability estimate."),
    defineParameter(".runInterval", "numeric", default = 1, NA, NA,
                    desc = paste("Years between burn-probability estimates. `NA` estimates once, at",
                                 "`.runInitialTime`. Each estimate uses that year's ignition/escape/spread",
                                 "probability surfaces (`pIg`, `p0`, `pSpread`), held fixed across all",
                                 "`nIterations` draws.")),
    defineParameter(".studyAreaName", "character", NA, NA, NA,
                    "Human-readable name for the study area used."),
    defineParameter(".useCache", "logical", TRUE, NA, NA,
                    desc = paste("Cache the Monte Carlo estimate (`runMonteCarlo()`), keyed on `pIg`,",
                                 "`p0`, `pSpread`, `flammableMap`, `fireRegimePolys`, `fireRegimeRas`,",
                                 "`nIterations` and the `neighbours` parameters. Repeated calls with",
                                 "unchanged inputs (e.g. the same year re-run, or an unchanged scenario)",
                                 "return the cached `burnProbability` instead of re-running `nIterations`",
                                 "Monte Carlo draws."))
  ),
  inputObjects = bindrows(
    expectsInput("flammableMap", "SpatRaster", sourceURL = NA,
                 paste("Binary/categorical flammability raster, as used by `scfmEscape` and",
                       "`scfmSpread`. Template for `burnProbability`; also passed straight to",
                       "`spread2()` inside `fastEscape()`/`fastBurn()`.")),
    expectsInput("fireRegimePolys", "SpatVector",
                 "Fire regime polygons with `maxBurnCells` and `PolyID`, from `scfmDataPrep`."),
    expectsInput("fireRegimeRas", "SpatRaster",
                 "Raster of fire regime polygon IDs, from `scfmDataPrep`."),
    expectsInput("pIg", "SpatRaster|numeric",
                 "This year's per-pixel ignition probability surface, from `scfmIgnition`."),
    expectsInput("p0", "SpatRaster|numeric",
                 "This year's per-pixel escape probability surface, from `scfmEscape`."),
    expectsInput("pSpread", "SpatRaster|numeric",
                 "This year's per-pixel spread probability surface, from `scfmSpread`.")
  ),
  outputObjects = bindrows(
    createsOutput("burnProbability", "SpatRaster",
                  desc = paste("Estimated probability each flammable pixel burns this interval:",
                               "fraction of `nIterations` Monte Carlo draws in which the pixel burned.",
                               "`NA` where not flammable."))
  )
))

#' Event dispatcher
#'
#' Events: `init` (schedules the first `estimate`), `estimate` (runs the Monte
#' Carlo loop and reschedules itself every `.runInterval`).
#'
#' @param sim A `simList`.
#' @param eventTime Time of the event.
#' @param eventType Name of the event.
#' @param debug Not used.
#'
#' @return The `simList`, invisibly.
doEvent.scfm_burnProbability <- function(sim, eventTime, eventType, debug = FALSE) {
  moduleName <- current(sim)$moduleName

  switch(
    eventType,
    init = {
      sim <- scheduleEvent(sim, eventTime = P(sim)$.runInitialTime, moduleName, "estimate",
                           eventPriority = 7.51) ## after scfmIgnition/scfmEscape/scfmSpread (7.5)
    },
    estimate = {
      sim <- estimateBurnProbability(sim)

      if (!is.na(P(sim)$.runInterval))
        sim <- scheduleEvent(sim, time(sim) + P(sim)$.runInterval, moduleName, "estimate",
                             eventPriority = 7.51)
    },
    warning(paste("Undefined event type: '", current(sim)[1, "eventType", with = FALSE],
                  "' in module '", current(sim)[1, "moduleName", with = FALSE], "'",
                  sep = ""
    ))
  )

  invisible(sim)
}

#' Estimate this interval's burn probability by Monte Carlo
#'
#' Thin wrapper around `runMonteCarlo()`: extracts this interval's inputs from
#' `sim` and, if `P(sim)$.useCache` is `TRUE`, runs it through `Cache()` so that
#' repeated calls with the same `pIg`, `p0`, `pSpread`, `flammableMap`,
#' `fireRegimePolys`, `fireRegimeRas`, `nIterations` and `neighbours` values
#' return the cached `burnProbability` instead of re-running `nIterations` Monte
#' Carlo draws.
#'
#' @param sim A `simList`.
#'
#' @return The `simList`, invisibly, with `burnProbability`.
estimateBurnProbability <- function(sim) {
  pIg <- sim$pIg
  p0 <- sim$p0
  pSpread <- sim$pSpread
  flammableMap <- sim$flammableMap
  fireRegimePolys <- sim$fireRegimePolys
  fireRegimeRas <- sim$fireRegimeRas
  nIter <- P(sim)$nIterations
  neighboursEscape <- P(sim, "neighbours", "scfmEscape")
  neighboursSpread <- P(sim, "neighbours", "scfmSpread")

  sim$burnProbability <- if (isTRUE(P(sim)$.useCache)) {
    runMonteCarlo(pIg = pIg, p0 = p0, pSpread = pSpread, flammableMap = flammableMap,
                  fireRegimePolys = fireRegimePolys, fireRegimeRas = fireRegimeRas, nIter = nIter,
                  neighboursEscape = neighboursEscape, neighboursSpread = neighboursSpread) |>
      Cache(cachePath = cachePath(sim), userTags = c(currentModule(sim), "runMonteCarlo"))
  } else {
    runMonteCarlo(pIg = pIg, p0 = p0, pSpread = pSpread, flammableMap = flammableMap,
                  fireRegimePolys = fireRegimePolys, fireRegimeRas = fireRegimeRas, nIter = nIter,
                  neighboursEscape = neighboursEscape, neighboursSpread = neighboursSpread)
  }

  invisible(sim)
}

#' Run the Monte Carlo burn-probability estimate
#'
#' Runs `nIter` independent draws of ignitions, escapes and fire spread and
#' tallies how many iterations burn each pixel. 
#'
#' For speed, each draw is simulated directly with `fastIgnite()`/`fastEscape()`/
#' `fastBurn()` rather than by calling `scfmIgnition::Ignite()`/
#' `scfmEscape::Escape()`/`scfmSpread::Burnemup()` through `sim$.mods` (as an
#' earlier version of this function did). Profiling showed those modules' own
#' event functions spend the large majority of their time on `terra` `SpatRaster`
#' indexed reads/writes and `scfmSpread::Burnemup()`'s reporting bookkeeping
#' (`countBurnedPixelsInSAR()`, `burnMap`/`timeSinceFire`/`burnSummary` updates,
#' optional plotting) -- all necessary for a real simulation year, but not for
#' estimating burn probability, which only needs the set of burned pixels per
#' draw. The `fastIgnite()`/`fastEscape()`/`fastBurn()` helpers below reproduce
#' the same stochastic draws (same `calcIgnitions()` logic, same `spread2()` calls
#' with the same arguments) using plain vectors instead of `sim$`/`SpatRaster`
#' access, and `burnCount` is accumulated as a plain integer vector and converted
#' to a `SpatRaster` only once, after the loop. **These helpers duplicate core
#' logic from `scfmIgnition::Ignite()`/`calcIgnitions()`, `scfmEscape::Escape()`
#' and `scfmSpread::Burnemup()`; if those modules' fire-process logic changes,
#' this function must be updated to match.**
#'
#' @param pIg `SpatRaster` or numeric; this year's per-pixel ignition probability surface.
#' @param p0 `SpatRaster` or numeric; this year's per-pixel escape probability surface.
#' @param pSpread `SpatRaster` or numeric; this year's per-pixel spread probability surface.
#' @param flammableMap `SpatRaster`; binary/categorical flammability raster. Template
#'   for the returned `burnProbability`, and passed straight to `spread2()`.
#' @param fireRegimePolys `SpatVector`; fire regime polygons with `maxBurnCells` and `PolyID`.
#' @param fireRegimeRas `SpatRaster`; raster of fire regime polygon IDs.
#' @param nIter integer; number of Monte Carlo iterations.
#' @param neighboursEscape integer; `scfmEscape`'s `neighbours` parameter.
#' @param neighboursSpread integer; `scfmSpread`'s `neighbours` parameter.
#'
#' @return A `SpatRaster`: estimated probability each flammable pixel burns this
#'   interval (`NA` where not flammable).
runMonteCarlo <- function(pIg, p0, pSpread, flammableMap, fireRegimePolys, fireRegimeRas,
                         nIter, neighboursEscape, neighboursSpread) {
  flamVals <- as.vector(flammableMap)
  nCell <- terra::ncell(flammableMap)

  fireRegimeVals <- as.vector(fireRegimeRas)
  fireRegimePixels <- which(!is.na(fireRegimeVals))
  pIgVals <- if (inherits(pIg, "SpatRaster")) as.vector(pIg) else pIg
  p0Vals <- if (inherits(p0, "SpatRaster")) as.vector(p0) else p0
  ## pSpread is left as-is: spread2() accepts either a SpatRaster or a scalar directly

  maxSizesByPoly <- if (!is.null(fireRegimePolys$maxBurnCells)) {
    fireRegimePolys$maxBurnCells[match(fireRegimeVals, fireRegimePolys$PolyID)]
  } else {
    NULL ## no per-polygon cap: fastEscape() falls back to NA_real_ (no cap)
  }

  burnCount <- integer(nCell)

  for (i in seq_len(nIter)) {
    ignitionLoci <- fastIgnite(fireRegimePixels, pIgVals, nCell)
    spreadState <- fastEscape(flammableMap, ignitionLoci, p0Vals, neighboursEscape,
                              maxSizesByPoly)
    burnedPixels <- fastBurn(flammableMap, spreadState, pSpread, neighboursSpread)

    if (length(burnedPixels) > 0) burnCount[burnedPixels] <- burnCount[burnedPixels] + 1L
  }

  burnProbability <- rast(flammableMap)
  burnProbability[] <- burnCount / nIter
  burnProbability[!(flamVals %in% 1)] <- NA ## flamVals may hold NaN (not NA) where not flammable
  burnProbability
}

#' Draw one Monte Carlo iteration's ignitions
#'
#' Mirrors `scfmIgnition::calcIgnitions()` (same `runif()` draw against `pIg`),
#' operating on a plain vector of candidate pixel indices and probabilities
#' instead of re-deriving them from `sim$` each call. Does not reproduce
#' `scfmIgnition::Ignite()`'s `SpaDES.tools:::resample()` permutation step, which
#' only randomizes draw order for cosmetic/display reasons and does not affect
#' which pixels ignite.
#'
#' @param fireRegimePixels integer; candidate pixel indices (non-`NA` in
#'   `fireRegimeRas`), precomputed once outside the Monte Carlo loop.
#' @param pIgVals numeric; either a per-pixel ignition-probability vector
#'   (same length as the full raster, indexed like `fireRegimePixels`) or a
#'   single scalar probability.
#' @param nCell integer; total number of cells in the template raster (used only
#'   when `pIgVals` is a per-pixel vector).
#'
#' @return Integer vector of ignition pixel indices.
fastIgnite <- function(fireRegimePixels, pIgVals, nCell) {
  if (length(pIgVals) > 1) {
    which(runif(nCell) < pIgVals) ## NAs in pIgVals are omitted, as in calcIgnitions()
  } else {
    fireRegimePixels[which(runif(length(fireRegimePixels)) <= pIgVals)]
  }
}

#' Draw one Monte Carlo iteration's escapes
#'
#' Mirrors `scfmEscape::Escape()`: a single-iteration `spread2()` call from this
#' draw's ignition loci, with the same per-fire size cap logic (`maxBurnCells`,
#' looked up per ignition pixel's fire-regime polygon).
#'
#' @param flammableMap `SpatRaster`; passed straight to `spread2()`.
#' @param ignitionLoci integer; this draw's ignition pixel indices, from `fastIgnite()`.
#' @param p0Vals `SpatRaster` or numeric; escape probability, passed straight to `spread2()`.
#' @param neighbours integer; `scfmEscape`'s `neighbours` parameter.
#' @param maxSizesByPoly numeric or `NULL`; per-pixel `maxBurnCells` lookup
#'   (indexed like the full raster), precomputed once outside the loop, or `NULL`
#'   if `fireRegimePolys` has no `maxBurnCells` column.
#'
#' @return A `data.table`, as returned by `spread2()` (`NULL` if no ignitions).
fastEscape <- function(flammableMap, ignitionLoci, p0Vals, neighbours, maxSizesByPoly) {
  if (length(ignitionLoci) == 0) return(NULL)

  maxSizes <- if (is.null(maxSizesByPoly)) NA_real_ else maxSizesByPoly[ignitionLoci]

  SpaDES.tools::spread2(landscape = flammableMap, start = ignitionLoci, iterations = 1,
                        spreadProb = p0Vals, directions = neighbours, maxSize = maxSizes,
                        asRaster = FALSE)
}

#' Draw one Monte Carlo iteration's fire spread
#'
#' Mirrors `scfmSpread::Burnemup()`'s core `spread2()` call, skipping the
#' reporting/bookkeeping steps (`countBurnedPixelsInSAR()`,
#' `rstCurrentBurn`/`burnMap`/`timeSinceFire`/`burnSummary` updates, plotting)
#' that are not needed for the burn-probability estimate. Matches `Burnemup()`'s
#' own `maxSizes <- NA_real_` fallback when `fireRegimePolys$maxBurnCells` is
#' set, since each fire's cap already travels with it via `spreadState` by the
#' time this step runs.
#'
#' @param flammableMap `SpatRaster`; passed straight to `spread2()`.
#' @param spreadState `data.table` or `NULL`; this draw's escape result, from `fastEscape()`.
#' @param pSpread `SpatRaster` or numeric; spread probability, passed straight to `spread2()`.
#' @param neighbours integer; `scfmSpread`'s `neighbours` parameter.
#'
#' @return Integer vector of burned pixel indices (empty if nothing burned).
fastBurn <- function(flammableMap, spreadState, pSpread, neighbours) {
  if (is.null(spreadState) || nrow(spreadState) == 0) return(integer(0))

  if (nrow(spreadState[spreadState$state == "activeSource", ]) > 0) {
    burnDT <- SpaDES.tools::spread2(flammableMap, start = spreadState, spreadProb = pSpread,
                                    directions = neighbours, maxSize = NA_real_, asRaster = FALSE)
  } else {
    ## No fire escaped: each ignition burned its own pixel (as in Burnemup()).
    burnDT <- spreadState
  }

  burnDT$pixels
}

#' Supply default inputs
#'
#' There are none.
#'
#' @param sim A `simList`.
#'
#' @return The `simList`, invisibly.
.inputObjects <- function(sim) {
  return(invisible(sim))
}
