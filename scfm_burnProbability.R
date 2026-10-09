# Everything in this file gets sourced during simInit, and all functions and objects
# are put into the simList. To use objects and functions, use sim$xxx.
defineModule(sim, list(
  name = "scfm_burnProbability",
  description = paste(
    "Estimates each pixel's probability of burning in a given year by Monte Carlo:",
    "repeatedly re-draws ignitions, escapes, and fire spread (via",
    "scfmIgnition, scfmEscape and scfmSpread) without touching the",
    "simulation's realized fire outputs, and tallies how often each pixel burns."),
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
  reqdPkgs = list("data.table", "terra", "SpaDES.core (>=3.0.4)"),
  loadOrder = list(after = c("scfmIgnition", "scfmEscape", "scfmSpread")),
  parameters = bindrows(
    defineParameter("nIterations", "integer", 10000L, 1L, NA,
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
                    "Human-readable name for the study area used.")
  ),
  inputObjects = bindrows(
    expectsInput("flammableMap", "SpatRaster", sourceURL = NA,
                 "Binary/categorical flammability raster, as used by `scfmEscape` and `scfmSpread`. Template for `burnProbability`."),
    expectsInput("rasterToMatch", "SpatRaster", sourceURL = NA,
                 "Template raster for the study area (as used by the scfm modules)."),
    expectsInput("fireRegimePolys", "SpatVector",
                 paste("Fire regime polygons with `pIgnition`, `p0`, `pSpread`, `maxBurnCells`,",
                       "`cellSize` and `PolyID`, from `scfmDataPrep`. Held fixed across Monte Carlo draws.")),
    expectsInput("fireRegimeRas", "SpatRaster",
                 "Raster of fire regime polygon IDs, from `scfmDataPrep`. Held fixed across Monte Carlo draws."),
    expectsInput("pIg", "SpatRaster|numeric",
                 "This year's per-pixel ignition probability surface, from `scfmIgnition`. Held fixed across Monte Carlo draws."),
    expectsInput("p0", "SpatRaster|numeric",
                 "This year's per-pixel escape probability surface, from `scfmEscape`. Held fixed across Monte Carlo draws."),
    expectsInput("pSpread", "SpatRaster|numeric",
                 "This year's per-pixel spread probability surface, from `scfmSpread`. Held fixed across Monte Carlo draws."),
    expectsInput("ignitionLoci", "numeric",
                 paste("This year's *realized* ignition pixels, from `scfmIgnition`.",
                       "Snapshotted and restored; never left modified by the Monte Carlo loop.")),
    expectsInput("spreadState", "data.table",
                 paste("This year's *realized* escaped-fire start pixels, from `scfmEscape`.",
                       "Snapshotted and restored; never left modified by the Monte Carlo loop.")),
    expectsInput("burnDT", "data.table",
                 "This year's *realized* burned-pixel table, from `scfmSpread`. Snapshotted and restored."),
    expectsInput("burnMap", "SpatRaster",
                 paste("Cumulative burn map, from `scfmSpread`.",
                       "Snapshotted and restored, never left modified by the Monte Carlo loop.")),
    expectsInput("rstCurrentBurn", "SpatRaster",
                 "This year's *realized* current-year burn indicator, from `scfmSpread`. Snapshotted and restored."),
    expectsInput("burnSummary", "data.table",
                 paste("Cumulative per-fire summary, from `scfmSpread`. Snapshotted and restored",
                       "(each Monte Carlo draw would otherwise append a spurious row).")),
    expectsInput("timeSinceFire", "SpatRaster",
                 paste("Time-since-fire raster, from `scfmSpread`; updated every year including years",
                       "with no fire. Snapshotted and restored."))
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
#' Snapshots every `sim$` object that `scfmIgnition::Ignite()`, `scfmEscape::Escape()`
#' and `scfmSpread::Burnemup()` write to, then runs `P(sim)$nIterations` independent
#' draws of ignitions, escapes and fire spread (re-calling those modules' own
#' functions, via `sim$.mods`, so the Monte Carlo draws use exactly the same
#' stochastic logic as a real simulation year), tallying how many iterations burn
#' each pixel. The snapshotted objects are restored before returning, so this module
#' never alters the "realized" simulation outputs other modules (e.g. vegetation or
#' carbon accounting) depend on -- it only adds `burnProbability`.
#'
#' Unlike `fireSense_burn`, scfm has no year-level random effect shared across that
#' year's fires: the only stochasticity is the per-pixel ignition draw
#' (`calcIgnitions()`, inside `Ignite()`) and the `spread2()` calls inside `Escape()`
#' and `Burnemup()`. The probability surfaces `pIg`, `p0` and `pSpread` are computed
#' once per year upstream and are *not* recomputed each iteration: they do not change
#' between draws within one interval, only the ignition/escape/spread realization does.
#'
#' `scfmSpread::Burnemup()` reads `mod$pixInSAR` from its own module environment, so
#' it must be called through `sim$.mods$scfmSpread$Burnemup()`, not reimplemented here.
#'
#' @param sim A `simList`.
#'
#' @return The `simList`, invisibly, with `burnProbability`.
estimateBurnProbability <- function(sim) {
  nIter <- P(sim)$nIterations

  snapshotNames <- c("ignitionLoci", "spreadState", "burnDT", "burnMap", "rstCurrentBurn",
                     "burnSummary", "timeSinceFire")
  snapshot <- snapshotSimObjects(sim, snapshotNames)

  burnCount <- rast(sim$flammableMap)
  burnCount[] <- 0
  flamVals <- as.vector(sim$flammableMap)
  burnCount[flamVals == 0] <- NA

  for (i in seq_len(nIter)) {
    sim <- sim$.mods$scfmIgnition$Ignite(sim)
    sim <- sim$.mods$scfmEscape$Escape(sim)
    sim <- sim$.mods$scfmSpread$Burnemup(sim)

    burned <- !is.na(terra::values(sim$rstCurrentBurn, mat = FALSE)) &
      terra::values(sim$rstCurrentBurn, mat = FALSE) == 1
    burnCount[which(burned)] <- burnCount[which(burned)] + 1
  }

  sim <- restoreSimObjects(sim, snapshot)

  sim$burnProbability <- burnCount / nIter

  invisible(sim)
}

#' Snapshot a set of `sim$` objects before they are overwritten
#'
#' `SpatRaster`s are deep-copied (`terra::deepcopy()`): otherwise they alias the same
#' underlying pointer, and `restoreSimObjects()` would "restore" a reference to data
#' the Monte Carlo loop has already overwritten. `data.table`s are copied with
#' `data.table::copy()` for the same reason (they are modified by reference).
#'
#' @param sim A `simList`.
#' @param names character; the `sim$` object names to snapshot.
#'
#' @return A named list holding each object's copy (or `NULL` if absent).
snapshotSimObjects <- function(sim, names) {
  setNames(lapply(names, function(n) {
    obj <- sim[[n]]
    if (is.null(obj)) return(NULL)
    if (inherits(obj, "SpatRaster")) return(terra::deepcopy(obj))
    if (is.data.table(obj)) return(data.table::copy(obj))
    obj
  }), names)
}

#' Restore `sim$` objects from a snapshot
#'
#' @param sim A `simList`.
#' @param snapshot A named list, as returned by `snapshotSimObjects()`.
#'
#' @return The `simList`, invisibly, with each snapshotted object restored.
restoreSimObjects <- function(sim, snapshot) {
  for (n in names(snapshot)) sim[[n]] <- snapshot[[n]]
  invisible(sim)
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
