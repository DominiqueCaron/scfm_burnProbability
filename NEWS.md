# scfm_burnProbability 0.1.0 (09 October 2026)

- initial module version
- rewrote the Monte Carlo loop to use fast internal helpers (`fastIgnite()`,
  `fastEscape()`, `fastBurn()`) that mirror `scfmIgnition`/`scfmEscape`/
  `scfmSpread`'s stochastic draws on plain vectors/`data.table`s, instead of
  calling those modules' own event functions through `sim$.mods`; profiling
  showed this removes most of the per-iteration `terra`/bookkeeping overhead
  (~50-60x faster per iteration in testing)
- trimmed `inputObjects` to only the objects actually read
  (`flammableMap`, `fireRegimePolys`, `fireRegimeRas`, `pIg`, `p0`, `pSpread`)
- added `.useCache` parameter (default `TRUE`): caches the Monte Carlo estimate
  via `reproducible::Cache()`, keyed on `pIg`, `p0`, `pSpread`, `flammableMap`,
  `fireRegimePolys`, `fireRegimeRas`, `nIterations`, and `scfmEscape`'s/
  `scfmSpread`'s `neighbours` parameters
