# scfm_burnProbability

Estimates each flammable pixel's probability of burning, by Monte Carlo, mirroring
the draw logic of the scfm ignition, escape, and spread modules
(`scfmIgnition`, `scfmEscape`, `scfmSpread`) via fast internal helpers rather than
calling those modules directly. Parallels `fireSense_burnProbability`. See
`scfm_burnProbability.Rmd` for module documentation.

<!-- Drafted with assistance from Claude (Posit Assistant). -->
