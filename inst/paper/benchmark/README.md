# Benchmark reproduction

Everything the benchmark tables rest on is simulated, so no external data is
needed. Run from the package root with R 4.6.0, which is the installation
carrying SGAT, FLightR, TwGeos and raster.

## Table 2 (seal scenarios)

The 40-day track and its three noise scenarios are regenerated from fixed
seeds, 999 for the movement and 555 for the light:

```sh
Rscript inst/paper/benchmark/01_simulate_seal_track.R   # -> simulated_seal_track_crw.rds
Rscript inst/paper/benchmark/02_generate_seal_light.R   # -> simulated_seal_light_scenarios.rds
```

The second script writes `light_ideal`, `light_shaded` and `light_alan`, which
are the Cloudy, Shaded and ALAN columns of Table 2. The intermediate `.rds`
files are not tracked: they are byte-identical on any machine from the seeds
above, and regenerating them takes seconds.

## Table 1 (geometric stress test)

The 180-day path is regenerated inside the driver itself (`simulate_geom()`),
seed 42 for the ideal scenario and 123 for the adverse one.

## Both tables

```sh
Rscript inst/paper/rerun_benchmarks.R
```

This measures every method back to back in a single process, with all thread
counts pinned to 1 and a fixed arithmetic probe run before and after so that
load drift is visible. It writes `inst/paper/benchmark_rerun.csv`, one row per
fit, with the Monte Carlo effort recorded next to each timing, and
`inst/paper/fig1_fits.rds` for Figure 1.

Note that the earlier drivers under `scratch/` recorded comparator runtimes as
literal constants (`time = c(elapsed = 131*60)`) rather than measuring them.
`rerun_benchmarks.R` supersedes them; the timings in Tables 1 and 2 come from
its output.

## Twilight settings

SGAT and FLightR need twilights, which `TwGeos::findTwilights` supplies. Under
noise it returns roughly ten times as many candidates as under cloud, and these
must be thinned by a minimum separation. The driver uses a threshold of 10 and
a 4 h separation for the seal scenarios, which is the only combination tested
that lets SGAT sample on all three; at a 2 h separation its sampler terminates
with an internal error on the shaded and ALAN records. Section 3.2 of the
manuscript reports the sensitivity.
