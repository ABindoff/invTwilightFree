# From topological holonomy to a class-aware estimator

Working note. Companion to `phase_residual_prototype.png` (the geometry) and
`holonomy_class_selection.png` (the estimator skeleton). This is the thread that
could become a second paper; it is not part of the current MEE submission.

## The reframing in one line

Light-level track reconstruction is a section of a circle-bundle over the known
time-torus. Its failure modes are not continuous estimation errors but two
global degrees of freedom:

- Hemisphere class: a discrete, topological **Z₂** (north vs south branch). The
  day-length signal that separates the branches is `D(-φ) = 24 - D(φ)`, so it
  vanishes exactly at the equinox (`D = 12 h`). A local or greedy estimator can
  flip there; this is the classic hemisphere swap.
- Longitude offset: a continuous **gauge** (U(1)-like). A constant clock error
  rigidly translates the whole track in longitude. Because time is known, this
  holonomy collapses to a single estimable constant, removed by one absolute fix.

The estimator that follows the topology is therefore: **select the discrete
basin (class), then slide to the continuous bottom (gauge), then smooth within
the class.** Current methods (FFBS, Viterbi, forward-backward) do the third step
well but handle the first two implicitly, which is why they can sit in the wrong
branch through an equinox.

## What the prototype shows (figures, all seeded)

Geometry (`phase_residual_prototype.png`):
- The diurnal-phase residual is latitude- and time-free; its zero-locus is the
  true meridian (longitude is a pure phase). `sd(phi)` across latitude = 0.
- Day-length sensitivity to latitude is `2.46 h/deg` at solstice and exactly
  `0.0000 h/deg` at the equinox: the latitude fiber pinch is an exact degeneration.
- A southern track and its northern mirror differ by `8.29 h` of day length at
  solstice and `0.08 h` at the equinox branch point.

Estimator skeleton (`holonomy_class_selection.png`), 120-day southern track
across the March equinox:
- Per-day class evidence (log Bayes factor, south vs north) averages `177 nats/day`
  over the season but only `0.8 nats/day` in the +/-5-day equinox window. The
  pinch is where the data go quiet.
- Cumulative class evidence is decisive (`~2.1e4 nats`): the **global** score
  selects the branch even though any single equinox day cannot.
- A `+18 min` clock bias appears as a rigid `+4.5 deg` longitude offset; one
  known fix recovers `4.44 deg` and cuts longitude RMSE from `4.46 -> 0.28 deg`.
- The score landscape over (class) x (longitude offset) is two basins; the south
  basin is lower by the class margin, and its minimum sits at the true offset.
  Discrete choice (vertical separation) then continuous refinement (horizontal).

## Concrete implementation path in invTwilightFree

The package already has the pieces; the missing primitive is the marginal
likelihood.

1. **Expose `logZ` from the grid HMM.** `run_grid_hmm` already runs the forward
   pass; return the log marginal likelihood `logZ = logsumexp(alpha[K-1])`
   alongside the path. This is a few lines in Rust and is generally useful for
   model comparison, not just this.
2. **`class_select()`** (Z₂): run `TwilightFreeGrid` twice, once with
   `hemisphere_prior("S")` and once with `hemisphere_prior("N")` (we already
   ship `hemisphere_prior`), and compare `logZ`. Report the Bayes factor and pick
   the branch. This turns the equinox flip from a local trap into an explicit,
   cheap discrete model choice over two options.
3. **Longitude gauge diagnostic.** The deployment fix already pins longitude; add
   a reported clock-bias estimate (the constant phase residual at the known
   start) so a track-wide longitude error is surfaced as one number rather than
   smeared across knots.
4. **Within-class smoothing** is the existing grid HMM / FFBS, unchanged.

The headline claim for a methods paper: credible-interval coverage and
hemisphere-recovery rate through the equinox, class-aware vs standard smoothing,
across replicate simulated tracks. The class-aware estimator should not flip.

## Honest caveats

- Longitude is a continuous gauge here, not an integer winding number: with a
  known calendar there is no 360-degree lap ambiguity. The only genuinely
  topological (integer) invariant is the hemisphere Z₂. Keep the language precise.
- Real day-length information is noisy and the branches are only weakly separated
  for a window around each equinox; the margin depends on how much non-equinox
  season the deployment covers. Quantifying the required season length is itself
  a useful result.
- SST or bathymetry break the Z₂ directly (independent coordinates the light
  pinch does not touch), which is the multi-sensor counterpart to class selection.
