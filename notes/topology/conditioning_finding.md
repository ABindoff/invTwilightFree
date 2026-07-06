# Conditioning: hypothesis, and what the demonstration actually showed

Working note. Figure: `conditioning_demo.png`. Primitives: `R/dual_time.R`
(`eval_logpt_loc`, `calibrate_clock`), tests in `tests/testthat/test-dual_time.R`.

## The hypothesis (and where it was wrong)

The strategic claim was: light geolocation is strongly anisotropic, latitude is
the soft direction, the condition number blows up at the equinox, so reduce the
problem to a 1D latitude search along a longitude line of position.

Running the information geometry refuted the strong form of this for the
**continuous-amplitude** likelihood the package uses:

- The per-day marginal Fisher information for latitude does NOT collapse at the
  equinox. It is `17.97` at solstice and `15.67` at equinox. Midday sun elevation
  (`max elev = 90 - |phi - decl|`) keeps constraining latitude magnitude even when
  day length is 12 h everywhere. The classic "latitude is hard near the equinox"
  is a property of THRESHOLD methods, which discard amplitude; the continuous
  model does not share it. This is a real advantage worth stating in the paper.
- There is no fixed soft eigen-direction to project out, so the 1D-reduction
  method does not follow. The condition number is mild (median ~2) and does not
  diverge at the equinox.

## What is real and exploitable

1. Clock / longitude calibration from the known endpoints. This held up exactly.
   The dual likelihood at a known fix is sharply peaked in time, so the two
   endpoints (location AND time known) pin a linear clock drift, and removing it
   takes interior longitude RMSE from `4.34 deg` to `~3e-15 deg`. This is the
   genuine payoff of knowing `deployed.at` and `retrieved.at`, and it is what the
   two new primitives deliver.
2. The residual equinox degeneracy is the discrete hemisphere SIGN, a Z2
   reflection. At the exact equinox instant `decl = 0`, the predicted light is
   identical for latitude `+phi` and `-phi`, so the latitude likelihood is a
   symmetric pair of modes at the true magnitude. The mode is exactly degenerate
   only at that instant; on the equinox day with clean data the wrong hemisphere
   is suppressed by ~26 nats, and realistic light noise washes that out over a
   window around the equinox, which is exactly when hemisphere flips occur. This
   is handled by the discrete machinery already built: `hemisphere_prior` and the
   `log_z` Bayes factor (`class_select`).

## Consequence for the method

There is no new 1D-latitude engine. The two handles that matter are both already
in hand:

- the clock-calibration primitives (`eval_logpt_loc`, `calibrate_clock`), which
  use the endpoints to make longitude exact across the interior; and
- hemisphere class selection (`log_z` / `hemisphere_prior`), which resolves the
  one genuine degeneracy, the sign, over the season.

So "best exploit the better conditioning" resolves to: calibrate the clock from
the endpoints (continuous, gauge), and select the hemisphere by evidence
(discrete, topological). The continuous model is otherwise well determined in
both coordinates, which is itself a selling point over threshold methods.

## Honest caveats

- The Fisher information here is the Gauss-Newton information of the expected
  light; it ignores the spike-kink and assumes Gaussian noise. The qualitative
  conclusion (latitude is well informed by amplitude; no equinox collapse) is
  robust, but exact numbers depend on noise model and sampling cadence.
- The hemisphere modes are exactly degenerate only at the equinox instant; the
  practical failure is a noise-driven window, not a point.
