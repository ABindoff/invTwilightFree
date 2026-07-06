"""
From holonomy to a class-aware estimator: a prototype.

Two global degrees of freedom organise the failure modes of light geolocation:

  * Hemisphere class  -- a DISCRETE, topological Z2 (north vs south branch).
    The day-length signal that separates the branches vanishes at the equinox
    (the pinch), so any local estimator can flip there. A GLOBAL class score
    integrated over the season selects the branch with a clear margin.

  * Longitude offset  -- a CONTINUOUS gauge (U(1)-like). A constant clock error
    rigidly translates the whole track in longitude: a holonomy that the known
    clock exposes as a single estimable constant, removed by one absolute fix.

The estimator skeleton is therefore: pick the discrete basin (class), then slide
to the continuous bottom (gauge). Panel D draws exactly that score landscape.
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

rng = np.random.default_rng(7)
J2000 = 946728000.0

def decl_of(t):
    d = (t - J2000) / 86400.0
    g = np.radians(357.529 + 0.98560028 * d)
    q = 280.459 + 0.98564736 * d
    Lsun = np.radians(q + 1.915 * np.sin(g) + 0.020 * np.sin(2 * g))
    e = np.radians(23.439 - 3.6e-7 * d)
    return np.arcsin(np.sin(e) * np.sin(Lsun))

def day_length_h(lat_deg, decl_rad):
    cosH0 = -np.tan(np.radians(lat_deg)) * np.tan(decl_rad)
    H0 = np.arccos(np.clip(cosH0, -1, 1))
    return 24.0 * H0 / np.pi

# ---- simulate a southern track across the March equinox --------------------
ndays = 120
t0 = (np.datetime64("2024-02-01") - np.datetime64("1970-01-01")) / np.timedelta64(1, "s")
t = float(t0) + np.arange(ndays) * 86400 + 43200
decl = decl_of(t)
true_lat = -50 + 6 * np.sin(2 * np.pi * np.arange(ndays) / 120)   # stays southern
true_lon = 150 + 12 * np.sin(2 * np.pi * np.arange(ndays) / 80)

sigma_D = 0.25                                   # day-length obs noise (hours)
D_obs = day_length_h(true_lat, decl) + rng.normal(0, sigma_D, ndays)

# equinox index (declination sign change)
eqx = int(np.where(np.diff(np.sign(decl)) != 0)[0][0])

# ---- hemisphere class log-likelihoods (south = true sign, north = mirror) --
D_south = day_length_h(true_lat, decl)
D_north = day_length_h(-true_lat, decl)          # = 24 - D_south exactly
LL_south = -0.5 * ((D_obs - D_south) / sigma_D) ** 2
LL_north = -0.5 * ((D_obs - D_north) / sigma_D) ** 2
dLL = LL_south - LL_north                          # per-day log Bayes factor (S vs N)
cumLL = np.cumsum(dLL)

# ---- longitude clock holonomy ---------------------------------------------
clock_bias_min = 18.0
bias_deg = 15.0 * (clock_bias_min / 60.0)          # deg of longitude per clock error
lon_biased = true_lon + bias_deg + rng.normal(0, 0.3, ndays)
# Known deployment fix at day 0 identifies the rigid offset (the gauge):
est_offset = np.mean((lon_biased - true_lon)[:3])  # use first 3 known days
lon_calibrated = lon_biased - est_offset
rmse_before = np.sqrt(np.mean((lon_biased - true_lon) ** 2))
rmse_after = np.sqrt(np.mean((lon_calibrated - true_lon) ** 2))

# ---- score landscape: class (Z2) x longitude offset (U(1)) -----------------
offsets = np.linspace(-10, 12, 221)
sigma_lon = 1.0
lon_term = np.array([0.5 * np.sum(((lon_biased - (true_lon + d)) / sigma_lon) ** 2) for d in offsets])
nll_south = (-np.sum(LL_south)) + lon_term
nll_north = (-np.sum(LL_north)) + lon_term
gmin_off = offsets[np.argmin(nll_south)]

# ---- figure ----------------------------------------------------------------
fig = plt.figure(figsize=(13.5, 9.5))
gs = GridSpec(2, 2, figure=fig, hspace=0.30, wspace=0.22)

axA = fig.add_subplot(gs[0, 0])
axA.axhline(0, color="grey", lw=0.8)
axA.bar(np.arange(ndays), dLL, color=np.where(dLL >= 0, "#1f77b4", "#d62728"), width=1)
axA.axvline(eqx, color="k", ls=":"); axA.annotate("equinox\n(pinch)", (eqx, axA.get_ylim()[1]*0.7), fontsize=9, ha="center")
axA.set_title("A. Per-day class evidence (south - north)\nsignal vanishes at the equinox", fontsize=10)
axA.set_xlabel("day"); axA.set_ylabel("log Bayes factor / day")

axB = fig.add_subplot(gs[0, 1])
axB.plot(cumLL, color="#1f77b4", lw=2)
axB.axhline(0, color="grey", lw=0.8); axB.axvline(eqx, color="k", ls=":")
axB.fill_between(range(max(0,eqx-7), min(ndays,eqx+8)), *axB.get_ylim(), color="0.85", zorder=0)
axB.annotate(f"final margin = {cumLL[-1]:.0f} nats\n(south selected)", (ndays*0.45, cumLL[-1]*0.5), fontsize=9)
axB.set_title("B. Cumulative class evidence\nglobal score selects the branch despite the pinch", fontsize=10)
axB.set_xlabel("day"); axB.set_ylabel("cumulative log Bayes factor")

axC = fig.add_subplot(gs[1, 0])
axC.plot(true_lon, color="k", lw=2, label="true")
axC.plot(lon_biased, color="#d62728", lw=1.5, ls="--", label=f"clock bias {clock_bias_min:.0f} min (+{bias_deg:.1f} deg)")
axC.plot(lon_calibrated, color="#2ca02c", lw=1.5, label="calibrated (one fix)")
axC.set_title(f"C. Longitude clock holonomy: a rigid offset\nRMSE {rmse_before:.2f} deg -> {rmse_after:.2f} deg after one fix", fontsize=10)
axC.set_xlabel("day"); axC.set_ylabel("longitude (deg)"); axC.legend(fontsize=8, loc="upper right")

axD = fig.add_subplot(gs[1, 1])
axD.plot(offsets, nll_south, color="#1f77b4", lw=2, label="south class")
axD.plot(offsets, nll_north, color="#d62728", lw=2, ls="--", label="north class")
axD.plot(gmin_off, np.min(nll_south), "*", color="#1f77b4", ms=14)
axD.axvline(bias_deg, color="grey", ls=":"); axD.annotate("true offset", (bias_deg, axD.get_ylim()[1]*0.9), fontsize=8, rotation=90, va="top")
axD.set_title("D. Score over class (Z2) x longitude offset (gauge)\npick the lower basin, then slide to its bottom", fontsize=10)
axD.set_xlabel("longitude offset (deg)"); axD.set_ylabel("negative log-likelihood"); axD.legend(fontsize=8)

fig.suptitle("A class-aware estimator: discrete hemisphere selection (topological) + continuous longitude calibration (gauge)",
             fontsize=12, y=0.99)
fig.savefig("holonomy_class_selection.png", dpi=140, bbox_inches="tight")
print("saved holonomy_class_selection.png")

# ---- diagnostics -----------------------------------------------------------
win = slice(max(0, eqx-5), min(ndays, eqx+6))
print(f"per-day class evidence: mean over season = {np.mean(dLL):.3f} nats, "
      f"mean in +/-5d equinox window = {np.mean(dLL[win]):.3f} nats")
print(f"cumulative class margin (south vs north) = {cumLL[-1]:.1f} nats "
      f"(Bayes factor ~ exp({cumLL[-1]:.0f}); decisive)")
print(f"clock holonomy: injected +{bias_deg:.2f} deg, recovered {est_offset:.2f} deg from a known fix")
print(f"longitude RMSE: {rmse_before:.3f} deg before, {rmse_after:.3f} deg after calibration")
print(f"score landscape: south NLL min = {np.min(nll_south):.1f}, north NLL min = {np.min(nll_north):.1f}, "
      f"gap = {np.min(nll_north)-np.min(nll_south):.1f}")
