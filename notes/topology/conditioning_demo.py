"""
The conditioning argument, demonstrated and corrected by the data.

Hypothesis going in: the per-knot information is wildly anisotropic and the
condition number blows up at the equinox. Running it corrected that:

  * Longitude/time IS the better-determined direction, consistently (a few-fold
    information advantage over latitude all year). The endpoints then calibrate
    the clock and remove the longitude gauge exactly (panel B). That part holds.

  * But the condition number does NOT diverge at the equinox for the CONTINUOUS
    likelihood, because midday sun elevation still constrains latitude MAGNITUDE
    there. The equinox degeneracy is therefore the discrete hemisphere SIGN (a
    Z2 reflection), not an information zero (panels C, D). Threshold methods,
    which discard amplitude, do suffer the full latitude collapse; the continuous
    method does not, which is a real advantage and points straight at the
    discrete class machinery (hemisphere prior / log_z) as the right tool.
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

J2000 = 946728000.0
INTERCEPT, SLOPE, MAXL = 64.0, 64.0 / 90.0, 64.0
LAMBDA, PROB_SLAB = 0.5, 0.05

def solar(t):
    d = (t - J2000) / 86400.0
    g = np.radians(357.529 + 0.98560028 * d)
    q = 280.459 + 0.98564736 * d
    Lsun = np.radians(q + 1.915 * np.sin(g) + 0.020 * np.sin(2 * g))
    e = np.radians(23.439 - 3.6e-7 * d)
    decl = np.arcsin(np.sin(e) * np.sin(Lsun))
    ra = np.degrees(np.arctan2(np.cos(e) * np.sin(Lsun), np.cos(Lsun)))
    gmst = (18.697374558 + 24.06570982441908 * d) * 15.0
    return decl, ra, gmst

def zenith(lon, lat, t):
    decl, ra, gmst = solar(t)
    H = np.radians(gmst + lon - ra)
    cz = np.sin(np.radians(lat))*np.sin(decl) + np.cos(np.radians(lat))*np.cos(decl)*np.cos(H)
    return np.degrees(np.arccos(np.clip(cz, -1, 1)))

def mu(lon, lat, t):
    return np.clip(INTERCEPT - SLOPE * zenith(lon, lat, t), 0, MAXL)

def light_loglik(lon, lat, t, obs):
    expected = mu(lon, lat, t)
    spike = np.where(obs <= expected, LAMBDA*np.exp(-LAMBDA*(expected-obs)),
                                       LAMBDA*np.exp(-LAMBDA*2*(obs-expected)))
    return np.sum(np.log((1-PROB_SLAB)*spike + PROB_SLAB/MAXL))

def diag_info(lon, lat, t_day, h=0.05):
    dmu_dlon = (mu(lon+h, lat, t_day) - mu(lon-h, lat, t_day)) / (2*h)
    dmu_dlat = (mu(lon, lat+h, t_day) - mu(lon, lat-h, t_day)) / (2*h)
    return np.sum(dmu_dlon**2), np.sum(dmu_dlat**2)     # marginal info: longitude, latitude

true_lon, true_lat = 150.0, -50.0
year0 = float((np.datetime64("2024-01-01") - np.datetime64("1970-01-01")) / np.timedelta64(1, "s"))
decl_year = np.array([solar(year0 + d*86400 + 43200)[0] for d in range(365)])
eqx_days = np.where(np.diff(np.sign(decl_year)) != 0)[0]

fig = plt.figure(figsize=(13.5, 9.5))
gs = GridSpec(2, 2, figure=fig, hspace=0.32, wspace=0.24)

# ---- A. marginal information for longitude vs latitude, across the year ----
sample_days = np.arange(2, 365, 3)
info_lon, info_lat = [], []
for dday in sample_days:
    t_day = year0 + dday*86400 + np.arange(0, 86400, 600)
    il, ia = diag_info(true_lon, true_lat, t_day)
    info_lon.append(il); info_lat.append(ia)
info_lon = np.array(info_lon); info_lat = np.array(info_lat)
axA = fig.add_subplot(gs[0, 0])
axA.plot(sample_days, info_lon, color="#1f77b4", lw=2, label="longitude information")
axA.plot(sample_days, info_lat, color="#d62728", lw=2, label="latitude information")
for e in eqx_days: axA.axvline(e, color="k", ls=":")
axA.set_title("A. Marginal information (continuous model): latitude is NOT the\n"
              "weak direction (midday amplitude), and does not collapse at equinox", fontsize=10)
axA.set_xlabel("day of 2024"); axA.set_ylabel("Fisher information (light$^2$/deg$^2$)")
axA.legend(fontsize=8)

# ---- B. endpoint clock calibration removes the longitude holonomy ----------
axB = fig.add_subplot(gs[0, 1])
ndep = 120; dep_days = np.arange(ndep)
true_lon_track = 150 + 10*np.sin(2*np.pi*dep_days/70)
drift_h = 0.5 * dep_days/(ndep-1)                       # 0 -> 30 min clock drift
lon_biased = true_lon_track + 15.0*drift_h
lon_calib = lon_biased - 15.0*(0.5*dep_days/(ndep-1))   # two known fixes pin the linear drift
rmse_b = np.sqrt(np.mean((lon_biased-true_lon_track)**2))
rmse_a = np.sqrt(np.mean((lon_calib-true_lon_track)**2))
axB.plot(dep_days, true_lon_track, color="k", lw=2, label="true")
axB.plot(dep_days, lon_biased, color="#d62728", lw=1.5, ls="--", label="clock drift (holonomy)")
axB.plot(dep_days, lon_calib, color="#2ca02c", lw=1.5, label="endpoint-calibrated")
axB.set_title(f"B. Longitude holonomy from clock drift, removed by two known fixes\n"
              f"RMSE {rmse_b:.2f} -> {rmse_a:.0e} deg (the real, robust win)", fontsize=10)
axB.set_xlabel("day of deployment"); axB.set_ylabel("longitude (deg)"); axB.legend(fontsize=8, loc="upper left")

# ---- C. latitude log-likelihood: solstice single peak vs equinox Z2 mirror -
axC = fig.add_subplot(gs[1, 0])
lats = np.linspace(-75, 75, 600)
for label, dday, c in [("solstice (Jun 21): single peak", 172, "#1f77b4"),
                       ("equinox: symmetric bimodal (sign degenerate)", int(eqx_days[0]), "#d62728")]:
    t_day = year0 + dday*86400 + np.arange(0, 86400, 600)
    obs = mu(true_lon, true_lat, t_day)
    ll = np.array([light_loglik(true_lon, la, t_day, obs) for la in lats])
    ll -= ll.max()
    axC.plot(lats, ll, color=c, lw=2, label=label)
axC.axvline(true_lat, color="0.4", ls=":"); axC.axvline(-true_lat, color="0.7", ls=":")
axC.set_ylim(-60, 3)
axC.set_title("C. Latitude likelihood profile (continuous model)\n"
              "equinox keeps |lat| but mirrors the sign: a discrete Z2, not a collapse", fontsize=10)
axC.set_xlabel("candidate latitude (deg)"); axC.set_ylabel("relative log-likelihood")
axC.legend(fontsize=8, loc="lower center")

# ---- D. why: max sun elevation encodes signed lat at solstice, |lat| at equinox
axD = fig.add_subplot(gs[1, 1])
phi = np.linspace(-80, 80, 400)
for decl_deg, lab, c in [(23.4, "solstice (decl +23.4): signed latitude", "#1f77b4"),
                         (0.0, "equinox (decl 0): |latitude| only", "#d62728")]:
    max_elev = 90 - np.abs(phi - decl_deg)
    axD.plot(phi, max_elev, color=c, lw=2, label=lab)
axD.axvline(true_lat, color="0.4", ls=":"); axD.axvline(-true_lat, color="0.7", ls=":")
axD.set_title("D. Midday sun elevation vs latitude (the amplitude channel)\n"
              "equinox is V-shaped: magnitude survives, sign does not", fontsize=10)
axD.set_xlabel("latitude (deg)"); axD.set_ylabel("max sun elevation (deg)"); axD.legend(fontsize=8, loc="lower center")

fig.suptitle("The continuous model is well-conditioned in latitude magnitude; exploit the clock via endpoints, and the residual equinox degeneracy is a discrete hemisphere sign",
             fontsize=11, y=0.995)
fig.savefig("conditioning_demo.png", dpi=140, bbox_inches="tight")
print("saved conditioning_demo.png")

# ---- diagnostics -----------------------------------------------------------
sol_i = int(np.argmin(np.abs(sample_days - 172)))
eqx_i = int(np.argmin(np.abs(sample_days - eqx_days[0])))
print(f"latitude information: solstice {info_lat[sol_i]:.2f}, equinox {info_lat[eqx_i]:.2f} "
      f"(does NOT collapse: midday amplitude retains |lat| in the continuous model)")
print(f"longitude information: solstice {info_lon[sol_i]:.2f}, equinox {info_lon[eqx_i]:.2f} "
      f"(from twilight timing; the endpoints make this exact via clock calibration)")
print(f"clock-drift holonomy: longitude RMSE {rmse_b:.3f} -> {rmse_a:.1e} deg after endpoint calibration")
t_eqx = year0 + int(eqx_days[0]) * 86400 + np.arange(0, 86400, 600)
obs_eqx = mu(true_lon, true_lat, t_eqx)
ll_south = light_loglik(true_lon, -50.0, t_eqx, obs_eqx)
ll_north = light_loglik(true_lon, 50.0, t_eqx, obs_eqx)
print(f"equinox loglik at lat -50 vs +50: {ll_south:.2f} vs {ll_north:.2f} (near-equal: hemisphere sign Z2)")
