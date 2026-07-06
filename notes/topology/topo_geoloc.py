"""
Topological recasting of light-level geolocation: a prototype.

Demonstrates four claims on a lon/lat grid using the same low-precision solar
almanac as the invTwilightFree Rust core (J2000 epoch, unix 946728000):

  A. Each observation at KNOWN time is a circle of position on the sphere
     (locus of equal solar zenith). Several circles intersect at the fix.
  B. The diurnal-phase residual phi(x) = wrap(H_true - H_x) is a circle-valued
     field whose zero set is exactly the true meridian. It is LATITUDE-FREE and
     TIME-FREE: longitude is a pure phase. Its winding = longitude in turns.
  C. Latitude rides on day length, whose sensitivity d(daylength)/d(lat) -> 0 at
     the equinox (declination 0). The fiber pinches; phi cannot see it.
  D. The hemisphere monodromy: a southern track and its northern mirror give
     near-identical day length, coinciding exactly at the equinoxes (branch
     points) and separating at the solstices.
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

J2000 = 946728000.0  # unix seconds at 2000-01-01 12:00 UTC

def solar(t):
    d = (t - J2000) / 86400.0
    g = np.radians(357.529 + 0.98560028 * d)
    q = 280.459 + 0.98564736 * d
    Lsun = np.radians(q + 1.915 * np.sin(g) + 0.020 * np.sin(2 * g))
    e = np.radians(23.439 - 3.6e-7 * d)
    decl = np.arcsin(np.sin(e) * np.sin(Lsun))
    ra = np.degrees(np.arctan2(np.cos(e) * np.sin(Lsun), np.cos(Lsun)))
    gmst = (18.697374558 + 24.06570982441908 * d) * 15.0  # degrees
    return decl, ra, gmst

def zenith(lon, lat, t):
    decl, ra, gmst = solar(t)
    H = np.radians(gmst + lon - ra)
    cz = np.sin(np.radians(lat)) * np.sin(decl) + np.cos(np.radians(lat)) * np.cos(decl) * np.cos(H)
    return np.degrees(np.arccos(np.clip(cz, -1, 1)))

def hour_angle_deg(lon, t):
    """Diurnal phase (0 at local solar noon), degrees, in (-180,180]."""
    _, ra, gmst = solar(t)
    H = (gmst + lon - ra + 180.0) % 360.0 - 180.0
    return H

def day_length_hours(lat_deg, decl_rad):
    cosH0 = -np.tan(np.radians(lat_deg)) * np.tan(decl_rad)
    H0 = np.arccos(np.clip(cosH0, -1, 1))      # radians
    return 24.0 * H0 / np.pi

# ---- grid ------------------------------------------------------------------
lons = np.arange(-180, 180.01, 1.0)
lats = np.arange(-80, 80.01, 1.0)
LON, LAT = np.meshgrid(lons, lats)

# true location and a strong-signal date (austral winter solstice)
true_lon, true_lat = 150.0, -50.0
t_solstice = (np.datetime64("2024-06-21T00:00:00") - np.datetime64("1970-01-01T00:00:00")) / np.timedelta64(1, "s")
t_solstice = float(t_solstice)

fig = plt.figure(figsize=(13.5, 10))
gs = GridSpec(2, 2, figure=fig, hspace=0.28, wspace=0.22)

# ---- Panel A: circles of position + fix ------------------------------------
axA = fig.add_subplot(gs[0, 0])
obs_times = [t_solstice + h * 3600 for h in (2, 9, 16)]
colors = ["#1f77b4", "#2ca02c", "#d62728"]
for t, c in zip(obs_times, colors):
    z_field = zenith(LON, LAT, t)
    z_obs = float(zenith(true_lon, true_lat, t))
    axA.contour(LON, LAT, z_field, levels=[z_obs], colors=[c], linewidths=2)
    # subsolar point (center of the circle)
    decl, ra, gmst = solar(t)
    sub_lon = ((ra - gmst + 180) % 360) - 180
    axA.plot(sub_lon, np.degrees(decl), "*", color=c, ms=9, alpha=0.7)
axA.plot(true_lon, true_lat, "ko", ms=9, mfc="yellow", mec="k", zorder=5)
axA.annotate("fix", (true_lon, true_lat), textcoords="offset points", xytext=(8, 6), fontsize=10)
axA.set_title("A. Circles of position (known time) intersect at the fix\n"
              "each colour = one observation; ★ = subsolar centre", fontsize=10)
axA.set_xlabel("longitude"); axA.set_ylabel("latitude")
axA.set_xlim(-180, 180); axA.set_ylim(-80, 80); axA.grid(alpha=0.2)

# ---- Panel B: S1 phase residual, zero-locus = true meridian ----------------
axB = fig.add_subplot(gs[0, 1])
H_true = hour_angle_deg(true_lon, t_solstice)
H_grid = hour_angle_deg(LON, t_solstice)
phi = (H_true - H_grid + 180) % 360 - 180          # circle-valued residual, degrees
im = axB.pcolormesh(LON, LAT, phi, cmap="twilight", shading="auto", vmin=-180, vmax=180)
axB.contour(LON, LAT, phi, levels=[0], colors="k", linewidths=2)
axB.axvline(true_lon, color="k", ls=":", lw=0.8)
axB.set_title("B. Diurnal-phase residual phi(x) = wrap(H_true - H_x)\n"
              "zero-locus (black) = true meridian; latitude- and time-free", fontsize=10)
axB.set_xlabel("longitude"); axB.set_ylabel("latitude")
cb = fig.colorbar(im, ax=axB, fraction=0.046, pad=0.04); cb.set_label("phi (deg)")

# ---- Panel C: latitude rides on day length; pinch at equinox ---------------
axC = fig.add_subplot(gs[1, 0])
lat_axis = np.linspace(-75, 75, 400)
for decl_deg, lab, c in [(23.4, "solstice (decl 23.4)", "#d62728"),
                         (12.0, "decl 12", "#ff7f0e"),
                         (5.0, "decl 5", "#9467bd"),
                         (0.0, "equinox (decl 0)", "#1f77b4")]:
    D = day_length_hours(lat_axis, np.radians(decl_deg))
    axC.plot(lat_axis, D, color=c, label=lab, lw=2)
axC.axvline(true_lat, color="k", ls=":", lw=0.8)
axC.axvline(-true_lat, color="grey", ls=":", lw=0.8)
axC.annotate("true", (true_lat, 3), fontsize=9)
axC.annotate("mirror", (-true_lat, 3), fontsize=9, color="grey")
axC.set_title("C. Day length vs latitude: the fiber flattens to 12 h at the\n"
              "equinox, so latitude (and its hemisphere mirror) become unidentified", fontsize=10)
axC.set_xlabel("latitude"); axC.set_ylabel("day length (h)")
axC.legend(fontsize=8, loc="upper left"); axC.grid(alpha=0.2)

# ---- Panel D: hemisphere monodromy across the year -------------------------
axD = fig.add_subplot(gs[1, 1])
days = np.arange(0, 365)
t_year = (np.datetime64("2024-01-01") - np.datetime64("1970-01-01")) / np.timedelta64(1, "s")
t_year = float(t_year) + days * 86400 + 43200
decl_year = np.array([solar(t)[0] for t in t_year])
D_true = day_length_hours(true_lat, decl_year)
D_mirror = day_length_hours(-true_lat, decl_year)
axD.plot(days, D_true, color="#1f77b4", lw=2, label=f"true track (lat {true_lat:g})")
axD.plot(days, D_mirror, color="#d62728", lw=2, ls="--", label=f"mirror (lat {-true_lat:g})")
# mark equinoxes (decl crossing zero)
sign = np.sign(decl_year)
cross = np.where(np.diff(sign) != 0)[0]
for ci in cross:
    axD.axvline(ci, color="k", ls=":", lw=1)
    axD.annotate("equinox\n(branch point)", (ci, 22), fontsize=8, ha="center")
axD.set_title("D. Hemisphere monodromy: true and mirror day-length signals\n"
              "coincide at the equinoxes (branch points), separate at solstices", fontsize=10)
axD.set_xlabel("day of 2024"); axD.set_ylabel("day length (h)")
axD.legend(fontsize=8, loc="lower center"); axD.grid(alpha=0.2); axD.set_ylim(0, 24)

fig.suptitle("Light-level geolocation as a circle-bundle over the time-torus: phase residual, fiber pinch, and monodromy",
             fontsize=12, y=0.99)
fig.savefig("phase_residual_prototype.png", dpi=140, bbox_inches="tight")
print("saved phase_residual_prototype.png")

# ---- quantitative sanity checks --------------------------------------------
# B: residual zero-locus should be the true meridian (latitude-independent)
col = np.argmin(np.abs(lons - true_lon))
print("B. |phi| on true meridian, max over latitudes (deg):",
      round(float(np.max(np.abs(phi[:, col]))), 4), "(expect ~0)")
print("B. phi std across latitude at lon=0 (deg):",
      round(float(np.std(phi[:, np.argmin(np.abs(lons-0))])), 6), "(expect 0: time/lat-free)")
# C: day-length sensitivity to latitude collapses at equinox
for decl_deg in (23.4, 5.0, 0.0):
    D = day_length_hours(lat_axis, np.radians(decl_deg))
    sens = np.max(np.abs(np.gradient(D, lat_axis)))
    print(f"C. max |d(daylength)/d(lat)| at decl {decl_deg:>5}: {sens:.4f} h/deg")
# D: true vs mirror gap at equinox vs solstice
print("D. |D_true - D_mirror| at first equinox (h):", round(float(abs(D_true[cross[0]]-D_mirror[cross[0]])), 4))
print("D. |D_true - D_mirror| at mid-year solstice (h):", round(float(abs(D_true[172]-D_mirror[172])), 4))
