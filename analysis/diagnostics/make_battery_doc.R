# Assemble the battery review document as a single self-contained local HTML file.
# Images are embedded as data URIs so the file can be moved or mailed as one piece.
# NOT published as a hosted artifact: these are tracks and light records from the
# restricted deliveries.
suppressMessages(library(data.table))
# run from the package root (setwd removed for portability)
DIR <- "scratch/nes_calibration/battery_fig"
OUT <- "scratch/nes_calibration/light_battery_review.html"
S <- fread("scratch/nes_calibration/battery_summary.csv")

b64 <- function(f) {
  raw <- readBin(f, "raw", file.info(f)$size)
  sprintf("data:image/png;base64,%s", jsonlite::base64_enc(raw))
}
img <- function(f, cap) sprintf(
  '<figure><img src="%s" alt="%s"><figcaption>%s</figcaption></figure>',
  b64(file.path(DIR, f)), cap, cap)

ids <- unique(S$id)
tbl <- paste0(
  "<table><thead><tr>",
  paste0("<th>", c("tag", "offset (d)", "n", "start", "days", "lat range", "lon range",
                   "median atten.", "frac blacked out", "median clear-sky",
                   "median shaded"), "</th>", collapse = ""),
  "</tr></thead><tbody>",
  paste0(apply(S, 1, function(r) sprintf(
    "<tr><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s&ndash;%s</td><td>%s&ndash;%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td></tr>",
    r[["id"]], r[["offset"]], r[["n"]], r[["start"]], r[["days"]],
    r[["lat_min"]], r[["lat_max"]], r[["lon_min"]], r[["lon_max"]],
    r[["atten_med"]], r[["atten_zero"]], r[["perfect_med"]], r[["degraded_med"]])),
    collapse = ""),
  "</tbody></table>")

html <- sprintf('<!doctype html><html><head><meta charset="utf-8">
<title>Light battery for review</title>
<style>
 body{font:15px/1.6 -apple-system,Segoe UI,Roboto,sans-serif;max-width:1180px;
      margin:2rem auto;padding:0 1.2rem;color:#1c2024;background:#fff}
 h1{font-size:1.7rem;margin-bottom:.2rem} h2{font-size:1.15rem;margin-top:2.4rem;
   border-bottom:1px solid #e3e6ea;padding-bottom:.3rem}
 .sub{color:#5b6470;margin-top:0}
 figure{margin:1.2rem 0} img{width:100%%;border:1px solid #e3e6ea;border-radius:6px}
 figcaption{color:#5b6470;font-size:.87rem;margin-top:.4rem}
 table{border-collapse:collapse;width:100%%;font-size:.83rem;margin:1rem 0}
 th,td{border:1px solid #e3e6ea;padding:.32rem .5rem;text-align:right}
 th{background:#f6f8fa;text-align:center} td:first-child,td:nth-child(4){text-align:left}
 .box{background:#f6f8fa;border-left:3px solid #1f6feb;padding:.8rem 1rem;margin:1.2rem 0;
      border-radius:0 5px 5px 0}
 .warn{border-left-color:#bf8700;background:#fff8e6}
 .ok{border-left-color:#1a7f37;background:#f2fbf5}
 code{background:#f0f3f6;padding:.1rem .3rem;border-radius:3px;font-size:.88em}
 ul{padding-left:1.2rem} li{margin:.3rem 0}
</style></head><body>
<h1>Light battery for review</h1>
<p class="sub">%d records &mdash; %d real Argos tracks &times; 4 date offsets. Built %s.</p>

<div class="box ok"><b>Why filtering mattered here and not elsewhere.</b> Rescoring
the campaign&rsquo;s own error metrics against filtered truth changed nothing at all
(median 253&nbsp;km either way, mean latitude bias &minus;0.137 either way, paired
p&nbsp;=&nbsp;0.92 and 0.70) &mdash; because those are scored at 12-hourly knots, where
interpolation steps over brief spikes. This battery samples truth every
<b>30&nbsp;minutes</b>, where it does not: unfiltered, the interpolated position moves
by 8&ndash;26&nbsp;km at q95 and up to 3634&nbsp;km. Clear-sky light computed at a
position the animal never occupied is not clear-sky light.</div>

<div class="box ok"><b>Acceptance control passed.</b> At offset 0 the constructed
&ldquo;degraded&rdquo; record reproduces the tag&rsquo;s real light record exactly
(median |difference| = 0.00 units on all six tags; correlation 0.981&ndash;0.998).
The degradation model is therefore a description of these data, not an invention.
An earlier version failed this control: it classified night from the original
geometry, so a date shift transplanted real night readings into what had become
daylight.</div>

<h2>What is real and what is synthetic</h2>
<ul>
<li><b>Movement is real, and speed-filtered.</b> True positions are Argos fixes
interpolated to the observation times, after deleting fixes implying more than 3 m/s.
That threshold is not free: the fvilches delivery arrived already filtered at exactly
3 m/s, while the 2021 delivery was raw &mdash; 15&ndash;25%% of its steps implied more
than 5 m/s and one reached 159,722 m/s. Filtering removes 23&ndash;39%% of the 2021
fixes and <b>nothing</b> from fvilches, which is the control: the threshold deletes
the impossible, not the merely fast.
<br><br>A speed filter alone is not sufficient, because it permits <i>vmax times the
gap</i>. On a sparse track that is a lot of rope &mdash; at 3&nbsp;m/s a 13-hour gap
allows 140&nbsp;km, so a fix displaced 65&nbsp;km passes the speed test and leaves a
triangular out-and-back detour. Every deployment carried such spikes at the
0.2&ndash;1.4%% level, worst on the sparsest tracks (2023041: 406 fixes, 13.7&nbsp;h
median gap, worst excursion 65&nbsp;km). A second, scale-free detour filter removes
them: <b>70 fixes across all 29 deployments</b>, and the spike rate falls to zero
everywhere.
<br><br>Neither filter could reject 2023041&rsquo;s large triangular excursion, and the
reason is instructive: that fix sits inside a 213-hour gap on one side and a 352-hour
gap on the other, implying 2.36 and 2.01&nbsp;m/s &mdash; both legal, because at
3&nbsp;m/s those gaps permit 3804&nbsp;km of travel. It is simply
<b>unverifiable</b>. A third filter drops fixes with large gaps on both sides (40
across all deployments), and, more importantly, the battery now refuses to claim a
true position more than 12 hours from a fix: removing an uncorroborated fix does not
make the truth known across the gap it sat in, it only stops pretending. No random walk, no invented step-length distribution &mdash; which
is what made the earlier simulation useless (it fitted two to three times worse than
anything real, so it could not have detected the effect being tested).</li>
<li><b>Shading is real.</b> Attenuation is measured per observation as
<code>observed sky / clear-sky expectation</code> against that tag&rsquo;s own fitted
response, then re-applied. Dive bouts, cloud and haul-outs come along with it.</li>
<li><b>Night contamination is real</b> and additive: dark current above the fitted
floor, plus whatever artificial light or moonlight is genuinely in the record.</li>
<li><b>Only the date is synthetic</b>, shifted in <i>whole days</i> so local time of
day is preserved exactly &mdash; a seal that dives at dawn still dives at dawn. A
fractional shift would slide behaviour against the solar cycle and manufacture an
artefact.</li>
</ul>

<div class="box warn"><b>Limits worth your judgement.</b>
<ul>
<li>A date shift puts a real track in a season it never experienced. The movement is
real; the <i>combination</i> of movement and season is not. That is deliberate &mdash;
it breaks the phenology confound that makes the 29 real deployments unable to separate
declination from latitude and behaviour &mdash; but these are not real animals.</li>
<li>Attenuation is transplanted onto a different solar geometry. Diel alignment is
preserved; solar-elevation alignment cannot be, for any date shift.</li>
<li>Argos carries tens of km of error and is irregularly sampled, so &ldquo;truth&rdquo;
is an interpolation of a noisy reference.</li>
<li>All six tracks are northern-hemisphere post-moult migrations. The battery varies
season, not hemisphere or species.</li>
<li>Tracks are now selected on <b>truth support</b> &mdash; the fraction of samples
within 12&nbsp;hours of an Argos fix &mdash; with latitude range breaking ties. The
first version scored latitude range times fix count, which rewards how many fixes a
tag has and ignores when they are; it selected a track with only 58%% support and a
1800&nbsp;km interpolated excursion. All six now sit at 97&ndash;99%%.</li>
</ul></div>

<h2>The tracks</h2>
%s

<h2>Measured attenuation</h2>
%s
<p class="sub">Median attenuation runs 0.78&ndash;0.87, so a typical sample loses
13&ndash;22%% of the sky. The mass at zero is real: 8&ndash;26%% of lit samples are
complete blackouts, which is the animal being deep. Values slightly above 1 are sensor
noise. The cap at 3 binds on less than 0.11%% of samples, so it is not doing the
work.</p>

<h2>Clear sky vs shaded, per track and season</h2>
<p class="sub">Blue is the clear-sky expectation at the true position; red is that
same curve after the measured attenuation. Dotted verticals mark equinoxes. Each panel
is the same animal on the same path, sampled at a different time of year.</p>
%s

<h2>Diel detail</h2>
<p class="sub">Six days near an equinox and near a solstice. This is where to judge
whether the shading looks like a diving seal rather than like noise.</p>
%s

<h2>Record summary</h2>
%s

<h2>What this is for</h2>
<p>Every record has an exactly known true position, a known clear-sky light curve and
a known attenuation. That makes it possible to ask where the latitude bias comes from
against a known answer, rather than against Argos. The immediate questions it can
settle: whether the bias survives when the light is <i>perfect</i> (if so, it is the
estimator, not the light); how much of it the measured shading reproduces; and whether
the optimal shading rate really tracks declination once season is varied while the
animal and its behaviour are held fixed.</p>
</body></html>',
  nrow(S), length(ids), format(Sys.Date(), "%%d %%B %%Y"),
  img("tracks.png", "All six Argos tracks, subsampled 1:20. Filled circles mark deployment."),
  img("attenuation.png", "Measured attenuation factor: pooled distribution and per deployment."),
  paste0(vapply(ids, function(i) img(sprintf("light_%s.png", i),
    sprintf("%s: clear-sky (blue) vs shaded (red) at each of four date offsets.", i)), ""),
    collapse = ""),
  paste0(vapply(ids[1:2], function(i) img(sprintf("zoom_%s.png", i),
    sprintf("%s: six-day windows near an equinox and near a solstice.", i)), ""),
    collapse = ""),
  tbl)

writeLines(html, OUT)
cat(sprintf("wrote %s (%.1f MB)\n", OUT, file.info(OUT)$size / 1e6))
