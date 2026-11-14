/* Plotly figures for the three result panels: a port of webapp/figures.py.

   Presentation only: values from the model are plotted unchanged. Default axis
   windows match the original app; a window is widened only when data would
   otherwise fall outside it. Keep this file in step with figures.py.
*/
(function () {
  // RP blue (as in the original app), J red, non-precessing reference (zero line) black dashed.
  const RP_COLOR = "#2457d6";
  const J_COLOR = "#d93a2b";
  const REF_COLOR = "#1d2433";
  const TEXT = "#1d2433";
  const MUTED = "#5d6677";
  const GRID = "#eceff4";
  const AXIS = "#c9cfd9";
  const FONT = "Inter, system-ui, -apple-system, 'Segoe UI', Roboto, sans-serif";

  const FREQ_TICKS = [20, 30, 40, 60, 80, 100, 150, 200, 300, 400];

  // B_RP/B_NP − 1 diverges where the non-precessing amplitude B_NP vanishes. Beyond this
  // value the axis stops growing and the panel is labelled instead.
  const AMP_DISPLAY_LIMIT = 20.0;
  const POINT_TOL = 1e-9;
  const STRAIN_VIEW = 0.5;  // s before t_c shown when the strain panel opens

  // NumPy semantics: min/max/ptp of an array containing NaN are NaN.
  function ptp(a) {
    let lo = Infinity, hi = -Infinity;
    for (const v of a) { if (isNaN(v)) return NaN; if (v < lo) lo = v; if (v > hi) hi = v; }
    return hi - lo;
  }
  function maxAbs(a) {
    let m = -Infinity;
    for (const v of a) { if (isNaN(v)) return NaN; m = Math.max(m, Math.abs(v)); }
    return m;
  }
  function nanMaxAbs(a) {  // np.nanmax(np.abs(a)): NaN only if every value is NaN
    let m = NaN;
    for (const v of a) { if (!isNaN(v)) m = isNaN(m) ? Math.abs(v) : Math.max(m, Math.abs(v)); }
    return m;
  }

  function isPoint(x, y) {
    return x.length > 0 && ptp(x) < POINT_TOL && ptp(y) < POINT_TOL;
  }

  function isFlatZero(...arrays) {
    return arrays.every((a) => a.length && maxAbs(a) < POINT_TOL);
  }

  function ampPeak(r) {
    const peak = r.f.length ? nanMaxAbs(r.amp) : 0;
    return peak > AMP_DISPLAY_LIMIT ? peak : null;
  }

  // Python's format(v, ".{p}g").
  function formatG(v, p) {
    if (v === 0) return "0";
    if (!isFinite(v)) return isNaN(v) ? "nan" : (v > 0 ? "inf" : "-inf");
    const [mant, e] = v.toExponential(p - 1).split("e");
    const exp = parseInt(e, 10);
    if (exp < -4 || exp >= p) {
      const m = mant.includes(".") ? mant.replace(/\.?0+$/, "") : mant;
      return `${m}e${exp < 0 ? "-" : "+"}${String(Math.abs(exp)).padStart(2, "0")}`;
    }
    const fixed = v.toFixed(Math.max(0, p - 1 - exp));
    return fixed.includes(".") ? fixed.replace(/\.?0+$/, "") : fixed;
  }

  // e.g. 3.9×10<sup>15</sup> for Plotly/HTML text.
  function sci(v) {
    if (v === 0 || (Math.abs(v) >= 1e-3 && Math.abs(v) < 1e4)) return formatG(v, 3);
    const exp = Math.floor(Math.log10(Math.abs(v)));
    return `${formatG(v / 10 ** exp, 2)}×10<sup>${exp}</sup>`;
  }

  function baseLayout(xTitle, yTitle) {
    const axis = (text) => ({
      showline: true, linecolor: AXIS, linewidth: 1, mirror: true,
      ticks: "outside", tickcolor: AXIS, ticklen: 4,
      gridcolor: GRID, zeroline: false,
      title: { text, font: { size: 13, color: TEXT }, standoff: 8 },
      tickfont: { size: 11.5, color: MUTED },
    });
    return {
      font: { family: FONT, size: 12, color: TEXT },
      paper_bgcolor: "rgba(0,0,0,0)",
      plot_bgcolor: "#ffffff",
      margin: { l: 62, r: 18, t: 14, b: 50 },
      showlegend: false,
      autosize: true,
      hoverlabel: { bgcolor: "#ffffff", bordercolor: AXIS, font: { family: FONT, size: 12, color: TEXT } },
      xaxis: axis(xTitle),
      yaxis: axis(yTitle),
      uirevision: "keep",
      shapes: [],
      annotations: [],
    };
  }

  // Default window, widened to include all finite data plus a small margin.
  function fit(defaultLo, defaultHi, arrays, padFrac = 0.04) {
    let lo = Infinity, hi = -Infinity, any = false;
    for (const a of arrays) {
      for (const v of a) { if (isFinite(v)) { any = true; if (v < lo) lo = v; if (v > hi) hi = v; } }
    }
    if (!any) return [defaultLo, defaultHi];
    const pad = padFrac * Math.max(hi - lo, 1e-9);
    return [Math.min(defaultLo, lo - pad), Math.max(defaultHi, hi + pad)];
  }

  function gcd(a, b) { return b ? gcd(b, a % b) : Math.abs(a); }

  function piLabel(num, den) {
    if (num === 0) return "0";
    const g = gcd(num, den);
    num /= g; den /= g;
    const sign = num < 0 ? "−" : "";
    num = Math.abs(num);
    const n = num === 1 ? "" : String(num);
    return den === 1 ? `${sign}${n}π` : `${sign}${n}π/${den}`;
  }

  function piTicks(lo, hi, maxTicks = 7) {
    let den, start, stop;
    for (den of [8, 4, 2, 1]) {
      start = Math.ceil(lo / (Math.PI / den));
      stop = Math.floor(hi / (Math.PI / den));
      if (stop - start + 1 <= maxTicks) break;
    }
    const vals = [], text = [];
    for (let i = start; i <= stop; i++) { vals.push(i / den * Math.PI); text.push(piLabel(i, den)); }
    return [vals, text];
  }

  // Insert gaps where φ_L jumps across ±π so lines don't streak across the panel.
  function breakWraps(phi, y) {
    const xs = [], ys = [];
    for (let i = 0; i < phi.length; i++) {
      if (i > 0 && Math.abs(phi[i] - phi[i - 1]) > Math.PI) { xs.push(NaN); ys.push(NaN); }
      xs.push(phi[i]); ys.push(y[i]);
    }
    return [xs, ys];
  }

  function modulationFigure(r) {
    const [phiJ, cosJ] = r.J;
    // φ_L comes from arctan2 and lies in (−π, π]; show J at the same angle in that convention.
    const phiJPlot = Math.atan2(Math.sin(phiJ), Math.cos(phiJ));
    const layout = baseLayout("φ<sub>L</sub>", "cos θ<sub>L</sub>");

    const point = isPoint(r.phi_L, r.cos_L);
    const [x, y] = point ? [[r.phi_L[0]], [r.cos_L[0]]] : breakWraps(r.phi_L, r.cos_L);
    const data = [
      { type: "scatter", x: [phiJPlot], y: [cosJ], mode: "markers", name: "J",
        marker: { color: J_COLOR, size: 14, line: { color: "#ffffff", width: 1.5 } },
        hovertemplate: "J<br>φ<sub>J</sub> = %{x:.4f}<br>cos θ<sub>J</sub> = %{y:.4f}<extra></extra>" },
      { type: "scatter", x, y, mode: point ? "markers" : "lines", name: "RP",
        line: { color: RP_COLOR, width: 2 }, marker: { color: RP_COLOR, size: 8 },
        hovertemplate: "RP<br>φ<sub>L</sub> = %{x:.4f}<br>cos θ<sub>L</sub> = %{y:.4f}<extra></extra>" },
    ];

    const xr = fit(phiJPlot - Math.PI / 4, phiJPlot + Math.PI / 4, [r.phi_L]);
    let yr = fit(cosJ - 0.8, cosJ + 0.8, [r.cos_L]);
    yr = [Math.max(yr[0], -1.05), Math.min(yr[1], 1.05)];
    const [tickvals, ticktext] = piTicks(xr[0], xr[1]);
    Object.assign(layout.xaxis, { range: xr, tickvals, ticktext });
    layout.yaxis.range = yr;
    layout.hovermode = "closest";
    return { data, layout };
  }

  function frequencyFigure(r, y, yTitle, defaultY, digits, limit) {
    const layout = baseLayout("f (Hz)", yTitle);
    // layer "below" keeps an RP curve that coincides with the reference visible on top of it.
    layout.shapes.push({ type: "line", xref: "x domain", x0: 0, x1: 1, yref: "y", y0: 0, y1: 0,
                         line: { color: REF_COLOR, width: 1, dash: "dash" }, opacity: 0.7, layer: "below" });
    const data = [{ type: "scatter", x: r.f, y, mode: "lines", name: "RP",
                    line: { color: RP_COLOR, width: 2 },
                    hovertemplate: `RP: %{y:.${digits}f}<extra></extra>` }];

    const fHi = r.f.length ? Math.max(150.0, r.f[r.f.length - 1] * 1.03) : 150.0;
    const ticks = FREQ_TICKS.filter((t) => t <= fHi);
    Object.assign(layout.xaxis, { type: "log", range: [Math.log10(19), Math.log10(fHi)],
                                  tickvals: ticks, ticktext: ticks.map(String), hoverformat: ".2f" });
    let yr = fit(defaultY[0], defaultY[1], [y]);
    if (limit !== undefined && (yr[0] < -limit || yr[1] > limit)) {
      yr = [Math.max(yr[0], -limit), Math.min(yr[1], limit)];
      layout.annotations.push({
        text: `Off scale: values reach ${sci(nanMaxAbs(y))}<br>(B<sub>NP</sub> ≈ 0 here, so the ratio diverges)`,
        xref: "paper", yref: "paper", x: 0.5, y: 0.97, xanchor: "center", yanchor: "top", showarrow: false,
        align: "center", font: { size: 11.5, color: "#6d4500" },
        bgcolor: "#fff8eb", bordercolor: "#f1d49a", borderwidth: 1, borderpad: 5 });
    }
    layout.yaxis.range = yr;
    layout.hovermode = "x unified";
    return { data, layout };
  }

  function amplitudeFigure(r) {
    return frequencyFigure(r, r.amp, "B<sub>RP</sub> / B<sub>NP</sub> − 1", [-1.05, 0.95], 4, AMP_DISPLAY_LIMIT);
  }

  function phaseFigure(r) {
    return frequencyFigure(r, r.phase, "Φ<sub>RP</sub> − Φ<sub>NP</sub> (rad)", [-3.25, 0.75], 4);
  }

  // Time-domain strain h(t) for RP and NP (webapp/time_domain.py). No Python counterpart yet:
  // the Dash app has no time-domain view.
  function strainFigure(r) {
    let peak = 0;
    for (const a of [r.h_RP, r.h_NP]) for (const v of a) if (isFinite(v)) peak = Math.max(peak, Math.abs(v));
    // Plot in units of a power of ten so the axis reads e.g. 4 rather than 4e−23.
    const exp = peak > 0 ? Math.floor(Math.log10(peak)) : 0;
    const scale = 10 ** -exp;
    const scaled = (a) => Array.from(a, (v) => v * scale);
    const unit = exp ? ` (10<sup>${exp < 0 ? "−" : ""}${Math.abs(exp)}</sup>)` : "";

    const layout = baseLayout("t − t<sub>c</sub> (s)", "strain h" + unit);
    layout.shapes.push({ type: "line", xref: "x domain", x0: 0, x1: 1, yref: "y", y0: 0, y1: 0,
                         line: { color: AXIS, width: 1 }, layer: "below" });
    const t = Array.from(r.t);
    const data = [
      { type: "scatter", x: t, y: scaled(r.h_NP), mode: "lines", name: "NP",
        line: { color: REF_COLOR, width: 1.2 }, opacity: 0.55,
        hovertemplate: "NP: %{y:.4f}<extra></extra>" },
      { type: "scatter", x: t, y: scaled(r.h_RP), mode: "lines", name: "RP",
        line: { color: RP_COLOR, width: 1.3 },
        hovertemplate: "RP: %{y:.4f}<extra></extra>" },
    ];
    const lim = 1.08 * peak * scale || 1;
    // Open on the last STRAIN_VIEW seconds, where RP and NP visibly separate; the range bar shows it all.
    Object.assign(layout.xaxis, { range: [Math.max(t[0], -STRAIN_VIEW), t[t.length - 1]], hoverformat: ".4f",
                                  rangeslider: { visible: true, thickness: 0.1, bgcolor: "#f8f9fb", bordercolor: AXIS, borderwidth: 1 } });
    layout.yaxis.range = [-lim, lim];
    // Fixed bottom margin for the range bar and axis title: with automatic expansion Plotly 3.1 grows
    // the range bar's margin a little on every relayout, which shrinks the plot during playback.
    layout.margin = { l: 62, r: 18, t: 14, b: 118, autoexpand: false };
    layout.hovermode = "x unified";
    return { data, layout };
  }

  // ---------- J-frame sphere (webapp/sphere.py) ----------
  // Unit sphere with J at the north pole: the band of polar angles L passes through, L's path,
  // the band of latitudes swept by the orbital plane (perpendicular to L), N̂, and a marker for
  // L and its orbital plane at one time. Needs the gl3d Plotly bundle. No Python counterpart.
  const SPHERE = {
    base: "#eef1f6", lBand: "#bfd1f6", orbitBand: "#f6dcae", overlap: "#d8cde6",
    orbit: "#b36200", N: REF_COLOR, equator: "#aab3c2",
  };
  const S_PATH = 1.004, S_MARK = 1.012;  // radii that keep lines and markers above the surface
  const ARROW = 1.4;                     // length of the J and N arrows (sphere radius 1)
  const EPS = 1e-6;
  const SPHERE_TRACE = [7, 8, 6];        // traces moved by sphereMoving: L marker, L arrow, orbital plane
  const TRAIL_TRACE = 11;                // path of L up to the current time (playback)
  const PATH_TRACE = 2;                  // full path of L, faded during playback

  function sphereBands(r) {
    let lo = Infinity, hi = -Infinity;
    for (const v of r.sp_theta) { if (isFinite(v)) { lo = Math.min(lo, v); hi = Math.max(hi, v); } }
    return { lo, hi, lBand: hi - lo > EPS, orbitBand: hi > EPS, overlap: hi > Math.PI / 2 - hi + EPS };
  }

  function unitVec(theta, phi) {
    return [Math.sin(theta) * Math.cos(phi), Math.sin(theta) * Math.sin(phi), Math.cos(theta)];
  }

  function scaled3(v, s) { return [v[0] * s, v[1] * s, v[2] * s]; }

  // Sphere surface coloured by band; rows are added at the band edges so the edges are sharp.
  function sphereSurface(b) {
    const rows = new Set();
    for (let i = 0; i <= 90; i++) rows.add(i * Math.PI / 90);
    const edges = [];
    if (b.lBand) edges.push(b.lo, b.hi);
    if (b.orbitBand) edges.push(Math.PI / 2 - b.hi, Math.PI / 2 + b.hi);
    for (const e of edges) for (const d of [-EPS, EPS]) rows.add(Math.min(Math.PI, Math.max(0, e + d)));
    const thetas = [...rows].sort((a, c) => a - c);
    const nPhi = 73;
    const x = [], y = [], z = [], c = [];
    for (const th of thetas) {
      const inL = b.lBand && th >= b.lo && th <= b.hi;
      const inO = b.orbitBand && Math.abs(Math.PI / 2 - th) <= b.hi;
      const cat = (inL ? 1 : 0) + (inO ? 2 : 0);
      const xr = [], yr = [], zr = [], cr = [];
      for (let j = 0; j < nPhi; j++) {
        const ph = j * 2 * Math.PI / (nPhi - 1);
        xr.push(Math.sin(th) * Math.cos(ph)); yr.push(Math.sin(th) * Math.sin(ph)); zr.push(Math.cos(th)); cr.push(cat);
      }
      x.push(xr); y.push(yr); z.push(zr); c.push(cr);
    }
    return {
      type: "surface", x, y, z, surfacecolor: c, cmin: 0, cmax: 3, showscale: false, hoverinfo: "skip",
      colorscale: [[0, SPHERE.base], [1 / 3, SPHERE.lBand], [2 / 3, SPHERE.orbitBand], [1, SPHERE.overlap]],
      lighting: { ambient: 0.8, diffuse: 0.35, specular: 0.05, roughness: 1, fresnel: 0 },
      contours: { x: { highlight: false }, y: { highlight: false }, z: { highlight: false } },
    };
  }

  function line3d(points, color, width, name) {
    return { type: "scatter3d", mode: "lines", name, hoverinfo: "skip", showlegend: false,
             x: points.map((p) => p[0]), y: points.map((p) => p[1]), z: points.map((p) => p[2]),
             line: { color, width } };
  }

  // Great circle perpendicular to the unit vector L (the orbital plane), at radius S_PATH.
  function orbitCircle(L) {
    let e1 = [L[1], -L[0], 0];  // L × z
    const n = Math.hypot(e1[0], e1[1]);
    e1 = n > 1e-9 ? [e1[0] / n, e1[1] / n, 0] : [1, 0, 0];
    const e2 = [L[1] * e1[2] - L[2] * e1[1], L[2] * e1[0] - L[0] * e1[2], L[0] * e1[1] - L[1] * e1[0]];
    const pts = [];
    for (let k = 0; k <= 120; k++) {
      const s = k * 2 * Math.PI / 120;
      pts.push([0, 1, 2].map((i) => S_PATH * (Math.cos(s) * e1[i] + Math.sin(s) * e2[i])));
    }
    return pts;
  }

  // State of L at time t: t(f) is increasing, so invert it by binary search and interpolate
  // f, θ_LJ and φ_LJ linearly between grid frequencies. null outside [t(20 Hz), t(f_cut)].
  function sphereAt(r, t) {
    const ts = r.sp_t, n = ts ? ts.length : 0;
    if (!n || !r.sp_t_monotonic || !(t >= ts[0] && t <= ts[n - 1])) return null;
    let lo = 0, hi = n - 1;
    while (hi - lo > 1) { const mid = (lo + hi) >> 1; if (ts[mid] <= t) lo = mid; else hi = mid; }
    const w = hi > lo ? (t - ts[lo]) / (ts[hi] - ts[lo]) : 0;
    const lerp = (a) => a[lo] + w * (a[hi] - a[lo]);
    return { t, f: lerp(r.f), theta: lerp(r.sp_theta), phi: lerp(r.sp_phi) };
  }

  // Coordinates for the traces in SPHERE_TRACE at state s (empty when s is null).
  function sphereMoving(s) {
    if (!s) return { x: [[], [], []], y: [[], [], []], z: [[], [], []] };
    const L = unitVec(s.theta, s.phi);
    const m = scaled3(L, S_MARK), a1 = scaled3(L, 1.3);
    const c = orbitCircle(L);
    return {
      x: [[m[0]], [L[0], a1[0]], c.map((p) => p[0])],
      y: [[m[1]], [L[1], a1[1]], c.map((p) => p[1])],
      z: [[m[2]], [L[2], a1[2]], c.map((p) => p[2])],
    };
  }

  function sphereFigure(r, s) {
    const b = sphereBands(r);
    const N = r.sp_N;
    // L's path; long signals are thinned to at most ~8000 points.
    const stride = Math.max(1, Math.ceil(r.sp_theta.length / 8000));
    const path = [], pathT = [];
    const add = (i) => { path.push(scaled3(unitVec(r.sp_theta[i], r.sp_phi[i]), S_PATH)); pathT.push(r.sp_t[i]); };
    for (let i = 0; i < r.sp_theta.length; i += stride) add(i);
    const last = r.sp_theta.length - 1;
    if (last % stride) add(last);

    const equator = [];
    for (let k = 0; k <= 120; k++) { const a = k * 2 * Math.PI / 120; equator.push([Math.cos(a), Math.sin(a), 0]); }

    const mv = sphereMoving(s);
    const tipN = scaled3(N, ARROW);
    const data = [
      sphereSurface(b),                                                         // 0
      line3d(equator, SPHERE.equator, 1.5, "Equator"),                          // 1
      line3d(path, RP_COLOR, 4, "Path of L"),                                   // 2
      line3d([[0, 0, 1], [0, 0, ARROW]], J_COLOR, 7, "J"),                      // 3
      line3d([N, scaled3(N, ARROW - 0.12)], SPHERE.N, 6, "N"),                  // 4
      { type: "cone", x: [tipN[0]], y: [tipN[1]], z: [tipN[2]], u: [N[0]], v: [N[1]], w: [N[2]],  // 5
        anchor: "tip", sizemode: "absolute", sizeref: 0.16, showscale: false, hoverinfo: "skip",
        colorscale: [[0, SPHERE.N], [1, SPHERE.N]], lighting: { ambient: 0.9, diffuse: 0.2 } },
      { type: "scatter3d", mode: "lines", hoverinfo: "skip", showlegend: false, name: "Orbital plane now",  // 6
        x: mv.x[2], y: mv.y[2], z: mv.z[2], line: { color: SPHERE.orbit, width: 5 } },
      { type: "scatter3d", mode: "markers+text", hoverinfo: "skip", showlegend: false, name: "L now",       // 7
        x: mv.x[0], y: mv.y[0], z: mv.z[0], text: ["L"], textposition: "top center",
        textfont: { color: RP_COLOR, size: 14, family: FONT },
        marker: { color: RP_COLOR, size: 6, line: { color: "#ffffff", width: 1 } } },
      { type: "scatter3d", mode: "lines", hoverinfo: "skip", showlegend: false, name: "L arrow",           // 8
        x: mv.x[1], y: mv.y[1], z: mv.z[1], line: { color: RP_COLOR, width: 7 } },
      { type: "scatter3d", mode: "markers+text", hoverinfo: "skip", showlegend: false, name: "Labels",     // 9
        x: [0, N[0] * S_MARK, tipN[0] * 1.1], y: [0, N[1] * S_MARK, tipN[1] * 1.1], z: [1, N[2] * S_MARK, tipN[2] * 1.1],
        text: ["", "", "N"], textposition: "middle center", textfont: { color: SPHERE.N, size: 14, family: FONT },
        marker: { color: [J_COLOR, SPHERE.N, "rgba(0,0,0,0)"], size: [6, 6, 1] } },
      { type: "scatter3d", mode: "text", hoverinfo: "skip", showlegend: false, name: "J label",            // 10
        x: [0], y: [0], z: [ARROW + 0.12], text: ["J"], textfont: { color: J_COLOR, size: 14, family: FONT } },
      { type: "scatter3d", mode: "lines", hoverinfo: "skip", showlegend: false, name: "Path of L so far",   // 11
        x: [], y: [], z: [], line: { color: RP_COLOR, width: 5 } },
    ];
    const hidden = { visible: false, showspikes: false, range: [-1.55, 1.55] };
    // Camera 35° in azimuth from N (which lies in the y–z plane) and 28° above the equator, J up.
    const az = Math.atan2(N[1], N[0]) - 35 * Math.PI / 180, el = 28 * Math.PI / 180, dist = 1.3;
    const layout = {
      font: { family: FONT, size: 12, color: TEXT },
      paper_bgcolor: "rgba(0,0,0,0)",
      margin: { l: 0, r: 0, t: 0, b: 0 },
      showlegend: false,
      autosize: true,
      uirevision: "keep",
      scene: {
        xaxis: hidden, yaxis: hidden, zaxis: hidden, aspectmode: "cube", hovermode: false, dragmode: "orbit",
        camera: { eye: { x: dist * Math.cos(el) * Math.cos(az), y: dist * Math.cos(el) * Math.sin(az), z: dist * Math.sin(el) },
                  up: { x: 0, y: 0, z: 1 }, center: { x: 0, y: 0, z: 0 } },
      },
    };
    return { data, layout, bands: b, path: { pts: path, t: pathT } };
  }

  // Path of L from the start of the signal up to state s (for TRAIL_TRACE), ending at L itself.
  function sphereTrail(path, s) {
    if (!s) return { x: [], y: [], z: [] };
    let lo = 0, hi = path.t.length;  // number of path points with t <= s.t
    while (lo < hi) { const mid = (lo + hi) >> 1; if (path.t[mid] <= s.t) lo = mid + 1; else hi = mid; }
    const pts = path.pts.slice(0, lo);
    pts.push(scaled3(unitVec(s.theta, s.phi), S_PATH));
    return { x: pts.map((p) => p[0]), y: pts.map((p) => p[1]), z: pts.map((p) => p[2]) };
  }

  // Axis-less figure used before the first computation and in empty/error states.
  function placeholderFigure(message) {
    const layout = { paper_bgcolor: "rgba(0,0,0,0)", plot_bgcolor: "rgba(0,0,0,0)",
                     margin: { l: 0, r: 0, t: 0, b: 0 }, autosize: true,
                     xaxis: { visible: false }, yaxis: { visible: false }, annotations: [] };
    if (message) {
      layout.annotations.push({ text: message, xref: "paper", yref: "paper", x: 0.5, y: 0.5, showarrow: false,
                                font: { family: FONT, size: 13, color: MUTED } });
    }
    return { data: [], layout };
  }

  function graphConfig(filename) {
    return {
      displaylogo: false,
      responsive: true,
      modeBarButtonsToRemove: ["select2d", "lasso2d"],
      toImageButtonOptions: { format: "png", filename, scale: 3 },
    };
  }

  window.RPFigures = {
    AMP_DISPLAY_LIMIT, isPoint, isFlatZero, ampPeak, formatG,
    modulationFigure, amplitudeFigure, phaseFigure, strainFigure, placeholderFigure, graphConfig,
    sphereFigure, sphereAt, sphereMoving, sphereTrail, SPHERE, SPHERE_TRACE, TRAIL_TRACE, PATH_TRACE,
  };
})();
