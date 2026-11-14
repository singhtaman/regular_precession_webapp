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
    modulationFigure, amplitudeFigure, phaseFigure, placeholderFigure, graphConfig,
  };
})();
