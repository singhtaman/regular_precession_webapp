/* Regular Precession explorer, browser build: controls, results and export.

   The page mirrors webapp/layout.py and webapp/callbacks.py; the calculation runs
   in worker.js (Pyodide) with the same Python model as the Dash app.
*/
(function () {
  const F = window.RPFigures;
  const $ = (sel, root = document) => root.querySelector(sel);

  const worker = new Worker("worker.js");
  const graphs = {
    mod: { el: $("#mod-graph"), file: "rp_modulation" },
    amp: { el: $("#amp-graph"), file: "rp_amplitude" },
    phase: { el: $("#phase-graph"), file: "rp_phase" },
    strain: { el: $("#strain-graph"), file: "rp_strain" },
  };

  let P = null;              // params.py, as JSON
  let grid = null;           // {fMin, deltaF} from physics.py
  const state = {};          // slider id -> current value
  const controls = {};       // slider id -> {range, number, readout, param}
  let domain = "frequency";  // "frequency" | "time": which result panels are shown
  let busy = false, pending = null, runId = 0;

  // The J-frame sphere needs Plotly's 3D traces, which plotly-basic lacks. The gl3d bundle (1.6 MB)
  // is loaded once the model is ready. It also defines window.Plotly, so the basic build is put
  // back afterwards and the two are kept apart: Plot3D draws only the sphere.
  const GL3D = "https://cdnjs.cloudflare.com/ajax/libs/plotly.js/3.1.0/plotly-gl3d.min.js";
  const sphereEl = $("#sphere-graph");
  let Plot3D = null, gl3dLoading = null;
  let sphereResult = null;   // last result with sphere data (either view)
  let sphereTime = null;     // time shown on the sphere (null: start of the signal)
  let sphereDrawn = false, sphereFrame = 0;

  // ---------- Small helpers ----------

  function h(tag, attrs = {}, ...children) {
    const el = document.createElement(tag);
    for (const [k, v] of Object.entries(attrs)) {
      if (v === null || v === undefined || v === false) continue;
      if (k === "class") el.className = v;
      else if (k === "html") el.innerHTML = v;
      else el.setAttribute(k, v === true ? "" : v);
    }
    for (const c of children.flat()) {
      if (c !== null && c !== undefined) el.append(c instanceof Node ? c : document.createTextNode(c));
    }
    return el;
  }
  const sub = (t) => h("sub", {}, t);

  function unitLabel(unit) {
    return unit === "M☉" ? ["M", sub("⊙")] : [unit];
  }

  // Number of decimals in a step, as rc-slider (Dash's slider) uses to round values.
  function decimals(step) {
    const s = String(step);
    if (s.includes("e-")) return parseInt(s.split("e-")[1], 10);
    return s.includes(".") ? s.split(".")[1].length : 0;
  }

  function valueAt(p, k) {
    return parseFloat((p.min + k * p.step).toFixed(decimals(p.step)));
  }

  function indexOf(p, v) {
    const n = Math.round((p.max - p.min) / p.step);
    return Math.min(n, Math.max(0, Math.round((v - p.min) / p.step)));
  }

  // Angle readout as a multiple of π (sliders step in π/16) plus its decimal value; from webapp/assets/app.js.
  function formatAngle(v) {
    if (v === undefined || v === null || isNaN(v)) return "—";
    const n = Math.round(v * 16 / Math.PI);
    if (Math.abs(n * Math.PI / 16 - v) > 1e-9) return v.toFixed(3) + " rad";
    if (n === 0) return "0";
    const gcd = (a, b) => (b ? gcd(b, a % b) : a);
    const g = gcd(Math.abs(n), 16), num = n / g, den = 16 / g;
    return (num === 1 ? "" : num) + "π" + (den === 1 ? "" : "/" + den) + " · " + v.toFixed(3);
  }

  // ---------- Controls (webapp/layout.py: control, group) ----------

  function control(p) {
    const n = Math.round((p.max - p.min) / p.step);
    const range = h("input", { type: "range", min: 0, max: n, step: 1, value: indexOf(p, p.value),
                               id: p.id, "aria-label": `${p.name}${p.unit ? " (" + p.unit + ")" : ""}`, disabled: true });
    const readout = p.angle ? h("span", { class: "control-value", "aria-hidden": "true" }) : null;
    const number = p.angle ? null : h("input", { type: "number", class: "control-number", step: "any",
                                                 value: p.value, "aria-label": p.name + " value", disabled: true });
    const marks = h("div", { class: "marks", "aria-hidden": "true" },
      p.marks.map(([v, label]) => h("span", { style: `--pos:${(v - p.min) / (p.max - p.min)}` }, label)));

    controls[p.id] = { range, number, readout, param: p };

    range.addEventListener("input", () => show(p.id, valueAt(p, +range.value)));
    range.addEventListener("change", () => { setValue(p.id, valueAt(p, +range.value)); schedule(); });
    if (number) {
      number.addEventListener("change", () => {
        const v = number.value === "" ? null : parseFloat(number.value);
        setValue(p.id, v);
        schedule();
      });
      number.addEventListener("keydown", (e) => { if (e.key === "Enter") number.blur(); });
    }

    return h("div", { class: "control" },
      h("div", { class: "control-head" },
        h("label", { class: "control-label", for: p.id },
          h("span", { class: "sym" }, p.symbol, p.sub ? sub(p.sub) : null),
          h("span", { class: "control-name" }, p.name),
          p.unit && !p.angle ? h("span", { class: "control-unit" }, unitLabel(p.unit)) : null),
        readout),
      h("div", { class: "slider-row" }, h("div", { class: "slider-track" }, range, marks), number),
      h("p", { class: "control-hint" }, p.hint));
  }

  function group(g) {
    let extra = null;
    if (g.key === "orientation") {
      extra = h("div", { class: "presets" },
        h("span", { class: "presets-label" }, "Presets"),
        Object.entries(P.presets).map(([name, v]) => {
          const b = h("button", { type: "button", class: "chip", disabled: true,
                                  title: `Set θ_S, φ_S, θ_J, φ_J to the ${v.label.toLowerCase()} system in systems_lib.py` },
            h("span", {}, v.label), h("small", {}, v.detail[0], sub(v.detail[1]), v.detail[2]));
          b.addEventListener("click", () => apply(v.values));
          return b;
        }));
    }
    const reset = h("button", { type: "button", class: "link-button", disabled: true }, "Reset " + g.title.toLowerCase());
    reset.addEventListener("click", () => apply(Object.fromEntries(g.params.map((p) => [p.id, p.value]))));
    return h("details", { class: "group", open: true },
      h("summary", { class: "group-summary" },
        h("span", { class: "group-title" }, g.title), h("span", { class: "group-desc" }, g.desc)),
      h("div", { class: "group-body" }, extra, g.params.map(control), reset));
  }

  function buildControls() {
    const aside = $("#controls");
    for (const g of P.groups) aside.append(group(g));
    for (const g of P.groups) for (const p of g.params) setValue(p.id, p.value);
  }

  function show(id, v) {
    const c = controls[id];
    if (c.readout) c.readout.textContent = formatAngle(v);
    if (c.number && v !== null) c.number.value = v;
  }

  function setValue(id, v) {
    const c = controls[id];
    state[id] = v;
    if (v !== null && !isNaN(v)) c.range.value = indexOf(c.param, v);
    show(id, v);
  }

  function apply(values) {
    for (const [id, v] of Object.entries(values)) setValue(id, v);
    schedule();
  }

  function enable() {
    document.querySelectorAll("#controls input, #controls button, #reset-all").forEach((el) => { el.disabled = false; });
  }

  // ---------- Running the model (one request at a time; the latest values win) ----------

  function currentValues() {
    return P.callOrder.map((id) => (state[id] === undefined || isNaN(state[id]) ? null : state[id]));
  }

  function schedule() {
    if (!grid) return;
    const request = { type: "run", values: currentValues(), domain };
    if (busy) { pending = request; return; }
    send(request);
  }

  function send(request) {
    busy = true;
    $("#plot-grid").classList.add("is-busy");
    worker.postMessage({ ...request, id: ++runId });
  }

  function done() {
    busy = false;
    if (pending) { const r = pending; pending = null; send(r); }
    else $("#plot-grid").classList.remove("is-busy");
  }

  function setDomain(d) {
    domain = d;
    document.querySelectorAll(".domain-switch button").forEach((b) => b.setAttribute("aria-pressed", String(b.dataset.domain === d)));
    document.querySelectorAll(".fd-only").forEach((el) => { el.hidden = d !== "frequency"; });
    document.querySelectorAll(".td-only").forEach((el) => { el.hidden = d !== "time"; });
    $("#plot-grid").classList.toggle("is-time", d === "time");
    pause();
    schedule();
  }

  // ---------- Results (webapp/callbacks.py: update) ----------

  function draw(key, fig) {
    const g = graphs[key];
    Plotly.react(g.el, fig.data, fig.layout, F.graphConfig(g.file));
  }

  function blank(message) {
    for (const key of Object.keys(graphs)) draw(key, F.placeholderFigure(message));
    sphereResult = null;
    clearSphere("");  // the notice above the plots explains the empty state
  }

  // ---------- J-frame sphere (both views) ----------

  function loadGl3d() {
    if (!gl3dLoading) {
      gl3dLoading = new Promise((resolve, reject) => {
        const basic = window.Plotly;
        const s = document.createElement("script");
        s.src = GL3D;
        s.onload = () => { Plot3D = window.Plotly; window.Plotly = basic; resolve(); };
        s.onerror = () => { window.Plotly = basic; reject(new Error("Could not download the 3D plotting library")); };
        document.head.append(s);
      }).then(() => { if (sphereResult) drawSphere(sphereResult); },
              (err) => { gl3dLoading = null; clearSphere(err.message + "; reload the page to try again."); });
    }
    return gl3dLoading;
  }

  function clearSphere(message) {
    if (Plot3D && sphereDrawn) Plot3D.purge(sphereEl);
    sphereDrawn = false;
    $("#sphere-readout").textContent = message || "";
  }

  function fmt(v, d) { return v.toFixed(d).replace("-", "−"); }

  function sphereReadout(r, s, requested) {
    const el = $("#sphere-readout");
    if (!r.sp_t_monotonic) {
      el.textContent = "The time map t(f) is not monotonic for these parameters, so the sphere is not linked to the other plots.";
    } else if (!s) {
      el.replaceChildren(`t − t`, h("sub", {}, "c"), ` = ${fmt(requested, 3)} s is outside the band `
        + `(${fmt(r.sp_t[0], 3)} to ${fmt(r.sp_t[r.sp_t.length - 1], 4)} s), so no L is marked.`);
    } else {
      const twoPi = 2 * Math.PI, phiMod = ((s.phi % twoPi) + twoPi) % twoPi;
      el.replaceChildren(`t − t`, h("sub", {}, "c"), ` = ${fmt(s.t, 3)} s · f = ${s.f.toFixed(2)} Hz · θ`, h("sub", {}, "LJ"),
        ` = ${s.theta.toFixed(3)} rad · φ`, h("sub", {}, "LJ"), ` = ${phiMod.toFixed(3)} rad (mod 2π)`);
    }
  }

  function stateAt(r) {
    const t = sphereTime === null ? r.sp_t[0] : sphereTime;
    return [F.sphereAt(r, t), t];
  }

  function drawSphere(r) {
    if (!Plot3D) return;
    const [s, t] = stateAt(r);
    const fig = F.sphereFigure(r, s);
    Plot3D.react(sphereEl, fig.data, fig.layout,
                 { ...F.graphConfig("rp_sphere"), modeBarButtonsToRemove: ["resetCameraLastSave3d", "hoverClosest3d"] });
    sphereDrawn = true;
    spherePath = fig.path;
    pathFaded = false;
    $("#sphere-key-overlap").hidden = !fig.bands.overlap;
    sphereReadout(r, s, t);
    $("#play").disabled = !r.sp_t_monotonic;
    $("#time-scrub").disabled = !r.sp_t_monotonic;
    paint();
  }

  // ---------- Time evolution: hover, scrubber and playback all set sphereTime ----------
  // While "traced" (after Play), the strain (time view) or the amplitude and phase (frequency view)
  // after the current time are faded behind a cursor, and L leaves a trail on the sphere; reaching
  // the end of the signal clears it.

  const player = { playing: false, traced: false, raf: 0, last: 0, follow: null };
  let spherePath = null, pathFaded = false, strainMarked = false, freqMarked = false;

  function bandRange(r) { return [r.sp_t[0], r.sp_t[r.sp_t.length - 1]]; }

  function paint() {
    const r = sphereResult;
    if (!r || !sphereDrawn) return;
    const [s, t] = stateAt(r);
    const mv = F.sphereMoving(s);
    const tr = player.traced && spherePath ? F.sphereTrail(spherePath, s) : { x: [], y: [], z: [] };
    Plot3D.restyle(sphereEl, { x: [...mv.x, tr.x], y: [...mv.y, tr.y], z: [...mv.z, tr.z] },
                   [...F.SPHERE_TRACE, F.TRAIL_TRACE]);
    if (pathFaded !== player.traced) {
      pathFaded = player.traced;
      Plot3D.restyle(sphereEl, { opacity: pathFaded ? 0.3 : 1 }, [F.PATH_TRACE]);
    }
    sphereReadout(r, s, t);
    if (domain === "time") markStrain(r, t); else markFrequency(s);
    const [t0, t1] = bandRange(r), scrub = $("#time-scrub");
    scrub.value = String(Math.round(1000 * Math.min(1, Math.max(0, (t - t0) / (t1 - t0 || 1)))));
    scrub.setAttribute("aria-valuetext", `t − t_c = ${t.toFixed(3)} s`);
  }

  // Cursor and fade on the strain plot while traced; when playing zoomed in, the view scrolls with it.
  function markStrain(r, t) {
    const el = graphs.strain.el;
    if (!r.t || !el.layout || (!player.traced && !strainMarked)) return;
    const shapes = [el.layout.shapes[0]];
    const upd = {};
    if (player.traced) {
      shapes.push({ type: "rect", xref: "x", yref: "paper", x0: t, x1: r.t[r.t.length - 1], y0: 0, y1: 1,
                    fillcolor: "rgba(255,255,255,0.82)", line: { width: 0 }, layer: "above" },
                  { type: "line", xref: "x", yref: "paper", x0: t, x1: t, y0: 0, y1: 1,
                    line: { color: "#1d2433", width: 1.2 } });
      const xr = el._fullLayout.xaxis.range, w = xr[1] - xr[0];
      if (player.playing && player.follow) upd["xaxis.range"] = [t - 0.8 * player.follow, t + 0.2 * player.follow];
      else if (!player.playing && (t < xr[0] || t > xr[1])) upd["xaxis.range"] = [t - 0.8 * w, t + 0.2 * w];  // scrubbed out of view
    }
    strainMarked = player.traced;
    upd.shapes = shapes;
    Plotly.relayout(el, upd);
  }

  // The same cursor and fade on the amplitude and phase plots, at the frequency f(t) of L's state s.
  // The frequency axis stays as it is: the whole band is normally in view.
  function markFrequency(s) {
    if (!player.traced && !freqMarked) return;
    const f = player.traced && s ? s.f : null;
    for (const key of ["amp", "phase"]) {
      const el = graphs[key].el;
      if (!el.layout || !el.layout.shapes) continue;
      const shapes = [el.layout.shapes[0]];
      if (f !== null) {
        shapes.push({ type: "rect", xref: "x", yref: "paper", x0: f, x1: 1e4, y0: 0, y1: 1,
                      fillcolor: "rgba(255,255,255,0.82)", line: { width: 0 }, layer: "above" },
                    { type: "line", xref: "x", yref: "paper", x0: f, x1: f, y0: 0, y1: 1,
                      line: { color: "#1d2433", width: 1.2 } });
      }
      Plotly.relayout(el, { shapes });
    }
    freqMarked = f !== null;
  }

  function setPlayButton() {
    const b = $("#play");
    b.textContent = player.playing ? "❚❚ Pause" : "▶ Play";
    b.setAttribute("aria-label", player.playing ? "Pause the time evolution" : "Play the time evolution");
  }

  function play() {
    const r = sphereResult;
    if (!r || !sphereDrawn || !r.sp_t_monotonic) return;
    const [t0, t1] = bandRange(r);
    if (!player.traced || sphereTime === null || sphereTime < t0 || sphereTime >= t1) sphereTime = t0;
    // Scroll the strain plot with the cursor only when zoomed in; keep a full view still.
    player.follow = null;
    if (domain === "time" && r.t && graphs.strain.el._fullLayout) {
      const xr = graphs.strain.el._fullLayout.xaxis.range, w = xr[1] - xr[0];
      if (w < 0.9 * (r.t[r.t.length - 1] - r.t[0])) player.follow = w;
    }
    player.playing = player.traced = true;
    player.last = performance.now();
    setPlayButton();
    const step = (now) => {
      if (!player.playing) return;
      const r2 = sphereResult;
      if (!r2) { pause(); return; }
      const dt = Math.min(0.1, (now - player.last) / 1000);
      player.last = now;
      const end = bandRange(r2)[1];
      sphereTime = Math.min(end, sphereTime + dt * parseFloat($("#speed").value));
      if (sphereTime >= end) { pause(); player.traced = false; paint(); return; }
      paint();
      player.raf = requestAnimationFrame(step);
    };
    player.raf = requestAnimationFrame(step);
  }

  function pause() {
    player.playing = false;
    cancelAnimationFrame(player.raf);
    setPlayButton();
  }

  // Hover, tap or scrub: move to time t (paused). Throttled to one update per frame.
  function moveSphere(t) {
    sphereTime = t;
    if (!sphereResult || !sphereDrawn || sphereFrame) return;
    sphereFrame = requestAnimationFrame(() => { sphereFrame = 0; paint(); });
  }

  // Hovering or clicking a plot places L there: x is a time on the strain plot and a frequency on
  // the amplitude and phase plots, converted with t(f).
  const linked = {};
  function link(key) {
    if (linked[key]) return;
    linked[key] = true;
    const el = graphs[key].el;
    const pick = (click) => (e) => {
      if (player.playing && !click) return;  // playback owns the time; a click pauses it there
      if (player.playing) pause();
      const x = e && (e.xvals ? e.xvals[0] : e.points && e.points[0] && e.points[0].x);
      if (typeof x !== "number" || !isFinite(x)) return;
      if (key === "strain") moveSphere(x);
      else if (sphereResult) { const t = F.timeAtFrequency(sphereResult, x); if (t !== null) moveSphere(t); }
    };
    el.on("plotly_hover", pick(false));
    el.on("plotly_click", pick(true));
  }

  function notice(kind, title, body) {
    const el = $("#notice");
    el.replaceChildren();
    if (!kind) { el.className = "notice is-empty"; return; }
    el.className = `notice notice-${kind}`;
    const icon = { warn: "!", error: "×", info: "i" }[kind];
    const content = Array.isArray(body[0]) ? h("ul", {}, body.map((item) => h("li", {}, item))) : h("p", {}, body);
    el.append(h("span", { class: "notice-icon", "aria-hidden": "true" }, icon), h("div", {}, h("strong", {}, title), content));
  }

  function sciNodes(v) {
    const exp = Math.floor(Math.log10(Math.abs(v)));
    return exp >= 4 ? [`${F.formatG(v / 10 ** exp, 2)}×10`, h("sup", {}, String(exp))] : [F.formatG(v, 3)];
  }

  function notes(r, time) {
    const out = [];
    let bad = 0;
    const arrays = time ? [r.phi_L, r.cos_L, r.h_RP, r.h_NP] : [r.phi_L, r.cos_L, r.amp, r.phase];
    for (const a of arrays) for (const v of a) if (!isFinite(v)) bad++;
    if (bad) out.push([`${bad.toLocaleString("en-US")} computed values are not finite (NaN or ∞) and are left out of the plots.`]);
    if (F.isPoint(r.phi_L, r.cos_L)) out.push(["L does not move for these parameters, so the modulation of L is a single point."]);
    if (time) return out;  // the remaining notes are about the amplitude and phase panels
    if (F.isFlatZero(r.amp, r.phase)) out.push(["RP coincides with the non-precessing signal here, so its curves lie on the zero line."]);
    const peak = F.ampPeak(r);
    if (peak) {
      out.push(["The non-precessing amplitude B", sub("NP"), " is ≈ 0 for this orientation, so B", sub("RP"), "/B", sub("NP"),
                " − 1 diverges (peak ", ...sciNodes(peak), `). The amplitude axis stops at ±${F.AMP_DISPLAY_LIMIT}; use autoscale in the `
                + "plot toolbar to see the full range."]);
    }
    return out;
  }

  function stat(label, value, unit, note) {
    return h("div", { class: "stat" },
      h("span", { class: "stat-label" }, label),
      h("span", {}, h("span", { class: "stat-value" }, value), unit ? h("span", { class: "stat-unit" }, unit) : null),
      h("span", { class: "stat-note" }, note));
  }

  function stats(d) {
    $("#stats").replaceChildren(
      stat("Total mass", d.total_mass.toFixed(2), unitLabel("M☉"), "source frame"),
      stat("Redshifted mass", d.total_mass_redshifted.toFixed(2), unitLabel("M☉"), ["M", sub("z"), " = M(1+z)"]),
      stat("Cut-off frequency", d.f_cut.toFixed(2), "Hz", ["f", sub("cut"), ", upper end of the band"]),
      stat(["cos ι", sub("JN")], d.cos_i_JN.toFixed(2), "", "angle between J and N"),
      stat("Luminosity distance", d.luminosity_distance.toFixed(3), "Gpc", "flat ΛCDM, Planck 2018"));
  }

  function render(r, time) {
    const exportBtn = $("#export-csv");
    const status = $("#status");
    if (r.invalid) {
      notice("warn", "Check the parameters", [r.invalid]);
      exportBtn.disabled = true;
      return;
    }
    if (r.derived) stats(r.derived);
    if (r.error) {
      blank("No plot for these parameters: see the message above");
      status.textContent = "";
      notice("error", "The model could not be evaluated for these parameters", [r.error]);
      exportBtn.disabled = true;
      return;
    }
    if (r.n === 0) {
      blank(`No frequencies in band: f<sub>cut</sub> = ${r.f_cut.toFixed(2)} Hz < f<sub>min</sub> = ${grid.fMin} Hz`);
      status.textContent = "";
      notice("warn", "No signal in the detector band",
             ["The cut-off frequency f", sub("cut"), ` = ${r.f_cut.toFixed(2)} Hz is below f`, sub("min"),
              ` = ${grid.fMin} Hz, so there are no frequencies to plot. Lower the chirp mass or redshift, or raise the mass ratio.`]);
      exportBtn.disabled = true;
      return;
    }
    status.textContent = time
      ? `${r.nt.toLocaleString("en-US")} samples · ${(r.t[r.nt - 1] - r.t[0]).toFixed(2)} s at ${r.fs} Hz · `
        + `from ${grid.fMin}–${r.f[r.n - 1].toFixed(1)} Hz · ${Math.round(r.ms)} ms`
      : `${r.n.toLocaleString("en-US")} samples · ${grid.fMin}–${r.f[r.n - 1].toFixed(1)} Hz, `
        + `Δf = ${grid.deltaF} Hz · ${Math.round(r.ms)} ms`;
    const n = notes(r, time);
    if (n.length) notice("info", "About these results", n); else notice(null);
    draw("mod", F.modulationFigure(r));
    if (time) {
      draw("strain", F.strainFigure(r));
      link("strain");
      strainMarked = false;
    } else {
      draw("amp", F.amplitudeFigure(r));
      draw("phase", F.phaseFigure(r));
      link("amp");
      link("phase");
      freqMarked = false;
    }
    sphereResult = r;
    if (Plot3D) drawSphere(r); else $("#sphere-readout").textContent = "Loading the 3D view…";
    exportBtn.disabled = false;
  }

  function exportCsv() {
    worker.postMessage({ type: "csv", id: ++runId, values: currentValues(), domain });
  }

  function download([name, content]) {
    const url = URL.createObjectURL(new Blob([content], { type: "text/csv" }));
    const a = h("a", { href: url, download: name });
    document.body.append(a);
    a.click();
    a.remove();
    setTimeout(() => URL.revokeObjectURL(url), 1000);
  }

  // ---------- Worker messages ----------

  worker.onmessage = (e) => {
    const m = e.data;
    if (m.type === "progress") {
      $("#loading-step").textContent = m.text + "…";
    } else if (m.type === "params") {
      P = m.params;
      buildControls();
    } else if (m.type === "ready") {
      grid = m.grid;
      $("#loading").hidden = true;
      enable();
      loadGl3d();
      schedule();
    } else if (m.type === "result") {
      if (!pending) render(m.data, Boolean(m.data.t));  // skip results already superseded
      done();
    } else if (m.type === "csv") {
      if (m.data) download(m.data);
    } else if (m.type === "fatal") {
      $("#loading").hidden = true;
      blank("The model could not be loaded");
      notice("error", "The model could not be loaded in this browser", [m.text]);
    }
  };

  $("#reset-all").addEventListener("click", () => {
    if (P) apply(Object.fromEntries(P.groups.flatMap((g) => g.params.map((p) => [p.id, p.value]))));
  });
  $("#export-csv").addEventListener("click", exportCsv);
  $("#play").addEventListener("click", () => (player.playing ? pause() : play()));
  $("#time-scrub").addEventListener("input", (e) => {
    const r = sphereResult;
    if (!r) return;
    pause();
    player.traced = true;
    const [t0, t1] = bandRange(r);
    moveSphere(t0 + (t1 - t0) * (+e.target.value) / 1000);
  });
  document.querySelectorAll(".domain-switch button").forEach((b) => b.addEventListener("click", () => setDomain(b.dataset.domain)));

  blank();
})();
