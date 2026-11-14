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
    schedule();
  }

  // ---------- Results (webapp/callbacks.py: update) ----------

  function draw(key, fig) {
    const g = graphs[key];
    Plotly.react(g.el, fig.data, fig.layout, F.graphConfig(g.file));
  }

  function blank(message) {
    for (const key of Object.keys(graphs)) draw(key, F.placeholderFigure(message));
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
    } else {
      draw("amp", F.amplitudeFigure(r));
      draw("phase", F.phaseFigure(r));
    }
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
  document.querySelectorAll(".domain-switch button").forEach((b) => b.addEventListener("click", () => setDomain(b.dataset.domain)));

  blank();
})();
