/* Runs the Python model in Pyodide, off the main thread so the page stays responsive.

   Messages in:  {type: "run", id, values, domain} | {type: "csv", id, values, domain}
   Messages out: {type: "progress", text} | {type: "params", params} | {type: "ready", grid}
                 {type: "result", id, data} | {type: "csv", id, data} | {type: "fatal", text}
*/
const PYODIDE = "https://cdn.jsdelivr.net/pyodide/v0.29.0/full/";

// Same paths as in the repository; build.sh copies them under py/.
const SOURCES = [
  "regular_precession_no_match.py",
  "systems_lib.py",
  "webapp/__init__.py",
  "webapp/params.py",
  "webapp/physics.py",
  "webapp/time_domain.py",
  "webapp/sphere.py",
  "rp_page.py",
];

importScripts(PYODIDE + "pyodide.js");

let bridge = null;

function f64(bytes) {
  // Copy into an aligned buffer: the view Pyodide hands over may start at any byte offset.
  return new Float64Array(bytes.slice().buffer);
}

async function fetchText(url, attempts = 3) {
  for (let i = 1; ; i++) {
    try {
      const res = await fetch(url);
      if (!res.ok) throw new Error(`HTTP ${res.status}`);
      return await res.text();
    } catch (err) {
      if (i >= attempts) throw new Error(`Could not download ${url} (${err.message})`);
      await new Promise((r) => setTimeout(r, 300 * i));
    }
  }
}

async function init() {
  const progress = (text) => postMessage({ type: "progress", text });

  progress("Starting Python");
  const pyodide = await loadPyodide({ indexURL: PYODIDE });

  progress("Loading the model code");
  pyodide.FS.mkdirTree("/home/pyodide/webapp");
  // One file at a time, with retries: simple local servers refuse bursts of parallel connections.
  for (const path of SOURCES) {
    pyodide.FS.writeFile("/home/pyodide/" + path, await fetchText("py/" + path));
  }

  // params.py needs only NumPy, so the controls can appear before SciPy and Astropy arrive.
  progress("Loading NumPy");
  await pyodide.loadPackage("numpy");
  pyodide.runPython("import sys; sys.path.insert(0, '/home/pyodide')");
  bridge = pyodide.pyimport("rp_page");
  postMessage({ type: "params", params: JSON.parse(bridge.params_json()) });

  progress("Loading SciPy and Astropy");
  await pyodide.loadPackage(["scipy", "astropy"]);
  progress("Preparing the model");
  postMessage({ type: "ready", grid: JSON.parse(bridge.load_model()) });
}

function run(values, domain) {
  const res = bridge.evaluate(values, domain).toJs({ dict_converter: Object.fromEntries });
  for (const key of ["f", "phi_L", "cos_L", "amp", "phase", "t", "h_RP", "h_NP", "sp_t", "sp_theta", "sp_phi"]) {
    if (res[key]) res[key] = f64(res[key]);
  }
  return res;
}

const ready = init().catch((err) => {
  postMessage({ type: "fatal", text: String(err && err.message || err) });
  throw err;
});

onmessage = async (e) => {
  const { type, id, values, domain } = e.data;
  await ready;
  try {
    if (type === "run") {
      const t0 = performance.now();
      const data = run(values, domain);
      data.ms = performance.now() - t0;
      postMessage({ type: "result", id, data });
    } else if (type === "csv") {
      const out = bridge.export_csv(values, domain);
      postMessage({ type: "csv", id, data: out ? out.toJs() : null });
    }
  } catch (err) {
    postMessage({ type: "result", id, data: { error: String(err && err.message || err) } });
  }
};
