/* Client-side helpers for the Regular Precession explorer (display formatting only). */
(function () {
  function gcd(a, b) { return b ? gcd(b, a % b) : a; }

  // Format an angle in radians as a multiple of π (sliders step in π/16) plus its decimal value.
  function formatAngle(v) {
    if (v === undefined || v === null || isNaN(v)) { return "—"; }
    var n = Math.round(v * 16 / Math.PI);
    if (Math.abs(n * Math.PI / 16 - v) > 1e-9) { return v.toFixed(3) + " rad"; }
    if (n === 0) { return "0"; }
    var g = gcd(Math.abs(n), 16), num = n / g, den = 16 / g;
    var frac = (num === 1 ? "" : num) + "π" + (den === 1 ? "" : "/" + den);
    return frac + " · " + v.toFixed(3);
  }

  // Dash slider number boxes commit on blur; make Enter commit too.
  document.addEventListener("keydown", function (e) {
    if (e.key === "Enter" && e.target.matches && e.target.matches(".slider input[type=number]")) {
      e.target.blur();
    }
  });

  window.dash_clientside = Object.assign({}, window.dash_clientside, {
    rp: {
      angle: function (dragValue, value) {
        var ctx = window.dash_clientside.callback_context;
        var fromValue = ctx && ctx.triggered && ctx.triggered.length &&
          /\.value$/.test(ctx.triggered[0].prop_id);
        var v = (fromValue || dragValue === undefined || dragValue === null) ? value : dragValue;
        return formatAngle(v);
      }
    }
  });
})();
