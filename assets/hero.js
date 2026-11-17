/* Hero animation: image -> latent representation -> GWAS.
   A diffusion-MRI axial slice (DTI-coloured glyphs) is swept by a scan line;
   sampled voxels flow through a latent bottleneck and land as SNPs on a
   Manhattan plot. Loops forever; renders a single still frame when the user
   prefers reduced motion. */
(function () {
  "use strict";

  var canvas = document.getElementById("hero-canvas");
  if (!canvas || !canvas.getContext) return;
  var ctx = canvas.getContext("2d");
  var reduceMotion = window.matchMedia("(prefers-reduced-motion: reduce)").matches;

  var C = {
    lr: [255, 93, 115],
    ap: [62, 230, 160],
    si: [90, 168, 255],
    gm: [150, 170, 205],
    chrA: [90, 168, 255],
    chrB: [184, 163, 255],
    hit: [62, 230, 160],
    text: "rgba(159,176,204,0.85)",
    line: "rgba(150,180,230,0.22)",
  };
  var rgba = function (c, a) { return "rgba(" + c[0] + "," + c[1] + "," + c[2] + "," + a + ")"; };

  /* deterministic RNG so the brain and plot look the same on every load */
  var seed = 20261117;
  function rnd() { seed = (seed * 16807) % 2147483647; return (seed - 1) / 2147483646; }

  /* ---------- brain slice (unit coordinates, x: L-R, y: A-P) ---------- */
  var glyphs = [];
  function rim(t) { return 1 + 0.035 * Math.sin(13 * t) + 0.025 * Math.sin(29 * t + 1) + 0.02 * Math.sin(7 * t + 2); }
  function smooth(e0, e1, v) { var t = Math.min(1, Math.max(0, (v - e0) / (e1 - e0))); return t * t * (3 - 2 * t); }
  function gauss2(x, y, cx, cy, sx, sy) { return Math.exp(-Math.pow((x - cx) / sx, 2) - Math.pow((y - cy) / sy, 2)); }
  (function buildBrain() {
    var target = 2300, guard = 0;
    while (glyphs.length < target && guard++ < 60000) {
      var x = rnd() * 2 - 1, y = rnd() * 2 - 1;
      var ex = x / 0.8, ey = y;
      var t = Math.atan2(ey, ex), rr = Math.hypot(ex, ey);
      if (rr > rim(t)) continue;
      if (Math.abs(x) < 0.025 && Math.abs(y) > 0.62) continue;                     // interhemispheric fissure
      if (Math.hypot((Math.abs(x) - 0.11) / 0.06, (y - 0.02) / 0.26) < 1) continue; // lateral ventricles
      var ax = Math.abs(x), sx = x < 0 ? -1 : 1, sy = y < 0 ? -1 : 1;
      // fibre field: callosal forceps (L-R, curving), association fibres (A-P), corona radiata (S-I)
      var wcc = gauss2(ax, Math.abs(y), 0, 0.45, 0.5, 0.15);
      var aCC = sx * sy * 1.1 * ax;                                  // forceps curve away from midline
      var aAP = Math.PI / 2 + 0.55 * x - 0.25 * sx * y;
      var vx = wcc * Math.cos(aCC) + (1 - wcc) * Math.cos(aAP);
      var vy = wcc * Math.sin(aCC) + (1 - wcc) * Math.sin(aAP);
      var z = Math.min(1, gauss2(ax, y, 0.3, 0.02, 0.1, 0.22) * 1.25);  // through-plane
      var cortex = smooth(0.74, 0.9, rr);
      var fa = (0.35 + 0.55 * Math.max(wcc, z, 0.55)) * (1 - 0.75 * cortex) + (rnd() - 0.5) * 0.12;
      var ang = Math.atan2(vy, vx);
      if (cortex > 0.5) ang = t + (rnd() - 0.5) * 0.6;              // grey matter: weak, radial
      var cx = Math.abs(Math.cos(ang)) * (1 - z), cy = Math.abs(Math.sin(ang)) * (1 - z);
      var col = [0, 1, 2].map(function (k) {
        var v = cx * C.lr[k] + cy * C.ap[k] + z * C.si[k];
        return Math.round(v / Math.max(cx + cy + z, 0.001) * (1 - cortex) + C.gm[k] * cortex);
      });
      glyphs.push({ x: x, y: y, ang: ang, len: (1 - 0.85 * z) * (1 - 0.5 * cortex), fa: Math.max(0.12, fa), col: col });
    }
  })();

  /* ---------- GWAS: chromosomes, SNPs, peaks ---------- */
  var chromLen = [248, 242, 198, 190, 181, 171, 159, 145, 138, 134, 135, 133, 114, 107, 102, 90, 83, 80, 59, 64, 47, 51];
  var genomeLen = chromLen.reduce(function (a, b) { return a + b; }, 0);
  var chromStart = [];
  chromLen.reduce(function (acc, l, i) { chromStart[i] = acc; return acc + l; }, 0);
  var peaks = [[1, 0.42, 11.5], [5, 0.7, 8.6], [9, 0.3, 13.2], [13, 0.55, 9.4], [16, 0.48, 7.9], [3, 0.2, 6.2]];
  var Y_MAX = 14, SIG = 7.3;
  var snps = [];
  (function buildSnps() {
    chromLen.forEach(function (l, ci) {
      var n = Math.round(l * 3.4);
      for (var k = 0; k < n; k++) {
        var u = rnd();
        snps.push({ g: (chromStart[ci] + u * l) / genomeLen, p: Math.min(-Math.log10(Math.max(rnd(), 1e-9)) * 0.9, 5.2), c: ci });
      }
    });
    peaks.forEach(function (pk) {
      var ci = pk[0], center = pk[1], top = pk[2];
      for (var k = 0; k < 34; k++) {
        var d = (rnd() - 0.5) * 0.12;
        var h = top * Math.exp(-Math.pow(d / 0.03, 2)) * (0.75 + rnd() * 0.25);
        snps.push({ g: (chromStart[ci] + (center + d) * chromLen[ci]) / genomeLen, p: Math.max(h, 1), c: ci });
      }
    });
    // shuffle so landing order is spread across the genome
    for (var i = snps.length - 1; i > 0; i--) { var j = Math.floor(rnd() * (i + 1)); var tmp = snps[i]; snps[i] = snps[j]; snps[j] = tmp; }
  })();

  /* ---------- layout ---------- */
  var W = 0, H = 0, L = null;
  function layout() {
    var dpr = Math.min(window.devicePixelRatio || 1, 2);
    var r = canvas.getBoundingClientRect();
    W = r.width; H = r.height;
    canvas.width = Math.round(W * dpr); canvas.height = Math.round(H * dpr);
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    var nodes = [], K = 10, i;
    if (W / H > 1.05) {
      var s = Math.min(W * 0.2, H * 0.4);
      L = { bx: W * 0.2, by: H * 0.47, s: s, plot: { x0: W * 0.55, x1: W * 0.985, y0: H * 0.84, y1: H * 0.14 }, vertical: false };
      for (i = 0; i < K; i++) nodes.push({ x: W * 0.425, y: H * 0.47 + (i - (K - 1) / 2) * Math.min(18, H * 0.03) });
    } else {
      var s2 = Math.min(W * 0.3, H * 0.19);
      L = { bx: W * 0.5, by: H * 0.22, s: s2, plot: { x0: W * 0.08, x1: W * 0.98, y0: H * 0.93, y1: H * 0.6 }, vertical: true };
      for (i = 0; i < K; i++) nodes.push({ x: W * 0.5 + (i - (K - 1) / 2) * Math.min(16, W * 0.04), y: H * 0.475 });
    }
    L.nodes = nodes;
  }
  function brainPos(g) { return { x: L.bx + g.x * L.s, y: L.by + g.y * L.s * 1.08 }; }
  function snpPos(sn) {
    var P = L.plot;
    return { x: P.x0 + sn.g * (P.x1 - P.x0), y: P.y0 - (sn.p / Y_MAX) * (P.y0 - P.y1) };
  }

  /* ---------- state ---------- */
  var particles = [];
  var nextTarget = 0;         // next SNP index to send a particle to
  var landedCount = 0;
  var landed = new Uint8Array(snps.length);
  var flash = new Float32Array(snps.length);
  var scanY = -1.1;           // unit coords
  var phase = "run", phaseT = 0, fade = 1;

  function spawn(dt) {
    var want = Math.min(snps.length - nextTarget, Math.ceil(dt * 0.11));
    for (var k = 0; k < want; k++) {
      // pick a glyph near the scan line
      var g = null, tries = 0;
      while (tries++ < 40) {
        var cand = glyphs[Math.floor(Math.random() * glyphs.length)];
        if (Math.abs(cand.y - scanY) < 0.12) { g = cand; break; }
      }
      if (!g) g = glyphs[Math.floor(Math.random() * glyphs.length)];
      g.lit = 1;
      particles.push({
        from: brainPos(g),
        node: L.nodes[Math.floor(Math.random() * L.nodes.length)],
        target: nextTarget++,
        col: g.col,
        t: 0,
        speed: 0.0006 + Math.random() * 0.00035,
      });
    }
  }

  function bez(a, c, b, t) {
    var u = 1 - t;
    return { x: u * u * a.x + 2 * u * t * c.x + t * t * b.x, y: u * u * a.y + 2 * u * t * c.y + t * t * b.y };
  }
  function ease(t) { return t < 0.5 ? 2 * t * t : 1 - Math.pow(-2 * t + 2, 2) / 2; }

  /* ---------- drawing ---------- */
  function drawBrain(alpha) {
    var unit = L.s / 78;
    ctx.lineCap = "round";
    for (var i = 0; i < glyphs.length; i++) {
      var g = glyphs[i], p = brainPos(g);
      var scan = Math.max(0, 1 - Math.abs(g.y - scanY) / 0.08);
      var lit = g.lit || 0;
      var a = Math.min(1, (g.fa * 0.75 + scan * 0.6 + lit * 0.6)) * alpha;
      if (g.len < 0.25) {
        ctx.fillStyle = rgba(g.col, a);
        ctx.beginPath(); ctx.arc(p.x, p.y, 1.5 * unit + scan * unit, 0, 6.283); ctx.fill();
      } else {
        var half = (2.8 * g.len + scan * 1.2) * unit;
        var dx = Math.cos(g.ang) * half, dy = Math.sin(g.ang) * half;
        ctx.strokeStyle = rgba(g.col, a);
        ctx.lineWidth = 1.3 * unit;
        ctx.beginPath(); ctx.moveTo(p.x - dx, p.y - dy); ctx.lineTo(p.x + dx, p.y + dy); ctx.stroke();
      }
      if (g.lit) g.lit = Math.max(0, g.lit - 0.03);
    }
    // scan line
    if (scanY > -1.05 && scanY < 1.05) {
      var sy = L.by + scanY * L.s * 1.08;
      var grad = ctx.createLinearGradient(L.bx - L.s, 0, L.bx + L.s, 0);
      grad.addColorStop(0, "rgba(233,239,250,0)");
      grad.addColorStop(0.5, "rgba(233,239,250," + 0.55 * alpha + ")");
      grad.addColorStop(1, "rgba(233,239,250,0)");
      ctx.strokeStyle = grad; ctx.lineWidth = 1;
      ctx.beginPath(); ctx.moveTo(L.bx - L.s * 0.95, sy); ctx.lineTo(L.bx + L.s * 0.95, sy); ctx.stroke();
    }
  }

  function drawNodes(t) {
    var n = L.nodes;
    ctx.strokeStyle = C.line; ctx.lineWidth = 1;
    ctx.beginPath();
    if (L.vertical) { ctx.moveTo(n[0].x - 10, n[0].y); ctx.lineTo(n[n.length - 1].x + 10, n[0].y); }
    else { ctx.moveTo(n[0].x, n[0].y - 10); ctx.lineTo(n[0].x, n[n.length - 1].y + 10); }
    ctx.stroke();
    for (var i = 0; i < n.length; i++) {
      var pulse = 0.5 + 0.5 * Math.sin(t * 0.003 + i * 0.9);
      ctx.fillStyle = "rgba(184,163,255," + (0.35 + pulse * 0.5) + ")";
      ctx.beginPath(); ctx.arc(n[i].x, n[i].y, 3 + pulse * 1.2, 0, 6.283); ctx.fill();
    }
  }

  function drawAxes() {
    var P = L.plot;
    ctx.strokeStyle = C.line; ctx.lineWidth = 1;
    ctx.beginPath(); ctx.moveTo(P.x0, P.y1); ctx.lineTo(P.x0, P.y0); ctx.lineTo(P.x1, P.y0); ctx.stroke();
    // genome-wide significance
    var sy = P.y0 - (SIG / Y_MAX) * (P.y0 - P.y1);
    ctx.setLineDash([4, 5]); ctx.strokeStyle = "rgba(255,93,115,0.5)";
    ctx.beginPath(); ctx.moveTo(P.x0, sy); ctx.lineTo(P.x1, sy); ctx.stroke(); ctx.setLineDash([]);
    ctx.fillStyle = C.text;
    ctx.font = "11px 'IBM Plex Sans', sans-serif";
    ctx.textAlign = "center";
    var step = (P.x1 - P.x0) < 380 ? 3 : 1;
    for (var ci = 0; ci < chromLen.length; ci += step) {
      var mid = P.x0 + ((chromStart[ci] + chromLen[ci] / 2) / genomeLen) * (P.x1 - P.x0);
      ctx.fillText(String(ci + 1), mid, P.y0 + 16);
    }
    ctx.save(); ctx.translate(P.x0 - 12, (P.y0 + P.y1) / 2); ctx.rotate(-Math.PI / 2);
    ctx.fillText("−log₁₀ p", 0, 0); ctx.restore();
  }

  function drawSnps(alpha) {
    for (var i = 0; i < snps.length; i++) {
      if (!landed[i]) continue;
      var sn = snps[i], p = snpPos(sn);
      var sig = sn.p >= SIG;
      var base = sig ? C.hit : (sn.c % 2 ? C.chrB : C.chrA);
      var f = flash[i];
      ctx.fillStyle = rgba(f > 0.05 ? [233, 239, 250] : base, (sig ? 0.95 : 0.62) * alpha + f * 0.3);
      ctx.beginPath(); ctx.arc(p.x, p.y, (sig ? 2.4 : 1.7) + f * 2.2, 0, 6.283); ctx.fill();
      if (f > 0) flash[i] = Math.max(0, f - 0.04);
    }
  }

  function drawParticles(dt) {
    for (var i = particles.length - 1; i >= 0; i--) {
      var pt = particles[i];
      pt.t += dt * pt.speed;
      var to = snpPos(snps[pt.target]);
      var pos, col = pt.col;
      if (pt.t < 0.5) {
        var a = ease(pt.t * 2);
        var c1 = L.vertical ? { x: pt.from.x, y: (pt.from.y + pt.node.y) / 2 } : { x: (pt.from.x + pt.node.x) / 2, y: pt.from.y };
        pos = bez(pt.from, c1, pt.node, a);
      } else {
        var b = ease((pt.t - 0.5) * 2);
        var c2 = L.vertical ? { x: pt.node.x, y: (pt.node.y + to.y) / 2 } : { x: (pt.node.x + to.x) / 2, y: pt.node.y };
        pos = bez(pt.node, c2, to, b);
        col = [Math.round(col[0] + (184 - col[0]) * b), Math.round(col[1] + (163 - col[1]) * b), Math.round(col[2] + (255 - col[2]) * b)];
      }
      ctx.fillStyle = rgba(col, 0.9);
      ctx.beginPath(); ctx.arc(pos.x, pos.y, 1.6, 0, 6.283); ctx.fill();
      if (pt.t >= 1) {
        particles.splice(i, 1);
        landed[pt.target] = 1; landedCount++;
        flash[pt.target] = 1;
      }
    }
  }

  /* ---------- loop ---------- */
  var last = 0, running = false, rafId = 0, clock = 0;
  function frame(now) {
    var dt = last ? Math.min(now - last, 50) : 16;
    last = now; clock += dt;
    ctx.clearRect(0, 0, W, H);

    if (phase === "run") {
      scanY += dt * 0.00042;
      if (scanY > 1.1) scanY = -1.1;
      spawn(dt);
      if (landedCount >= snps.length && particles.length === 0) { phase = "hold"; phaseT = 0; }
    } else if (phase === "hold") {
      phaseT += dt;
      scanY += dt * 0.00042; if (scanY > 1.1) scanY = -1.1;
      if (phaseT > 3200) { phase = "fade"; phaseT = 0; }
    } else if (phase === "fade") {
      phaseT += dt; fade = Math.max(0, 1 - phaseT / 900);
      if (fade === 0) { nextTarget = 0; landedCount = 0; landed.fill(0); flash.fill(0); fade = 1; phase = "run"; }
    }

    drawBrain(1);
    drawNodes(clock);
    drawAxes();
    drawSnps(phase === "fade" ? fade : 1);
    drawParticles(dt);
    if (running) rafId = requestAnimationFrame(frame);
  }

  function still() {
    landed.fill(1); landedCount = snps.length; scanY = 2;
    ctx.clearRect(0, 0, W, H);
    drawBrain(1); drawNodes(0); drawAxes(); drawSnps(1);
  }

  function start() { if (running || reduceMotion) return; running = true; last = 0; rafId = requestAnimationFrame(frame); }
  function stop() { running = false; cancelAnimationFrame(rafId); }

  layout();
  if (reduceMotion) still();

  var resizeTimer = 0;
  window.addEventListener("resize", function () {
    clearTimeout(resizeTimer);
    resizeTimer = setTimeout(function () { layout(); if (reduceMotion) still(); }, 120);
  });

  if ("IntersectionObserver" in window) {
    new IntersectionObserver(function (entries) {
      entries.forEach(function (e) { e.isIntersecting ? start() : stop(); });
    }).observe(canvas);
  } else {
    start();
  }
  document.addEventListener("visibilitychange", function () { document.hidden ? stop() : start(); });
})();
