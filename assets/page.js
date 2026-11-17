/* Page behaviour: nav state, agent log typing, publication filters. */
(function () {
  "use strict";
  var reduceMotion = window.matchMedia("(prefers-reduced-motion: reduce)").matches;

  /* nav background once scrolled */
  var nav = document.getElementById("nav");
  function onScroll() { nav.classList.toggle("is-scrolled", window.scrollY > 24); }
  window.addEventListener("scroll", onScroll, { passive: true });
  onScroll();

  /* agent log: types the pipeline once */
  var LOG = [
    ["cmd", "$ agent run imaging-genetics"],
    ["", "<span class='k'>load</span>    diffusion MRI → FA maps"],
    ["", "<span class='k'>encode</span>  images → latent endophenotypes"],
    ["", "<span class='k'>gwas</span>    torchgwas --phenotypes all --gpu"],
    ["", "<span class='k'>explain</span> loci → genes → cell types"],
    ["ok", "✓ genetic architecture mapped"],
  ];
  var logEl = document.getElementById("agent-log");
  if (logEl) {
    if (reduceMotion) {
      logEl.innerHTML = LOG.map(function (l) { return "<div class='ln " + l[0] + "'>" + l[1] + "</div>"; }).join("");
    } else {
      typeLines(logEl, LOG);
    }
  }

  function typeLines(el, lines) {
    var li = 0;
    var caret = document.createElement("span");
    caret.className = "caret";
    function nextLine() {
      if (li >= lines.length) { el.lastChild.appendChild(caret); return; }
      var row = document.createElement("div");
      row.className = "ln " + lines[li][0];
      el.appendChild(row);
      // split the line's HTML into tags (inserted whole) and characters (typed)
      var tokens = lines[li][1].match(/<[^>]+>[^<]*<\/[^>]+>|./g) || [];
      var ti = 0;
      (function step() {
        if (ti < tokens.length) {
          var tok = tokens[ti++];
          row.insertAdjacentHTML("beforeend", tok);
          row.appendChild(caret);
          setTimeout(step, tok.length > 1 ? 40 : 18);
        } else {
          li++;
          setTimeout(nextLine, li === 1 ? 500 : 380);
        }
      })();
    }
    setTimeout(nextLine, 400);
  }

  /* publications: filter + search */
  var list = document.getElementById("pub-list");
  if (!list) return;
  var pubs = Array.prototype.slice.call(list.querySelectorAll(".pub"));
  var groups = Array.prototype.slice.call(list.querySelectorAll(".year-group"));
  var buttons = Array.prototype.slice.call(document.querySelectorAll(".seg button"));
  var search = document.getElementById("pub-search");
  var count = document.getElementById("pub-count");
  var empty = document.getElementById("pub-empty");
  var filter = "all";

  function apply() {
    var q = (search.value || "").trim().toLowerCase();
    var shown = 0;
    pubs.forEach(function (p) {
      var ok = (filter === "all")
        || (filter === "first" && p.dataset.first === "true")
        || (filter === "preprint" && p.dataset.kind === "preprint");
      if (ok && q) ok = p.dataset.search.indexOf(q) !== -1;
      p.hidden = !ok;
      if (ok) shown++;
    });
    groups.forEach(function (g) { g.hidden = !g.querySelector(".pub:not([hidden])"); });
    count.textContent = shown + (shown === 1 ? " paper" : " papers");
    empty.hidden = shown !== 0;
  }

  buttons.forEach(function (b) {
    b.addEventListener("click", function () {
      filter = b.dataset.filter;
      buttons.forEach(function (x) { x.setAttribute("aria-pressed", String(x === b)); });
      apply();
    });
  });
  search.addEventListener("input", apply);
  apply();
})();
