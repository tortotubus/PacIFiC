// docs/js/codewrap.js
(function () {
  function inferLang(pre) {
    const code = pre.querySelector("code");
    if (!code) return "";

    // Common patterns:
    //   <code class="language-html">...
    //   <code class="lang-html">...
    //   <code class="language-html hljs">...
    //   Pygments often doesn't put language on <code>, but may on parent
    const classes = (code.className || "") + " " + (pre.className || "");
    const m =
      classes.match(/\blanguage-([a-z0-9_-]+)\b/i) ||
      classes.match(/\blang-([a-z0-9_-]+)\b/i);

    return m ? m[1] : "";
  }

  function wrapPre(pre) {
    // Avoid double-wrapping
    if (pre.closest(".codewrap")) return;

    const lang = inferLang(pre);
    const code = pre.querySelector("code");

    // Outer wrapper
    const outer = document.createElement("div");
    outer.className = "codewrap bg-dark text-white rounded border mb-4";

    // Header bar
    const header = document.createElement("div");
    header.className = "d-flex align-items-center px-2 py-2 border-bottom";

    const label = document.createElement("small");
    label.className = "font-monospace text-uppercase px-2";
    label.textContent = lang || ""; // leave blank if unknown
    header.appendChild(label);

    const right = document.createElement("div");
    right.className = "d-flex ms-auto";

    const btn = document.createElement("button");
    btn.type = "button";
    btn.className = "btn btn-primary mt-0 me-0 codewrap-copy";
    btn.innerHTML = '<i class="bi bi-clipboard"></i>';
    btn.setAttribute("aria-label", "Copy code");
    btn.setAttribute("title", "Copy");

    right.appendChild(btn);
    header.appendChild(right);

    // Body nesting exactly like your desired HTML
    const body = document.createElement("div");
    const bodyInner1 = document.createElement("div");
    const bodyInner2 = document.createElement("div");

    // Move the <pre> into the new structure
    const parent = pre.parentNode;
    parent.insertBefore(outer, pre);

    outer.appendChild(header);
    outer.appendChild(body);
    body.appendChild(bodyInner1);
    bodyInner1.appendChild(bodyInner2);
    bodyInner2.appendChild(pre);

    // Clipboard behavior
    btn.addEventListener("click", async () => {
      try {
        const text = code ? code.innerText : pre.innerText;
        await navigator.clipboard.writeText(text);

        // quick visual feedback
        btn.classList.add("btn-success");
        btn.classList.remove("btn-primary");
        btn.innerHTML = '<i class="bi bi-check2"></i>';

        window.setTimeout(() => {
          btn.classList.add("btn-primary");
          btn.classList.remove("btn-success");
          btn.innerHTML = '<i class="bi bi-clipboard"></i>';
        }, 900);
      } catch (e) {
        // fallback feedback (no alerts)
        btn.classList.add("btn-danger");
        btn.classList.remove("btn-primary");
        btn.innerHTML = '<i class="bi bi-exclamation-triangle"></i>';

        window.setTimeout(() => {
          btn.classList.add("btn-primary");
          btn.classList.remove("btn-danger");
          btn.innerHTML = '<i class="bi bi-clipboard"></i>';
        }, 1200);
      }
    });
  }

  function run() {
    // Wrap plain Markdown code blocks:
    document.querySelectorAll("pre").forEach(wrapPre);

    // If you're using Pygments, you may get <div class="highlight"><pre>...</pre></div>
    // The loop above still catches the inner <pre>, so this is usually enough.
  }

  if (document.readyState === "loading") {
    document.addEventListener("DOMContentLoaded", run);
  } else {
    run();
  }
})();
