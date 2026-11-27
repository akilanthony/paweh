const notice = document.querySelector("#paweh-loading");
const status = document.querySelector("#paweh-status");
const requirements = document.querySelector("#paweh-requirements");

if (!isSecureContext || !("serviceWorker" in navigator) ||
    typeof WebAssembly !== "object" || typeof Worker !== "function") {
  status.textContent = "This browser cannot start the dashboard.";
  requirements.textContent = "Use a browser with WebAssembly, workers, and service workers over HTTPS or localhost.";
} else {
  function readyIn(frame) {
    try {
      const doc = frame.document;
      if (doc.querySelector("#home-case_control") &&
          frame.Shiny?.shinyapp?.isConnected?.()) return true;
      return Array.from(doc.querySelectorAll("iframe")).some(
        child => child.contentWindow && readyIn(child.contentWindow)
      );
    } catch {
      return false;
    }
  }
  const started = performance.now();
  const timer = setInterval(() => {
    if (readyIn(window)) {
      notice.hidden = true;
      document.querySelector("#root").setAttribute("aria-busy", "false");
      clearInterval(timer);
    } else if (performance.now() - started > 45000) {
      status.textContent = "Still preparing the R workspace…";
      requirements.textContent = "The first visit downloads the R runtime and several packages. Slow connections can take a few minutes.";
    }
  }, 250);
}
