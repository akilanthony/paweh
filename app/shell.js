const notice = document.querySelector("#paweh-loading");
const heading = document.querySelector("#paweh-status");
const message = document.querySelector("#paweh-requirements");
const retry = document.querySelector("#paweh-retry");
const banner = document.querySelector("#paweh-build");
const root = document.querySelector("#root");

function alignNotice() {
  document.body.style.setProperty("--paweh-banner-height", banner.offsetHeight + "px");
}
alignNotice();
window.addEventListener("resize", alignNotice);

function show(title, detail, canRetry = false) {
  if (heading.textContent !== title) heading.textContent = title;
  if (message.textContent !== detail) message.textContent = detail;
  retry.hidden = !canRetry;
}

retry.addEventListener("click", () => location.reload());

function readyIn(frame) {
  try {
    const doc = frame.document;
    if (doc.querySelector("#home-case_control") &&
        frame.Shiny?.shinyapp?.isConnected?.()) return true;
    return Array.from(doc.querySelectorAll("iframe")).some(
      child => child.contentWindow && readyIn(child.contentWindow));
  } catch {
    return false;
  }
}

if (!isSecureContext || !("serviceWorker" in navigator) ||
    typeof WebAssembly !== "object" || typeof Worker !== "function") {
  show("This browser cannot start the dashboard.",
       "Use an up-to-date desktop browser with JavaScript, WebAssembly, and service workers enabled over HTTPS.");
} else {
  const initialMessage = "The first visit downloads about 83 MB and can take several minutes on slower connections. The interactive 3D view requires WebGL.";
  const started = performance.now();
  let failedResource = false;
  show("Preparing the R runtime…", initialMessage);
  window.addEventListener("error", event => {
    const url = event.target?.src || event.target?.href;
    if (url?.startsWith(new URL("./", location.href).href)) failedResource = true;
  }, true);
  function update() {
    if (readyIn(window)) {
      notice.hidden = true;
      root.setAttribute("aria-busy", "false");
      clearInterval(timer);
      return;
    }
    const elapsed = performance.now() - started;
    if (!navigator.onLine) {
      show("You appear to be offline.", "Reconnect to the internet, then retry loading the dashboard.", true);
    } else if (failedResource || elapsed > 180000) {
      show("The dashboard has not finished loading.",
           "Check your connection and browser settings, then retry. A slow first visit may still finish if you wait.", true);
    } else if (elapsed > 45000) {
      show("Still preparing the R workspace…", initialMessage, true);
    } else if (root.querySelector("iframe")) {
      show("Starting the dashboard session…", initialMessage);
    } else {
      show("Preparing the R runtime…", initialMessage);
    }
  }
  const timer = setInterval(update, 250);
  update();
}
