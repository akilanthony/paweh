# Shared presentation helpers for the paweh dashboard.

.paweh_dashboard_theme <- function() {
  theme <- bslib::bs_theme(
    version = 5,
    bg = "#F7F8FA",
    fg = "#20252B",
    primary = "#355C7D",
    secondary = "#6F879A",
    base_font = bslib::font_collection(
      "system-ui", "-apple-system", "BlinkMacSystemFont", "Segoe UI",
      "sans-serif"
    )
  )
  bslib::bs_add_rules(
    theme,
    "
    :root {
      --paweh-bg: #F7F8FA;
      --paweh-panel: #FFFFFF;
      --paweh-text: #20252B;
      --paweh-secondary: #5F6B76;
    --paweh-muted: #66727D;
      --paweh-border: #D9DEE3;
      --paweh-divider: #E9ECEF;
      --paweh-primary: #355C7D;
      --paweh-caution: #A8844F;
      --paweh-error: #9A4D4D;
    }
    body { font-size: .925rem; color: var(--paweh-text); }
    .bslib-page-navbar > .navbar {
      min-height: 49px; padding-top: 0; padding-bottom: 0;
      background: #FFFFFF; border-bottom: 1px solid var(--paweh-border);
      box-shadow: none;
    }
    .navbar-brand { font-size: 1.2rem; font-weight: 500; }
    .navbar .nav-link, .navbar .navbar-nav > li > a {
      color: var(--paweh-secondary); padding: .82rem .75rem .68rem;
      border-bottom: 2px solid transparent; font-size: .925rem;
    }
    .navbar .nav-link.active, .navbar .navbar-nav > li.active > a {
      color: var(--paweh-text); background: transparent;
      border-bottom-color: var(--paweh-primary);
    }
    body.bslib-page-navbar > .container-fluid { padding: 0 20px 20px; }
    .paweh-page-heading { margin: 20px 0 22px; }
    .paweh-page-heading h2 {
      margin: 0; font-size: 1.75rem; line-height: 1.2; font-weight: 500;
    }
    .paweh-page-heading p {
      margin: 7px 0 0; color: var(--paweh-secondary); max-width: 700px;
      font-size: .975rem;
    }
    .paweh-home-hero { margin: 24px 0; padding: 0; }
    .paweh-home-hero h1 { margin: 0 0 5px; font-size: 1.8rem; font-weight: 500; }
    .paweh-home-hero .lead { margin-bottom: 6px; font-size: 1rem; color: var(--paweh-text); }
    .paweh-home-hero p:last-child { margin-bottom: 0; }
    .paweh-home-grid { margin-bottom: 24px; }
    .card {
      background: var(--paweh-panel); border-color: var(--paweh-border);
      border-radius: 7px; box-shadow: none;
    }
    .paweh-study-card { min-height: 0; height: auto; }
    .paweh-study-card .card-header {
      padding: 16px 18px 7px; border-bottom: 0; background: #FFFFFF;
    }
    .paweh-study-card .card-header h3 { font-size: 1.05rem; font-weight: 500; }
    .paweh-study-card .card-body {
      display: flex; flex-direction: column; gap: 12px; padding: 7px 18px 18px;
    }
    .paweh-study-card p { margin: 0; color: var(--paweh-secondary); }
    .paweh-study-card .btn { margin: 2px 0 0; align-self: flex-start; }
    .btn { border-radius: 5px; font-size: .9rem; }
    .btn-primary { background: var(--paweh-primary); border-color: var(--paweh-primary); }
    .btn:focus-visible, .form-control:focus-visible, .form-select:focus-visible,
    .nav-link:focus-visible, summary:focus-visible, input:focus-visible {
      outline: 3px solid rgba(53, 92, 125, .35); outline-offset: 2px;
    }
    .paweh-workspace .bslib-sidebar-layout > .sidebar {
      width: 280px; background: #FFFFFF; border: 1px solid var(--paweh-border);
      border-radius: 7px 0 0 7px;
    }
    .paweh-workspace .bslib-sidebar-layout > .sidebar > .sidebar-content { padding: 16px; }
    .paweh-workspace .bslib-sidebar-layout > .sidebar .sidebar-title {
      margin-bottom: 12px; font-size: 1.05rem; font-weight: 500;
    }
    .paweh-workspace .form-group, .paweh-workspace .shiny-input-container {
      margin-bottom: 10px;
    }
    .paweh-workspace label, .paweh-workspace .control-label {
      margin-bottom: 4px; color: #3F4850; font-size: .91rem; font-weight: 500;
    }
    .paweh-workspace .form-control, .paweh-workspace .form-select {
      min-height: 36px; padding-top: 5px; padding-bottom: 5px; font-size: .9rem;
    }
    .paweh-workspace .radio, .paweh-workspace .checkbox { margin-bottom: 4px; }
    .paweh-sidebar-section {
      margin: 4px 0 12px; padding: 10px 0; border-top: 1px solid var(--paweh-divider);
    }
    .paweh-sidebar-section summary { cursor: pointer; font-size: .93rem; }
    .paweh-sidebar-section summary::after {
      content: 'Misclassification, heterogeneity, and related assumptions';
      display: block; margin-top: 3px; color: var(--paweh-muted);
      font-size: .8rem; font-weight: 400;
    }
    .paweh-qtl-workspace .paweh-sidebar-section summary::after {
      content: 'Test and sampling options';
    }
    .paweh-sidebar-section[open] summary { margin-bottom: 10px; }
    .paweh-sidebar-section p, .paweh-placeholder p { color: var(--paweh-secondary); }
    .paweh-calculate { width: 100%; min-height: 40px; margin-top: 2px; }
    .paweh-workspace .bslib-sidebar-layout > .main {
      padding: 0 0 0 12px; background: transparent; border: 0;
    }
    .paweh-workspace, .paweh-workspace .bslib-sidebar-layout,
    .paweh-workspace .bslib-sidebar-layout > .main,
    .paweh-workspace .bslib-card {
      min-height: 0 !important; height: auto !important; flex: 0 0 auto;
    }
    .paweh-workspace .card-header-tabs { padding: 0 16px; background: #FFFFFF; }
    .paweh-workspace .nav-tabs { border-bottom-color: var(--paweh-border); }
    .paweh-workspace .nav-tabs .nav-link,
    .paweh-workspace .nav-tabs > li > a {
      padding: 11px 12px 9px; color: var(--paweh-secondary); font-size: .92rem;
      border: 0; border-bottom: 2px solid transparent; border-radius: 0;
    }
    .paweh-workspace .nav-tabs .nav-link.active,
    .paweh-workspace .nav-tabs > li.active > a {
      color: var(--paweh-text); border-bottom-color: var(--paweh-primary);
      background: #FFFFFF;
    }
    .paweh-workspace .card-body { padding: 16px; }
    .paweh-workspace .html-widget-output { max-width: 100%; overflow: hidden; }
    .paweh-placeholder { max-width: 560px; padding: 16px 8px; }
    .paweh-placeholder h3 { margin-bottom: 6px; font-size: 1.08rem; font-weight: 500; }
    .paweh-placeholder p { margin-bottom: 0; }
    .paweh-design-summary { margin-top: 12px; }
    .paweh-design-summary .card-header { padding: 10px 12px; font-weight: 500; }
    .paweh-design-summary .card-body { padding: 8px 12px 10px; }
    .paweh-design-summary ul { margin: 0; padding: 0; list-style: none; }
    .paweh-design-summary li { padding: 5px 0; border-bottom: 1px solid var(--paweh-divider); font-size: .84rem; }
    .paweh-design-summary li:last-child { border-bottom: 0; }
    .paweh-summary-grid { display: grid; gap: 0; }
    .paweh-summary-row {
      display: grid; grid-template-columns: minmax(72px, .8fr) minmax(0, 1.2fr);
      gap: 8px; padding: 5px 0; border-bottom: 1px solid var(--paweh-divider);
      font-size: .82rem; line-height: 1.3;
    }
    .paweh-summary-row:last-child { border-bottom: 0; }
    .paweh-summary-label { color: var(--paweh-secondary); }
    .paweh-summary-value { color: var(--paweh-text); text-align: right; }
    .paweh-changed-notice, .paweh-caution {
      margin-top: 10px; padding: 9px 10px; border-left: 3px solid var(--paweh-caution);
      background: #F7F3EB; color: #665334; font-size: .84rem;
    }
    .paweh-error { padding: 10px; background: #F7EEEE; color: var(--paweh-error); }
    .paweh-interpretation { margin-top: 12px; padding: 14px; background: #F7F8FA; border: 1px solid var(--paweh-divider); }
    .paweh-interpretation h4 { font-size: 1.05rem; font-weight: 500; }
    .paweh-result-card .card-header { padding: 10px 14px; font-weight: 500; }
    .paweh-result-card .card-body { padding: 8px 14px 12px; }
    .paweh-result-card .table { margin-bottom: 0; font-size: .88rem; }
    .paweh-result-card .table > :not(caption) > * > * { padding: 7px 8px; }
    .paweh-model-specification {
      margin: 0 0 12px; padding: 9px 11px; color: var(--paweh-secondary);
      background: #F7F8FA; border: 1px solid var(--paweh-divider); font-size: .86rem;
    }
    .paweh-advanced-details, .paweh-advanced-visualization {
      margin-top: 12px; border: 1px solid var(--paweh-border); border-radius: 6px;
      background: #FFFFFF;
    }
    .paweh-advanced-details > summary, .paweh-advanced-visualization > summary {
      cursor: pointer; padding: 11px 13px; color: var(--paweh-text); font-weight: 500;
    }
    .paweh-advanced-details[open] > summary, .paweh-advanced-visualization[open] > summary {
      border-bottom: 1px solid var(--paweh-divider);
    }
    .paweh-advanced-body { padding: 10px 13px 13px; }
    .paweh-advanced-subtitle { margin: 0 0 10px; color: var(--paweh-secondary); font-size: .84rem; }
    .paweh-detail-section { margin-top: 12px; }
    .paweh-detail-section:first-of-type { margin-top: 0; }
    .paweh-detail-section h5 { margin: 0 0 4px; font-size: .92rem; font-weight: 500; }
    .paweh-reproduce { margin-top: 14px; }
    .paweh-reproduce h5 { margin-bottom: 5px; font-size: .92rem; font-weight: 500; }
    .paweh-reproduce pre {
      margin: 0; padding: 10px 12px; max-height: 280px; overflow: auto;
      color: #303840; background: #F7F8FA; border: 1px solid var(--paweh-divider);
      border-radius: 4px; font-size: .78rem; line-height: 1.45; user-select: text;
    }
    .paweh-sensitivity-controls {
      display: grid; grid-template-columns: minmax(240px, 320px) minmax(280px, 380px);
      gap: 0 18px; align-items: end; max-width: 720px; margin-bottom: 12px;
    }
    .paweh-sensitivity-controls .btn { width: auto; padding-left: 16px; padding-right: 16px; }
    .paweh-zoom-note { margin: 2px 0 8px; color: var(--paweh-muted); font-size: .82rem; }
    .paweh-visual-intro h3 { margin: 0 0 5px; font-size: 1.1rem; font-weight: 500; }
    .paweh-visual-intro p { margin-bottom: 8px; color: var(--paweh-secondary); }
    .text-muted { color: var(--paweh-muted) !important; }
    @media (max-width: 767.98px) {
      body.bslib-page-navbar > .container-fluid { padding-left: 12px; padding-right: 12px; }
      .paweh-page-heading { margin-top: 16px; }
      .paweh-workspace .bslib-sidebar-layout > .main { padding: 12px 0 0; }
      .paweh-workspace .bslib-sidebar-layout > .sidebar { width: auto; border-radius: 7px; }
      .paweh-sensitivity-controls { display: block; max-width: 100%; }
      .paweh-workspace .card-header-tabs { padding: 0 8px; overflow-x: auto; flex-wrap: nowrap; }
      .paweh-summary-row { grid-template-columns: minmax(0, 1fr); gap: 2px; }
      .paweh-summary-value { text-align: left; overflow-wrap: anywhere; }
      .paweh-reproduce pre { max-width: 100%; white-space: pre; overflow-x: auto; }
    }
    "
  )
}

.paweh_plot_colors <- function() {
  c(
    cases = "#3F4850", controls = "#6F879A",
    genotype = "#3F4850", trend = "#355C7D",
    tdt_baseline = "#3F4850", tdt_misclassification = "#355C7D",
    tdt_heterogeneity = "#6F879A", transmitted = "#3F4850",
    nontransmitted = "#8FA1AF",
    baseline = "#C7CDD2", adjusted = "#3F4850",
    reference = "#7A848C"
  )
}

.paweh_qtl_genotype_colors <- function() {
  c(`Genotype 0` = "#3F4850", `Genotype 1` = "#6F879A", `Genotype 2` = "#355C7D")
}

.paweh_plot_theme <- function() {
  ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 12.5, face = "plain"),
      axis.title = ggplot2::element_text(size = 10.5, color = "#3F4850"),
      axis.text = ggplot2::element_text(size = 9.5, color = "#5F6B76"),
      legend.title = ggplot2::element_text(size = 9.5),
      legend.text = ggplot2::element_text(size = 9.5),
      panel.grid.major = ggplot2::element_line(color = "#E9ECEF", linewidth = .35),
      panel.grid.minor = ggplot2::element_blank(),
      plot.background = ggplot2::element_rect(fill = "white", color = NA),
      panel.background = ggplot2::element_rect(fill = "white", color = NA),
      legend.background = ggplot2::element_rect(fill = "white", color = NA),
      plot.margin = ggplot2::margin(8, 8, 8, 8)
    )
}

.paweh_page_heading <- function(title, description) {
  shiny::div(
    class = "paweh-page-heading",
    shiny::h2(title),
    shiny::p(description)
  )
}

.paweh_summary_row <- function(label, value) {
  shiny::div(
    class = "paweh-summary-row",
    shiny::span(class = "paweh-summary-label", label),
    shiny::span(class = "paweh-summary-value", value)
  )
}

.paweh_detail_section <- function(title, rows) {
  shiny::div(
    class = "paweh-detail-section", shiny::h5(title),
    shiny::div(class = "paweh-summary-grid", rows)
  )
}

.paweh_repro_call <- function(function_name, args) {
  stopifnot(is.character(function_name), length(function_name) == 1L)
  call <- as.call(c(list(as.name(function_name)), unname(args)))
  names(call) <- c("", names(args))
  call
}

.paweh_call_text <- function(call) {
  paste(deparse(call, width.cutoff = 72L), collapse = "\n")
}

.paweh_format_percent <- function(x, digits = 1L) {
  vapply(x, function(value) {
    if (!is.finite(value)) return("not defined")
    paste0(formatC(100 * value, format = "f", digits = digits), "%")
  }, character(1), USE.NAMES = FALSE)
}

.paweh_format_count <- function(x, round_up = FALSE) {
  vapply(x, function(value) {
    if (is.na(value) || is.nan(value)) return("not defined")
    if (is.infinite(value)) return(if (value > 0) "Inf" else "-Inf")
    if (round_up) value <- ceiling(value)
    formatC(value, format = "d", big.mark = ",")
  }, character(1), USE.NAMES = FALSE)
}

.paweh_reproduce_ui <- function(call) {
  shiny::div(
    class = "paweh-reproduce", shiny::h5("Reproduce in R"),
    shiny::tags$pre(shiny::tags$code(.paweh_call_text(call)))
  )
}

.paweh_advanced_details_ui <- function(...) {
  shiny::tags$details(
    class = "paweh-advanced-details",
    shiny::tags$summary("Advanced calculation details"),
    shiny::div(
      class = "paweh-advanced-body",
      shiny::p(
        class = "paweh-advanced-subtitle",
        "Inspect model probabilities and intermediate quantities returned by the canonical calculation."
      ),
      ...
    )
  )
}

.paweh_power_axis_zoomed <- function(sensitivity) {
  if (is.null(sensitivity) || !identical(sensitivity$objective, "power")) return(FALSE)
  values <- sensitivity$data$y[is.finite(sensitivity$data$y)]
  length(values) > 1L && diff(range(values)) < 0.5
}

.paweh_placeholder_ui <- function(title, message) {
  shiny::div(
    class = "paweh-placeholder",
    shiny::h3(title),
    shiny::p(message)
  )
}

.paweh_empty_ui <- function(section) {
  messages <- c(
    Results = "Configure a design and select Calculate study design.",
    Sensitivity = "Calculate a design before exploring sensitivity.",
    Visualize = "Calculate a design to generate study-specific visualizations.",
    Methods = "Calculate a design to record its analysis specification."
  )
  stopifnot(section %in% names(messages))
  .paweh_placeholder_ui(section, unname(messages[[section]]))
}
