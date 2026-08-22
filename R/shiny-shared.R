# Shared presentation helpers for the pawh dashboard.

.pawh_dashboard_theme <- function() {
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
      --pawh-bg: #F7F8FA;
      --pawh-panel: #FFFFFF;
      --pawh-text: #20252B;
      --pawh-secondary: #5F6B76;
      --pawh-muted: #7B8792;
      --pawh-border: #D9DEE3;
      --pawh-divider: #E9ECEF;
      --pawh-primary: #355C7D;
      --pawh-caution: #A8844F;
      --pawh-error: #9A4D4D;
    }
    body { font-size: .925rem; color: var(--pawh-text); }
    .bslib-page-navbar > .navbar {
      min-height: 49px; padding-top: 0; padding-bottom: 0;
      background: #FFFFFF; border-bottom: 1px solid var(--pawh-border);
      box-shadow: none;
    }
    .navbar-brand { font-size: 1.2rem; font-weight: 500; }
    .navbar .nav-link, .navbar .navbar-nav > li > a {
      color: var(--pawh-secondary); padding: .82rem .75rem .68rem;
      border-bottom: 2px solid transparent; font-size: .925rem;
    }
    .navbar .nav-link.active, .navbar .navbar-nav > li.active > a {
      color: var(--pawh-text); background: transparent;
      border-bottom-color: var(--pawh-primary);
    }
    body.bslib-page-navbar > .container-fluid { padding: 0 20px 20px; }
    .pawh-page-heading { margin: 20px 0 22px; }
    .pawh-page-heading h2 {
      margin: 0; font-size: 1.75rem; line-height: 1.2; font-weight: 500;
    }
    .pawh-page-heading p {
      margin: 7px 0 0; color: var(--pawh-secondary); max-width: 700px;
      font-size: .975rem;
    }
    .pawh-home-hero { margin: 24px 0; padding: 0; }
    .pawh-home-hero h1 { margin: 0 0 5px; font-size: 1.8rem; font-weight: 500; }
    .pawh-home-hero .lead { margin-bottom: 6px; font-size: 1rem; color: var(--pawh-text); }
    .pawh-home-hero p:last-child { margin-bottom: 0; }
    .pawh-home-grid { margin-bottom: 24px; }
    .card {
      background: var(--pawh-panel); border-color: var(--pawh-border);
      border-radius: 7px; box-shadow: none;
    }
    .pawh-study-card { min-height: 0; height: auto; }
    .pawh-study-card .card-header {
      padding: 16px 18px 7px; border-bottom: 0; background: #FFFFFF;
    }
    .pawh-study-card .card-header h3 { font-size: 1.05rem; font-weight: 500; }
    .pawh-study-card .card-body {
      display: flex; flex-direction: column; gap: 12px; padding: 7px 18px 18px;
    }
    .pawh-study-card p { margin: 0; color: var(--pawh-secondary); }
    .pawh-study-card .btn { margin: 2px 0 0; align-self: flex-start; }
    .btn { border-radius: 5px; font-size: .9rem; }
    .btn-primary { background: var(--pawh-primary); border-color: var(--pawh-primary); }
    .pawh-workspace .bslib-sidebar-layout > .sidebar {
      width: 280px; background: #FFFFFF; border: 1px solid var(--pawh-border);
      border-radius: 7px 0 0 7px;
    }
    .pawh-workspace .bslib-sidebar-layout > .sidebar > .sidebar-content { padding: 16px; }
    .pawh-workspace .bslib-sidebar-layout > .sidebar .sidebar-title {
      margin-bottom: 12px; font-size: 1.05rem; font-weight: 500;
    }
    .pawh-workspace .form-group, .pawh-workspace .shiny-input-container {
      margin-bottom: 10px;
    }
    .pawh-workspace label, .pawh-workspace .control-label {
      margin-bottom: 4px; color: #3F4850; font-size: .91rem; font-weight: 500;
    }
    .pawh-workspace .form-control, .pawh-workspace .form-select {
      min-height: 36px; padding-top: 5px; padding-bottom: 5px; font-size: .9rem;
    }
    .pawh-workspace .radio, .pawh-workspace .checkbox { margin-bottom: 4px; }
    .pawh-sidebar-section {
      margin: 4px 0 12px; padding: 10px 0; border-top: 1px solid var(--pawh-divider);
    }
    .pawh-sidebar-section summary { cursor: pointer; font-size: .93rem; }
    .pawh-sidebar-section summary::after {
      content: 'Misclassification, heterogeneity, and related assumptions';
      display: block; margin-top: 3px; color: var(--pawh-muted);
      font-size: .8rem; font-weight: 400;
    }
    .pawh-sidebar-section[open] summary { margin-bottom: 10px; }
    .pawh-sidebar-section p, .pawh-placeholder p { color: var(--pawh-secondary); }
    .pawh-calculate { width: 100%; min-height: 40px; margin-top: 2px; }
    .pawh-workspace .bslib-sidebar-layout > .main {
      padding: 0 0 0 12px; background: transparent; border: 0;
    }
    .pawh-workspace, .pawh-workspace .bslib-sidebar-layout,
    .pawh-workspace .bslib-sidebar-layout > .main,
    .pawh-workspace .bslib-card {
      min-height: 0 !important; height: auto !important; flex: 0 0 auto;
    }
    .pawh-workspace .card-header-tabs { padding: 0 16px; background: #FFFFFF; }
    .pawh-workspace .nav-tabs { border-bottom-color: var(--pawh-border); }
    .pawh-workspace .nav-tabs .nav-link,
    .pawh-workspace .nav-tabs > li > a {
      padding: 11px 12px 9px; color: var(--pawh-secondary); font-size: .92rem;
      border: 0; border-bottom: 2px solid transparent; border-radius: 0;
    }
    .pawh-workspace .nav-tabs .nav-link.active,
    .pawh-workspace .nav-tabs > li.active > a {
      color: var(--pawh-text); border-bottom-color: var(--pawh-primary);
      background: #FFFFFF;
    }
    .pawh-workspace .card-body { padding: 16px; }
    .pawh-placeholder { max-width: 560px; padding: 16px 8px; }
    .pawh-placeholder h3 { margin-bottom: 6px; font-size: 1.08rem; font-weight: 500; }
    .pawh-placeholder p { margin-bottom: 0; }
    .pawh-design-summary { margin-top: 12px; }
    .pawh-design-summary .card-header { padding: 10px 12px; font-weight: 500; }
    .pawh-design-summary .card-body { padding: 8px 12px 10px; }
    .pawh-design-summary ul { margin: 0; padding: 0; list-style: none; }
    .pawh-design-summary li { padding: 5px 0; border-bottom: 1px solid var(--pawh-divider); font-size: .84rem; }
    .pawh-design-summary li:last-child { border-bottom: 0; }
    .pawh-summary-grid { display: grid; gap: 0; }
    .pawh-summary-row {
      display: grid; grid-template-columns: minmax(72px, .8fr) minmax(0, 1.2fr);
      gap: 8px; padding: 5px 0; border-bottom: 1px solid var(--pawh-divider);
      font-size: .82rem; line-height: 1.3;
    }
    .pawh-summary-row:last-child { border-bottom: 0; }
    .pawh-summary-label { color: var(--pawh-secondary); }
    .pawh-summary-value { color: var(--pawh-text); text-align: right; }
    .pawh-changed-notice, .pawh-caution {
      margin-top: 10px; padding: 9px 10px; border-left: 3px solid var(--pawh-caution);
      background: #F7F3EB; color: #665334; font-size: .84rem;
    }
    .pawh-error { padding: 10px; background: #F7EEEE; color: var(--pawh-error); }
    .pawh-interpretation { margin-top: 12px; padding: 14px; background: #F7F8FA; border: 1px solid var(--pawh-divider); }
    .pawh-interpretation h4 { font-size: 1.05rem; font-weight: 500; }
    .pawh-result-card .card-header { padding: 10px 14px; font-weight: 500; }
    .pawh-result-card .card-body { padding: 8px 14px 12px; }
    .pawh-result-card .table { margin-bottom: 0; font-size: .88rem; }
    .pawh-result-card .table > :not(caption) > * > * { padding: 7px 8px; }
    .pawh-model-specification {
      margin: 0 0 12px; padding: 9px 11px; color: var(--pawh-secondary);
      background: #F7F8FA; border: 1px solid var(--pawh-divider); font-size: .86rem;
    }
    .pawh-advanced-details, .pawh-advanced-visualization {
      margin-top: 12px; border: 1px solid var(--pawh-border); border-radius: 6px;
      background: #FFFFFF;
    }
    .pawh-advanced-details > summary, .pawh-advanced-visualization > summary {
      cursor: pointer; padding: 11px 13px; color: var(--pawh-text); font-weight: 500;
    }
    .pawh-advanced-details[open] > summary, .pawh-advanced-visualization[open] > summary {
      border-bottom: 1px solid var(--pawh-divider);
    }
    .pawh-advanced-body { padding: 10px 13px 13px; }
    .pawh-advanced-subtitle { margin: 0 0 10px; color: var(--pawh-secondary); font-size: .84rem; }
    .pawh-detail-section { margin-top: 12px; }
    .pawh-detail-section:first-of-type { margin-top: 0; }
    .pawh-detail-section h5 { margin: 0 0 4px; font-size: .92rem; font-weight: 500; }
    .pawh-reproduce { margin-top: 14px; }
    .pawh-reproduce h5 { margin-bottom: 5px; font-size: .92rem; font-weight: 500; }
    .pawh-reproduce pre {
      margin: 0; padding: 10px 12px; max-height: 280px; overflow: auto;
      color: #303840; background: #F7F8FA; border: 1px solid var(--pawh-divider);
      border-radius: 4px; font-size: .78rem; line-height: 1.45; user-select: text;
    }
    .pawh-sensitivity-controls {
      display: grid; grid-template-columns: minmax(240px, 320px) minmax(280px, 380px);
      gap: 0 18px; align-items: end; max-width: 720px; margin-bottom: 12px;
    }
    .pawh-sensitivity-controls .btn { width: auto; padding-left: 16px; padding-right: 16px; }
    .pawh-zoom-note { margin: 2px 0 8px; color: var(--pawh-muted); font-size: .82rem; }
    .pawh-visual-intro h3 { margin: 0 0 5px; font-size: 1.1rem; font-weight: 500; }
    .pawh-visual-intro p { margin-bottom: 8px; color: var(--pawh-secondary); }
    .text-muted { color: var(--pawh-muted) !important; }
    @media (max-width: 767.98px) {
      body.bslib-page-navbar > .container-fluid { padding-left: 12px; padding-right: 12px; }
      .pawh-page-heading { margin-top: 16px; }
      .pawh-workspace .bslib-sidebar-layout > .main { padding: 12px 0 0; }
      .pawh-workspace .bslib-sidebar-layout > .sidebar { width: auto; border-radius: 7px; }
      .pawh-sensitivity-controls { display: block; max-width: 100%; }
    }
    "
  )
}

.pawh_plot_colors <- function() {
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

.pawh_plot_theme <- function() {
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

.pawh_page_heading <- function(title, description) {
  shiny::div(
    class = "pawh-page-heading",
    shiny::h2(title),
    shiny::p(description)
  )
}

.pawh_summary_row <- function(label, value) {
  shiny::div(
    class = "pawh-summary-row",
    shiny::span(class = "pawh-summary-label", label),
    shiny::span(class = "pawh-summary-value", value)
  )
}

.pawh_detail_section <- function(title, rows) {
  shiny::div(
    class = "pawh-detail-section", shiny::h5(title),
    shiny::div(class = "pawh-summary-grid", rows)
  )
}

.pawh_repro_call <- function(function_name, args) {
  stopifnot(is.character(function_name), length(function_name) == 1L)
  call <- as.call(c(list(as.name(function_name)), unname(args)))
  names(call) <- c("", names(args))
  call
}

.pawh_call_text <- function(call) {
  paste(deparse(call, width.cutoff = 72L), collapse = "\n")
}

.pawh_reproduce_ui <- function(call) {
  shiny::div(
    class = "pawh-reproduce", shiny::h5("Reproduce in R"),
    shiny::tags$pre(shiny::tags$code(.pawh_call_text(call)))
  )
}

.pawh_advanced_details_ui <- function(...) {
  shiny::tags$details(
    class = "pawh-advanced-details",
    shiny::tags$summary("Advanced calculation details"),
    shiny::div(
      class = "pawh-advanced-body",
      shiny::p(
        class = "pawh-advanced-subtitle",
        "Inspect model probabilities and intermediate quantities returned by the canonical calculation."
      ),
      ...
    )
  )
}

.pawh_power_axis_zoomed <- function(sensitivity) {
  if (is.null(sensitivity) || !identical(sensitivity$objective, "power")) return(FALSE)
  values <- sensitivity$data$y[is.finite(sensitivity$data$y)]
  length(values) > 1L && diff(range(values)) < 0.5
}

.pawh_placeholder_ui <- function(title, message) {
  shiny::div(
    class = "pawh-placeholder",
    shiny::h3(title),
    shiny::p(message)
  )
}

.pawh_design_sidebar_ui <- function(objective, phase) {
  bslib::sidebar(
    title = "Design setup",
    open = "desktop",
    shiny::div(
      class = "pawh-sidebar-section",
      shiny::h4("Objective"),
      shiny::p(objective)
    ),
    shiny::div(
      class = "pawh-sidebar-section",
      shiny::h4("Core assumptions"),
      shiny::p("Study inputs will be added in a focused implementation phase.")
    ),
    shiny::tags$details(
      class = "pawh-sidebar-section",
      shiny::tags$summary(shiny::strong("Advanced assumptions")),
      shiny::p("Advanced error, heterogeneity, and modeling options will appear here.")
    ),
    shiny::tags$button(
      type = "button",
      class = "btn btn-primary w-100 mt-3",
      disabled = NA,
      "Calculate"
    ),
    shiny::tags$small(class = "text-muted", phase),
    bslib::card(
      class = "pawh-design-summary",
      bslib::card_header("Your design"),
      bslib::card_body(
        shiny::p("Your selected assumptions will be summarized here.")
      )
    )
  )
}

.pawh_study_workspace_ui <- function(id, objective, phase) {
  ns <- shiny::NS(id)
  bslib::layout_sidebar(
    sidebar = .pawh_design_sidebar_ui(objective, phase),
    bslib::navset_card_tab(
      id = ns("section"),
      bslib::nav_panel(
        "Results",
        .pawh_placeholder_ui(
          "Results",
          paste(phase, "No calculations are performed by this skeleton.")
        )
      ),
      bslib::nav_panel(
        "Sensitivity",
        .pawh_placeholder_ui(
          "Sensitivity",
          "Sensitivity controls and parameter sweeps will be implemented later."
        )
      ),
      bslib::nav_panel(
        "Visualize",
        .pawh_placeholder_ui(
          "Visualize",
          "Study-design visualizations will be connected to canonical plotting functions later."
        )
      ),
      bslib::nav_panel(
        "Methods",
        .pawh_placeholder_ui(
          "Methods",
          "Model assumptions, supported methods, and interpretation guidance will appear here."
        )
      )
    )
  )
}
