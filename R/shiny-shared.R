# Shared presentation helpers for the pawh dashboard.

.pawh_dashboard_theme <- function() {
  theme <- bslib::bs_theme(
    version = 5,
    bg = "#f5f7f8",
    fg = "#24313a",
    primary = "#245b78",
    secondary = "#647681"
  )
  bslib::bs_add_rules(
    theme,
    "
    .pawh-page-heading { margin-bottom: 1rem; }
    .pawh-page-heading p { color: #52636d; max-width: 70ch; }
    .pawh-study-card { min-height: 15rem; }
    .pawh-study-card .card-body {
      display: flex;
      flex-direction: column;
      gap: 0.75rem;
    }
    .pawh-study-card .btn { margin-top: auto; align-self: flex-start; }
    .pawh-sidebar-section + .pawh-sidebar-section { margin-top: 1.25rem; }
    .pawh-sidebar-section p, .pawh-placeholder p { color: #52636d; }
    .pawh-design-summary { margin-top: 1.25rem; }
    "
  )
}

.pawh_page_heading <- function(title, description) {
  shiny::div(
    class = "pawh-page-heading",
    shiny::h2(title),
    shiny::p(description)
  )
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
