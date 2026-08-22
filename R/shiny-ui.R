# Top-level dashboard user interface.

.pawh_home_card_ui <- function(input_id, title, description, button_label) {
  bslib::card(
    class = "pawh-study-card",
    bslib::card_header(shiny::h3(class = "h5 mb-0", title)),
    bslib::card_body(
      shiny::p(description),
      shiny::actionButton(input_id, button_label, class = "btn-primary")
    )
  )
}

.pawh_home_ui <- function(id) {
  ns <- shiny::NS(id)
  shiny::tagList(
    shiny::div(
      class = "pawh-home-hero",
      shiny::h1("pawh"),
      shiny::p(
        class = "lead",
        "Power and Sample Size for Genetic Association Studies"
      ),
      shiny::p(
        class = "text-muted",
        "Choose a study type. Case-Control study design is available now; additional workflows are forthcoming."
      )
    ),
    shiny::div(
      class = "pawh-home-grid",
      bslib::layout_column_wrap(
        width = "280px",
        .pawh_home_card_ui(
          ns("case_control"),
          "Case-Control",
          "Compare genetic variation between cases and controls.",
          "Open Case-Control"
        ),
        .pawh_home_card_ui(
          ns("tdt"),
          "TDT / Family",
          "Design affected-child trio / transmission disequilibrium studies.",
          "Open TDT / Family"
        ),
        .pawh_home_card_ui(
          ns("qtl"),
          "Quantitative Trait",
          "Design studies involving continuous phenotypes.",
          "Open Quantitative Trait"
        )
      )
    )
  )
}

.pawh_home_server <- function(id, navigate) {
  shiny::moduleServer(id, function(input, output, session) {
    shiny::observeEvent(input$case_control, navigate("case_control"))
    shiny::observeEvent(input$tdt, navigate("tdt"))
    shiny::observeEvent(input$qtl, navigate("qtl"))
  })
}

.pawh_app_ui <- function() {
  bslib::page_navbar(
    title = "pawh",
    id = "main_nav",
    selected = "home",
    theme = .pawh_dashboard_theme(),
    window_title = "pawh - Power and Sample Size for Genetic Association Studies",
    bslib::nav_panel("Home", value = "home", .pawh_home_ui("home")),
    bslib::nav_panel(
      "Case-Control",
      value = "case_control",
      .pawh_case_control_ui("case_control")
    ),
    bslib::nav_panel(
      "TDT / Family",
      value = "tdt",
      .pawh_tdt_ui("tdt")
    ),
    bslib::nav_panel(
      "Quantitative Trait",
      value = "qtl",
      .pawh_qtl_ui("qtl")
    )
  )
}
