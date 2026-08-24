# Top-level dashboard user interface.

.paweh_home_card_ui <- function(input_id, title, description, button_label) {
  bslib::card(
    class = "paweh-study-card",
    bslib::card_header(shiny::h3(class = "h5 mb-0", title)),
    bslib::card_body(
      shiny::p(description),
      shiny::actionButton(input_id, button_label, class = "btn-primary")
    )
  )
}

.paweh_home_ui <- function(id) {
  ns <- shiny::NS(id)
  shiny::tagList(
    shiny::div(
      class = "paweh-home-hero",
      shiny::h1("paweh"),
      shiny::p(
        class = "lead",
        "Power and Sample Size for Genetic Association Studies"
      ),
      shiny::p(
        class = "text-muted",
        "Choose a study type. Case-Control, TDT / Family, and Quantitative Trait study-design workflows are available now."
      )
    ),
    shiny::div(
      class = "paweh-home-grid",
      bslib::layout_column_wrap(
        width = "280px",
        .paweh_home_card_ui(
          ns("case_control"),
          "Case-Control",
          "Compare genetic variation between cases and controls.",
          "Open Case-Control"
        ),
        .paweh_home_card_ui(
          ns("tdt"),
          "TDT / Family",
          "Design affected-child trio / transmission disequilibrium studies.",
          "Open TDT / Family"
        ),
        .paweh_home_card_ui(
          ns("qtl"),
          "Quantitative Trait",
          "Design continuous, extreme-phenotype, or joint multiple-trait studies.",
          "Open Quantitative Trait"
        )
      )
    )
  )
}

.paweh_home_server <- function(id, navigate) {
  shiny::moduleServer(id, function(input, output, session) {
    shiny::observeEvent(input$case_control, navigate("case_control"))
    shiny::observeEvent(input$tdt, navigate("tdt"))
    shiny::observeEvent(input$qtl, navigate("qtl"))
  })
}

.paweh_app_ui <- function() {
  bslib::page_navbar(
    title = "paweh",
    id = "main_nav",
    selected = "home",
    theme = .paweh_dashboard_theme(),
    window_title = "paweh - Power and Sample Size for Genetic Association Studies",
    bslib::nav_panel("Home", value = "home", .paweh_home_ui("home")),
    bslib::nav_panel(
      "Case-Control",
      value = "case_control",
      .paweh_case_control_ui("case_control")
    ),
    bslib::nav_panel(
      "TDT / Family",
      value = "tdt",
      .paweh_tdt_ui("tdt")
    ),
    bslib::nav_panel(
      "Quantitative Trait",
      value = "qtl",
      .paweh_qtl_ui("qtl")
    )
  )
}
