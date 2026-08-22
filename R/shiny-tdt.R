# TDT and family-based dashboard shell.

.pawh_tdt_ui <- function(id) {
  ns <- shiny::NS(id)
  shiny::tagList(
    .pawh_page_heading(
      "TDT / Family",
      "Design affected-child trio / transmission disequilibrium studies."
    ),
    .pawh_study_workspace_ui(
      ns("workspace"),
      objective = "Estimate power or required affected-child trios for a family-based study.",
      phase = "TDT / Family study design will be implemented in a later development phase."
    )
  )
}
