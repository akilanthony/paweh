# Case-control dashboard shell.

.pawh_case_control_ui <- function(id) {
  ns <- shiny::NS(id)
  shiny::tagList(
    .pawh_page_heading(
      "Case-Control",
      "Compare genetic variation between cases and controls."
    ),
    .pawh_study_workspace_ui(
      ns("workspace"),
      objective = "Estimate power or required sample size for a case-control association study.",
      phase = "Case-Control study design will be implemented in P3B."
    )
  )
}
