# Quantitative-trait dashboard shell.

.pawh_qtl_subtype_ui <- function(id, objective, phase) {
  .pawh_study_workspace_ui(id, objective = objective, phase = phase)
}

.pawh_qtl_ui <- function(id) {
  ns <- shiny::NS(id)
  shiny::tagList(
    .pawh_page_heading(
      "Quantitative Trait",
      "Design studies involving continuous phenotypes."
    ),
    bslib::navset_pill_list(
      id = ns("subtype"),
      widths = c(3, 9),
      well = TRUE,
      bslib::nav_panel(
        "Full continuous trait",
        .pawh_qtl_subtype_ui(
          ns("continuous"),
          objective = "Design a study that measures a continuous phenotype in all participants.",
          phase = "Full continuous-trait inputs will be implemented in a later development phase."
        )
      ),
      bslib::nav_panel(
        "Extreme phenotype sampling",
        .pawh_qtl_subtype_ui(
          ns("extreme"),
          objective = "Design a study that recruits participants from selected phenotype tails.",
          phase = "Extreme-phenotype inputs will be implemented in a later development phase."
        )
      ),
      bslib::nav_panel(
        "Multiple quantitative traits",
        .pawh_qtl_subtype_ui(
          ns("multivariate"),
          objective = "Design a study that analyzes multiple quantitative phenotypes jointly.",
          phase = "Multiple-trait inputs will be implemented in a later development phase."
        )
      )
    )
  )
}
