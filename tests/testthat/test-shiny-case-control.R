test_that("Case-Control UI exposes the production workflow", {
  html <- paste(as.character(pawh:::.pawh_case_control_ui("cc")), collapse = "\n")
  expect_match(html, "Estimate power")
  expect_match(html, "Minimum sample size")
  expect_match(html, "Direct genotype probabilities")
  expect_match(html, "Advanced assumptions")
  expect_match(html, "Sensitivity")
  expect_match(html, "Visualize")
  expect_match(html, "Methods")
})

test_that("model and direct snapshots reproduce canonical power and MSSN", {
  v <- pawh:::.pawh_cc_defaults()
  calc <- pawh:::.pawh_cc_calculate(pawh:::.pawh_cc_snapshot(v))
  direct <- cc_power(N_case = 500, alpha = .05, input_mode = "model_based",
    prev = .1, pd = .3, R2 = 2, MOI = "M", k = 1, verbose = FALSE)
  expect_equal(calc$adjusted$tests, direct$tests)

  v$objective <- "mssn"; v$input_mode <- "model_free"
  calc <- pawh:::.pawh_cc_calculate(pawh:::.pawh_cc_snapshot(v))
  direct <- cc_mssn(power = .8, alpha = .05, input_mode = "model_free",
    g1 = c(.35,.45,.2), g0 = c(.49,.42,.09), k = 1, verbose = FALSE)
  expect_equal(calc$adjusted$tests, direct$tests)
})

test_that("advanced modifiers and zero-error narration follow canonical metadata", {
  v <- pawh:::.pawh_cc_defaults()
  v$locus_het <- TRUE; v$pi <- 75; v$pheno_misclass <- TRUE
  v$theta <- 3; v$phi <- 1; v$genotype_error <- TRUE; v$e <- 2
  calc <- pawh:::.pawh_cc_calculate(pawh:::.pawh_cc_snapshot(v))
  direct <- cc_power(N_case=500,alpha=.05,input_mode="model_based",prev=.1,
    pd=.3,R2=2,MOI="M",k=1,locus_het=TRUE,pi=.75,
    pheno_misclass=TRUE,theta=.03,phi=.01,geno_misclass="1p",e=.02,
    verbose=FALSE)
  expect_equal(calc$adjusted$tests, direct$tests)
  expect_true(all(unlist(calc$active)))

  for (model in c("1p","2p","3p","diff3p")) {
    z <- pawh:::.pawh_cc_defaults(); z$genotype_error <- TRUE; z$geno_misclass <- model
    zero <- pawh:::.pawh_cc_calculate(pawh:::.pawh_cc_snapshot(z))
    expect_false(zero$active$genotype)
    expect_equal(zero$adjusted$tests, zero$baseline$tests)
  }
})

test_that("validation is friendly and calculated state remains frozen", {
  v <- pawh:::.pawh_cc_defaults(); v$input_mode <- "model_free"; v$g1_2 <- .3
  expect_error(pawh:::.pawh_cc_snapshot(v), "sum to 1")
  shiny::testServer(pawh:::.pawh_case_control_server, {
    session$setInputs(objective="power",input_mode="model_based",N_case=500,
      alpha=.05,k=1,prev=10,pd=30,MOI="M",R2=2)
    session$setInputs(calculate=1)
    first <- session$returned$calculation()$adjusted$tests$genotypes$power
    session$setInputs(alpha=.01)
    session$flushReact()
    expect_true(session$returned$changed())
    expect_equal(session$returned$calculation()$adjusted$tests$genotypes$power, first)
  })
})

test_that("sensitivity points are direct canonical calls of frozen input", {
  v <- pawh:::.pawh_cc_defaults()
  calc <- pawh:::.pawh_cc_calculate(pawh:::.pawh_cc_snapshot(v))
  sens <- pawh:::.pawh_cc_sensitivity(calc, "N_case", c(400,600), n=5)
  expect_equal(nrow(sens$data), 10)
  expected <- cc_power(N_case=500,alpha=.05,input_mode="model_based",prev=.1,
    pd=.3,R2=2,MOI="M",k=1,verbose=FALSE)
  middle <- subset(sens$data, x == 500)
  expect_equal(middle$y, c(expected$tests$genotypes$power, expected$tests$trend$power))
  expect_equal(sens$baseline_x, 500)
})
