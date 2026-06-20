# tests/testthat/test-optimize-methods.R
# Regression tests for the optimize_CS() method engine, especially the Kneedle
# elbow detector that fixes the plateau-then-cliff asymmetry where the old
# 'dynamic' method stopped *before* the main drop.

test_that(".ocs_knee_cs finds the elbow after the cliff, not the plateau", {
    cs <- c(0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100)
    pct <- c(100, 97, 67, 28, 9, 7, 6, 6, 5, 5, 5)
    cdf <- data.frame(CS = cs, Pct_Retained = pct)
    cdf$Step_Loss_Pct <- c(0, head(pct, -1) - tail(pct, -1))

    knee <- karioCaS:::.ocs_knee_cs(cdf)
    # The elbow must sit at/after the cliff (>= CS20), never on the CS10 plateau.
    expect_true(knee >= 20)
    expect_true(knee <= 50)
})

test_that(".ocs_postcliff_cs lands at or after the steepest drop", {
    cs <- c(0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100)
    pct <- c(100, 97, 67, 28, 9, 7, 6, 6, 5, 5, 5)
    cdf <- data.frame(CS = cs, Pct_Retained = pct)
    cdf$Step_Loss_Pct <- c(0, head(pct, -1) - tail(pct, -1))
    toll <- karioCaS:::.ocs_dyn_toll(cdf)
    pc <- karioCaS:::.ocs_postcliff_cs(cdf, toll, max_cs = 100)
    # Steepest single drop is into CS30; stabilisation must be >= that.
    expect_true(pc >= 30)
})

test_that("optimize_CS defaults to kneedle and accepts the new methods", {
    expect_equal(eval(formals(optimize_CS)$method)[1], "kneedle")
    expect_true(all(
        c("kneedle", "postcliff", "segmented", "dynamic", "manual") %in%
            eval(formals(optimize_CS)$method)
    ))
})
