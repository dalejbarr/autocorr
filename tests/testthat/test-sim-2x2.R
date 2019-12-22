test_that("randomization is correct", {
  dat <- sim_2x2(24, 24, 0, 0, 0, 0, 1, 1, .5, 1, TRUE)
  err_r <- dat[["Y_r"]] - dat[["Y_fit"]]
  err_b <- dat[["Y_b"]] - dat[["Y_fit"]]
  expect_equal(
    err_r[order(dat[["subj_id"]], dat[["tnum_r"]])],
    err_b[order(dat[["subj_id"]], dat[["tnum_b"]])])
})

test_that("intercept is correct", {
  dat <- sim_2x2(24, 24, 5, 0, 0, 0, 0, 0, 0, 1, TRUE)
  expect_equal(dat[["Y_fit"]],
               rep(5, length(dat[["Y_fit"]])))
})

test_that("A is correct", {
  dat <- sim_2x2(24, 24, 0, 5, 0, 0, 0, 0, 0, 1, TRUE)
  agg <- aggregate(Y_fit ~ A, dat, mean)
  expect_equal(agg[["Y_fit"]][2] - agg[["Y_fit"]][1], 5)
})

test_that("B is correct", {
  dat <- sim_2x2(24, 24, 0, 0, 5, 0, 0, 0, 0, 1, TRUE)
  agg <- aggregate(Y_fit ~ B, dat, mean)
  expect_equal(agg[["Y_fit"]][2] - agg[["Y_fit"]][1], 5)
})

test_that("AB is correct", {
  dat <- sim_2x2(24, 24, 0, 0, 0, 5, 0, 0, 0, 1, TRUE)
  agg <- aggregate(Y_fit ~ A + B, dat, mean)[["Y_fit"]]
  expect_equal(agg[1] - agg[2] - (agg[3] - agg[4]), 5)
})

test_that("all residual SDs are 1", {
  ff <- lapply(1:9, function(.i) {
    dat <- sim_2x2(12, 12, 0, 0, 0, 0, 1, 1, .5, .i, TRUE)
    err <- dat[["Y_r"]] - dat[["Y_fit"]]
    ds <- split(err, dat[["subj_id"]])
    vv <- sapply(ds, sd)
    names(vv) <- NULL
    vv
  })
  expect_true(all(sapply(ff, function(.x) {
    all.equal(.x, rep(1, 12))
  })))
})

test_that("all residual means are 1", {
  ff <- lapply(1:9, function(.i) {
    dat <- sim_2x2(12, 12, 0, 0, 0, 0, 1, 1, .5, .i, TRUE)
    err <- dat[["Y_r"]] - dat[["Y_fit"]]
    ds <- split(err, dat[["subj_id"]])
    vv <- sapply(ds, mean)
    names(vv) <- NULL
    vv
  })
  expect_true(all(sapply(ff, function(.x) {
    all.equal(.x, rep(0, 12), tolerance = 1e-10)
  })))
})
