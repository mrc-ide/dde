context("events")

test_that("change variables", {
  target_r <- function(t, y, p) {
    -y * p
  }
  event_r <- function(t, y, p) {
    y * 2
  }
  target_c <- "exponential"
  event_c <- "double_variables"
  dllname <- "dde_growth"

  set.seed(1)
  y0 <- runif(5)
  p <- rep(1, length(y0))
  tol <- 1e-8
  tt <- seq(0, 5, length.out = 51) # 5001

  cmp1_r <- dopri(y0, tt, target_r, p, atol = tol, rtol = tol)
  cmp2_r <- dopri(y0, tt, target_r, p, tcrit = 1:5, atol = tol, rtol = tol)
  cmp1_c <- dopri(y0, tt, target_c, p, dllname = dllname,
                  atol = tol, rtol = tol)
  cmp2_c <- dopri(y0, tt, target_c, p, dllname = dllname, tcrit = 1:5,
                  atol = tol, rtol = tol)

  expect_equal(cmp1_r, cmp2_r, tolerance = 1e-6)
  expect_equal(cmp1_r, cmp1_c, tolerance = 1e-12)
  expect_equal(cmp2_r, cmp2_c, tolerance = 1e-12)

  te <- 1:5

  ## First, the dummy event:
  res1_r <- dopri(y0, tt, target_r, p,
                  event_time = te, event_function = function(t, y, p) y,
                  atol = tol, rtol = tol)
  expect_equal(res1_r, cmp2_r, tolerance = 1e-12)

  ## Then the real event where we double the output variable:
  res2_r <- dopri(y0, tt, target_r, p,
                  event_time = te, event_function = event_r,
                  atol = tol, rtol = tol)

  ## This does seem to throw the tolerances out of whack quite badly,
  ## but then they're being computed on totally different quantities.
  m <- 2^findInterval(tt, te + 1e-8,
                      rightmost.closed = TRUE)
  expect_equal(res2_r[, -1], cmp2_r[, -1] * m, tolerance = 1e-5)

  ## Confirm that the C cases agree fairly well:
  res1_c <- dopri(y0, tt, target_c, p, dllname = dllname,
                  event_time = te, event_function = "identity",
                  atol = tol, rtol = tol)
  expect_equal(res1_c, res1_r, tolerance = 1e-12)
  res2_c <- dopri(y0, tt, target_c, p, dllname = dllname,
                  event_time = te, event_function = event_c,
                  atol = tol, rtol = tol)
  expect_equal(res2_c, res2_r, tolerance = 1e-12)
})

test_that("change parameters", {
  target_r <- function(t, y, p) {
    p
  }
  event_r <- function(t, y, p) {
    attr(y, "parms") <- p * 2
    y
  }

  target_c <- "linear"
  event_c <- "double_parameters"
  dllname <- "dde_growth"

  set.seed(1)
  y0 <- runif(5)
  p <- runif(5)
  p0 <- p
  tol <- 1e-8
  tt <- seq(0, 5, length.out = 51) # 5001

  cmp1_r <- dopri(y0, tt, target_r, p, atol = tol, rtol = tol)
  cmp2_r <- dopri(y0, tt, target_r, p, tcrit = 1:5, atol = tol, rtol = tol)
  cmp1_c <- dopri(y0, tt, target_c, p, dllname = dllname,
                  atol = tol, rtol = tol)
  cmp2_c <- dopri(y0, tt, target_c, p, dllname = dllname, tcrit = 1:5,
                  atol = tol, rtol = tol)

  expect_equal(cmp1_r, cmp2_r, tolerance = 1e-6)
  expect_equal(cmp1_r, cmp1_c, tolerance = 1e-12)
  expect_equal(cmp2_r, cmp2_c, tolerance = 1e-12)

  ## TODO: This should surely imply that the final event should fire
  ## too?  But it's not firing here.  I think that this way around is
  ## correct though (but we *should* fire on starting events).  It's
  ## as if events fire at t + eps, with eps lim -> 0
  te <- 1:5

  ## First, the dummy event (duplicates the test above really)
  res1_r <- dopri(y0, tt, target_r, p,
                  event_time = te, event_function = function(t, y, p) y,
                  atol = tol, rtol = tol)
  expect_equal(res1_r, cmp2_r, tolerance = 1e-12)

  ## Then the real event where we double the parameters
  res2_r <- dopri(y0, tt, target_r, p,
                  event_time = te, event_function = event_r,
                  atol = tol, rtol = tol)
  expect_identical(p, p0) # unchanged -- this is important!

  cmp <- matrix(NA, length(tt), length(y0))
  cmp[1, ] <- y0
  ## Then start the process of computing the trajectory.  This is
  ## actually quite a faff; surprisingly more than the above because I
  ## don't know how to do this sort of thing very effectively.
  times <- c(tt[1], te)
  for (i in seq_along(te)) {
    j <- which(tt >= times[i] & tt <= times[i + 1])
    cmp[j, ] <- t(cmp[j[1], ] + outer(p * 2 ^ (i - 1), tt[j] - tt[j][1]))
  }

  expect_equal(res2_r[, -1], cmp, tolerance = 1e-5)

  ## Confirm that the C cases agree fairly well:
  res1_c <- dopri(y0, tt, target_c, p, dllname = dllname,
                  event_time = te, event_function = "identity",
                  atol = tol, rtol = tol)
  expect_equal(res1_c, res1_r, tolerance = 1e-12)
  res2_c <- dopri(y0, tt, target_c, p, dllname = dllname,
                  event_time = te, event_function = event_c,
                  atol = tol, rtol = tol)
  expect_equal(res2_c, res2_r, tolerance = 1e-12)
})

test_that("vector of events", {
  target_r <- function(t, y, p) {
    -y * p
  }
  event_r1 <- function(t, y, p) {
    y * 2
  }
  event_r2 <- function(t, y, p) {
    y / 2
  }
  target_c <- "exponential"
  event_c1 <- "double_variables"
  event_c2 <- "halve_variables"
  dllname <- "dde_growth"

  set.seed(1)
  y0 <- runif(5)
  p <- rep(1, length(y0))
  tol <- 1e-8
  tt <- seq(0, 3, length.out = 31) # 5001

  te <- c(1, 2)
  events_r <- list(event_r1, event_r2)
  events_c <- list(event_c1, event_c2)

  res_r <- dopri(y0, tt, target_r, p,
                 event_time = te, event_function = events_r,
                 atol = tol, rtol = tol)
  res_c <- dopri(y0, tt, target_c, p, dllname = dllname,
                 event_time = te, event_function = events_c,
                 atol = tol, rtol = tol)

  cmp <- t(y0 * exp(outer(-p, tt)))
  m <- 1 + findInterval(tt, te + 1e-8, rightmost.closed = TRUE) %% 2

  expect_equal(res_r[, -1], cmp * m, tolerance = 1e-5)
  expect_equal(res_r, res_c, tolerance = 1e-12)
})

test_that("interleave events and critical", {
  ## We'll use the same trick as earlier to detect if tcrit was used;
  ## check the number of steps taken as it should be smaller.  This
  ## means that our target is going to be a step function.
  target <- function(t, y, p) {
    if (t > 2) -p else p
  }
  event1 <- function(t, y, p) {
    y * 2
  }
  event2 <- function(t, y, p) {
    y / 2
  }

  te <- c(1, 3)
  tcrit <- 2
  events <- list(event1, event2)

  set.seed(1)
  y0 <- runif(5)
  p <- runif(length(y0))
  tol <- 1e-8
  tt <- seq(0, 4, length.out = 41) # 5001

  res1 <- dopri(y0, tt, target, p,
                event_time = te, event_function = events,
                atol = tol, rtol = tol, return_statistics = TRUE)

  res2 <- dopri(y0, tt, target, p,
                event_time = te, event_function = events, tcrit = 2,
                atol = tol, rtol = tol, return_statistics = TRUE)

  s1 <- attr(res1, "statistics")
  s2 <- attr(res2, "statistics")
  expect_lt(s2[["n_step"]], s1[["n_step"]])
  expect_lt(s2[["n_reject"]], s1[["n_reject"]])

  cmp <- outer(p, tt) + y0
  ## The first double
  j <- which(tt >= 1 & tt <= 2)
  cmp[, j[-1]] <- cmp[, j[-1]] + cmp[, j[[1]]]
  ## Flips around:
  j <- which(tt >= 2)
  cmp[, j[-1]] <- outer(-p, tt[j[-1]] - tt[j[[1]]]) + cmp[, j[[1]]]
  ## The halve
  j <- which(tt >= 3)
  cmp[, j[-1]] <- cmp[, j[-1]] - cmp[, j[[1]]] / 2

  expect_equal(res2[, -1], t(cmp), tolerance = 1e-5)
  expect_equal(res1[, -1], res2[, -1], tolerance = 1e-5)
})

test_that("event ordering when stacked", {
  target <- function(t, y, p) {
    0
  }
  called <- numeric(0)
  event1 <- function(t, y, p) {
    called <<- c(called, c(1, t))
    y
  }
  event2 <- function(t, y, p) {
    called <<- c(called, c(2, t))
    y
  }

  te <- c(1, 1)
  events <- list(event1, event2)

  set.seed(1)
  y0 <- runif(5)
  p <- runif(length(y0))
  tol <- 1e-8
  tt <- seq(0, 2, length.out = 21) # 5001

  res1 <- dopri(y0, tt, target, p,
                event_time = te, event_function = events,
                atol = tol, rtol = tol, return_statistics = TRUE)
  m <- matrix(called, 2)
  expect_equal(m[1, ], c(1, 2))
  expect_equal(m[2, ], te)

  ## And again, with more events:
  called <- numeric(0)
  te <- rep(c(0.8, 1.1), each = 4)
  i <- sample(2, length(te), replace = TRUE)
  res2 <- dopri(y0, tt, target, p,
                event_time = te, event_function = events[i],
                atol = tol, rtol = tol, return_statistics = TRUE)
  m <- matrix(called, 2)
  expect_equal(m[1, ], i)
  expect_equal(m[2, ], te)
})

test_that("events colliding with tcrit", {
  target <- function(t, y, p) {
    0
  }
  called <- numeric(0)
  event1 <- function(t, y, p) {
    called <<- c(called, c(1, t))
    y
  }
  event2 <- function(t, y, p) {
    called <<- c(called, c(2, t))
    y
  }
  events <- rep(list(event1, event2), 5)

  te <- rep(1:5, each = 2)
  tcrit <- c(2, 3.5, 4)
  ans <- check_events(te, events, tcrit)
  expect_equal(ans$tcrit, sort(c(te, setdiff(tcrit, te))))

  set.seed(1)
  y0 <- 0
  p <- NULL
  tt <- seq(0, 6)

  res <- dopri(y0, tt, target, p,
               tcrit = tcrit, event_time = te, event_function = events)
  m <- matrix(called, 2)
  expect_equal(m[1, ], rep(1:2, 5))
  expect_equal(m[2, ], te)
})

test_that("single event", {
  target <- function(t, y, p) {
    0
  }
  event <- function(t, y, p) {
    message("this is an event")
    y
  }
  expect_message(dopri(0, 0:4, target, NULL,
                       event_time = 1, event_function = event),
                 "this is an event")
})

test_that("no crash after serialisation", {
  dllname <- "dde_growth"
  target <- "exponential"
  ## TODO: passing in $address here does not work:
  event <- getNativeSymbolInfo("double_variables", dllname)

  y0 <- 1
  p <- 0
  tt <- seq(0, 2, length.out = 21)
  res1 <- dopri(y0, tt, target, p, dllname = dllname,
                event_time = 1, event_function = event)
  expect_equal(res1[, 2], ifelse(tt > 1.0, 2, 1))

  event0 <- make_null_pointer(event)
  expect_error(dopri(y0, tt, target, p, dllname = dllname,
                     event_time = 1, event_function = event0),
               "Was passed null pointer for events[1]", fixed = TRUE)
})

test_that("error cases", {
  target <- function(t, y, p) {
    0
  }
  event <- function(t, y, p) {
    message("this is an event")
    y
  }

  expect_error(dopri(0, 0:4, target, NULL,
                     event_time = NULL, event_function = event),
               "'event_function' given without 'event_time'")
  expect_error(dopri(0, 0:4, target, NULL,
                     event_time = 1, event_function = NULL),
               "'event_time' given without 'event_function'")
  expect_error(dopri(0, 0:4, target, NULL,
                     event_time = 1, event_function = "double_variables"),
               "'event_function' must be an R function")
  expect_error(dopri(0, 0:4, "exponential", NULL, dllname = "dde_growth",
                     event_time = 1, event_function = event),
               "'event_function' must be a compiled function")
  expect_error(dopri(0, 0:4, target, NULL,
                     event_time = 1:3, event_function = list(event, event)),
               "'event_function' must be a single event or a list of length 3")
})
