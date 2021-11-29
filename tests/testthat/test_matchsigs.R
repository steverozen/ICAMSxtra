test_that("MatchSigs2Diretions", {
  seta <- matrix(c(1, 3,   4, 1, 2, 4), ncol = 2)
  setb <- matrix(c(1, 3.1, 4, 5, 1, 1, 1, 2.8, 4), ncol = 3)
  colnames(seta) <- c("A.1", "A.2")
  colnames(setb) <- c("B.1", "B.2", "B.3")
  expected <-
    list(
      averCosSim = 0.890140132539195,
      match1 = structure(list(
        to = c("B.1", "B.3"),
        sim = c(0.999877135226674, 0.989516398944649)),
        class = "data.frame", row.names = c("A.1", "A.2")),
      match2 = structure(list(
        to = c("A.1", "A.2", "A.1"),
        sim = c(0.999877135226674, 0.461956578122389, 0.999473415175591)),
        class = "data.frame", row.names = c("B.1",  "B.2", "B.3")))

  tout <- MatchSigs2Directions(seta, setb)
  expect_equal(tout, expected)

})

test_that("MatchSigsAndRelabel 1", {
  gt.sigs <- matrix(c(1, 3,   4, 1, 2, 4), ncol = 2)
  ex.sigs <- matrix(c(1, 3.1, 4, 5, 1, 1, 1, 2.8, 4), ncol = 3)
  colnames(gt.sigs) <- c("gt.1", "gt.2")
  colnames(ex.sigs) <- c("ex.1", "ex.2", "ex.3")
  expected <-
    list(averCosSim = 0.890140132539195,
         match1 =
           structure(
             list(to = c("gt.1", "gt.2", "gt.1"),
                  sim = c(0.999877135226674,
                          0.461956578122389,
                          0.999473415175591)),
             class = "data.frame", row.names = c("ex.1",  "ex.2", "ex.3")),
         match2 =
           structure(
             list(to = c("ex.1", "ex.3" ),
                  sim = c(0.999877135226674, 0.989516398944649)),
             class = "data.frame", row.names = c("gt.1",  "gt.2")),
         extracted.with.no.best.match = "ex.2",
         ground.truth.with.no.best.match = "gt.2",
         ex.sigs = structure(
           c(1, 3.1, 4, 1, 2.8, 4, 5, 1, 1),
           .Dim = c(3L,      3L),
           .Dimnames = list(NULL, c(ex.1 = "ex.1 (gt.1 0.9999)",
                                    ex.3 = "ex.3 (gt.1 0.9995)",
                                    ex.2 = "ex.2 (gt.2 0.462)"))),
         gt.sigs = structure(
           c(1, 3, 4, 1, 2, 4),
           .Dim = 3:2,
           .Dimnames = list(         NULL, c("gt.1", "gt.2"))),
         gt.mean.cos.sim = list(gt.1 = 0.999675275201132,
                                gt.2 = 0.461956578122389))

  tout <- MatchSigsAndRelabel(gt.sigs = gt.sigs, ex.sigs = ex.sigs)
  expect_equal(tout, expected)

})

test_that("MatchSigsAndRelabel 2", {
  ex.sigs <- matrix(c(1, 3,   4, 1, 2, 4), ncol = 2)
  gt.sigs <- matrix(c(1, 3.1, 4, 5, 1, 1, 1, 2.8, 4), ncol = 3)
  colnames(ex.sigs) <- c("ex.1", "ex.2")
  colnames(gt.sigs) <- c("gt.1", "gt.2", "gt.3")
  expected <-
    list(averCosSim = 0.890140132539195,
         match1 =
           structure(list(to = c("gt.1", "gt.3"),
                          sim = c(0.999877135226674, 0.989516398944649     )),
                     class = "data.frame", row.names = c("ex.1", "ex.2")),
         match2 =
           structure(
             list(to = c("ex.1", "ex.2", "ex.1"),
                  sim = c(0.999877135226674,
                          0.461956578122389,
                          0.999473415175591)),
             class = "data.frame", row.names = c("gt.1",      "gt.2", "gt.3")),
         extracted.with.no.best.match = "ex.2",
         ground.truth.with.no.best.match = "gt.2",
         ex.sigs = structure(
           c(1,      3, 4, 1, 2, 4),
           .Dim = 3:2,
           .Dimnames =
             list(NULL, c(ex.1 = "ex.1 (gt.1 0.9999)",
                          ex.2 = "ex.2 (gt.3 0.9895) (gt.2 0.462)"))),
         gt.sigs = structure(
           c(1,      3.1, 4, 5, 1, 1, 1, 2.8, 4),
           .Dim = c(3L, 3L), .Dimnames = list(         NULL, c("gt.1", "gt.2", "gt.3"))),
         gt.mean.cos.sim = list(gt.1 = 0.999877135226674,
                                gt.2 = 0.461956578122389,
                                gt.3 = 0.989516398944649))

  tout <- MatchSigsAndRelabel(gt.sigs = gt.sigs, ex.sigs = ex.sigs)
  expect_equal(tout, expected)

})

if (FALSE) {
test_that("MatchSigsAndRelabel 3: 2 identical extracted signatures", {
  gt.sigs <- matrix(c(1, 3,   4, 
                      1, 2,   4, 
                      1, 3.1, 4), ncol = 3)
  ex.sigs <- matrix(c(1, 3.1, 4,
                      5, 1,   1,
                      1, 3.1, 4), ncol = 3)
  colnames(gt.sigs) <- c("gt.1", "gt.2", "gt.3")
  colnames(ex.sigs) <- c("ex.1", "ex.2", "ex.3")
  
  expected <-
    list(averCosSim = 0.907209327422042, 
         match1 = structure(list(
           to = c("gt.3", "gt.2", "gt.3"), 
           sim = c(1, 0.461956578122389,      1)), 
           class = "data.frame", row.names = c("ex.1", "ex.2",  "ex.3")), 
         match2 = structure(
           list(to = c("ex.1", "ex.1", "ex.1" ), 
                sim = c(0.999877135226674, 0.981422251183191, 1)),
           class = "data.frame", row.names = c("gt.1",  "gt.2", "gt.3")),
         extracted.with.no.best.match = c("ex.2", "ex.3" ),
         ground.truth.with.no.best.match = c("gt.1", "gt.2"),
         ex.sigs = structure(c(5,  1, 1, 1, 3.1, 4, 1, 3.1, 4),
                             .Dim = c(3L, 3L), 
                             .Dimnames = list(
                               NULL,
                               c(ex.2 = "ex.2 (gt.2 0.462)",
                                       ex.1 = "ex.1 (gt.3 1) (gt.1 0.9999)",   
                                       ex.3 = "ex.3 (gt.3 1)"))),
         gt.sigs = structure(c(1, 3, 4,  1, 2, 4, 1, 3.1, 4),
                             .Dim = c(3L, 3L), 
                             .Dimnames = list(
                               NULL,  
                               c("gt.1", "gt.2", "gt.3"))), 
         gt.mean.cos.sim = list(gt.1 = 0.999877135226674, 
                                gt.2 = 0.461956578122389,
                                gt.3 = 1))
  tout <- MatchSigsAndRelabel(gt.sigs = gt.sigs, ex.sigs = ex.sigs)
  expect_equal(tout, expected)
  
})


test_that("MatchSigsAndRelabel 4: 2 identical ground truth signatures", {
  gt.sigs <- matrix(c(1, 3, 4,
                      1, 2, 4,
                      1, 3, 4), ncol = 3)
  ex.sigs <- matrix(c(1, 3.1,  4,
                      5, 1,    1,
                      1, 3.11, 4), ncol = 3)
  colnames(gt.sigs) <- c("gt.1", "gt.2", "gt.3")
  colnames(ex.sigs) <- c("ex.1", "ex.2", "ex.3")
  expected <-
    list(averCosSim = 0.907143652143019,
         match1 = 
           structure(list(     to = c("gt.1", "gt.2", "gt.1"),
                               sim = c(0.999877135226674,      
                                       0.461956578122389,
                                       0.999851677872512)), 
                     class = "data.frame", 
                     row.names = c("ex.1",  "ex.2", "ex.3")),
         match2 = 
           structure(list(to = c("ex.1", "ex.1",  "ex.1"),
                          sim = c(0.999877135226674, 
                                  .981422251183191, 
                                  0.999877135226674 )),
                     class = "data.frame", 
                     row.names = c("gt.1", "gt.2", "gt.3" )), 
         extracted.with.no.best.match = c("ex.2", "ex.3"),
         ground.truth.with.no.best.match = c("gt.2",  "gt.3"),
         ex.sigs = structure(c(1, 3.1, 4, 1, 3.11, 4, 5, 1, 1 ),
                             .Dim = c(3L, 3L), 
                             .Dimnames = 
                               list(NULL,
                                    c(ex.1 = "ex.1 (gt.1 0.9999) (gt.3 0.9999)", 
                                      ex.3 = "ex.3 (gt.1 0.9999)", 
                                      ex.2 = "ex.2 (gt.2 0.462)"))),
         gt.sigs =
           structure(c(1,  3, 4, 1, 2, 4, 1, 3, 4), .Dim = c(3L, 3L),
                     .Dimnames =
                       list(NULL,      c("gt.1", "gt.2", "gt.3"))), 
         gt.mean.cos.sim = 
           list(gt.1 = 0.999864406549593, 
                gt.2 = 0.461956578122389,
                gt.3 = 0.999877135226674))
  tout <- MatchSigsAndRelabel(gt.sigs = gt.sigs, ex.sigs = ex.sigs)
  expect_equal(tout, expected)
  
})


test_that("MatchSigsAndRelabel 5: another stress test", {
  gt.sigs <- matrix(c(1, 3,
                      1, 3.001), ncol = 2)
  ex.sigs <- matrix(c(1, 3,
                      1, 3.1), ncol = 2)
  colnames(gt.sigs) <- c("g1", "g2")
  colnames(ex.sigs) <- c("x1", "x2")
  expected <- 1
  tout <- MatchSigsAndRelabel(gt.sigs = gt.sigs, ex.sigs = ex.sigs)
  # expect_equal(tout, expected)
  
})

test_that("MatchSigsAndRelabel 6: another stress test", {
  gt.sigs <- matrix(c(1, 3,
                      1, 2.9), ncol = 2)
  ex.sigs <- matrix(c(1, 3,
                      1, 3.1), ncol = 2)
  colnames(gt.sigs) <- c("g1", "g2")
  colnames(ex.sigs) <- c("x1", "x2")
  expected <- 1 
  tout <- MatchSigsAndRelabel(gt.sigs = gt.sigs, ex.sigs = ex.sigs)
  tout
  # expect_equal(tout, expected)
  
})

test_that("MatchSigsAndRelabel 7: another stress test", {
  gt.sigs <- matrix(c(1, 3,
                      1, 3.3), ncol = 2)
  ex.sigs <- matrix(c(1, 3,
                      1, 3.1), ncol = 2)
  colnames(gt.sigs) <- c("g1", "g2")
  colnames(ex.sigs) <- c("x1", "x2")
  expected <- 1 
  tout <- MatchSigsAndRelabel(gt.sigs = gt.sigs, ex.sigs = ex.sigs)
  tout
  # expect_equal(tout, expected)
  
})
}

