context("test-roverlaps.R")
library(data.table)
test_that("unit test 1", {
  k=1e1
  verbose=FALSE
  cores=1
  index_only=FALSE
  o1 <- data.table(seqnames=factor(rep(c(1, "X"), each=k)), start=seq(k*2), end=seq(k*2)+2)
  k=1e1
  o2 <- data.table(seqnames=factor(rep(c(1, "X"), each=k)), start=seq(k*2), end=seq(k*2)+2)

  oc <- roverlaps(o1, o2, cores=cores, verbose=verbose, index_only=index_only)

  expect_equal(as.character(oc$seqnames)[1], as.character("1"))
  expect_equal(nrow(oc), 88)
  expect_equal(oc$start[1], 1)
  expect_equal(oc$end[1], 3)
  expect_equal(oc$query.id[1], 1)
  expect_equal(oc$subject.id[1], 1)
})
