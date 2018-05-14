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
  expect_equal(nrow(oc), 20)
  expect_equal(oc$start[1], 1)
  expect_equal(oc$end[1], 3)
  expect_equal(oc$query.id[1], 1)
  expect_equal(oc$subject.id[1], 1)
})

test_that("check sort fail", {
  o1 <- data.table(seqnames=factor(c(1,1)), start=c(2,1), end=c(3,2))
  o2 <- data.table(seqnames=factor(c(2,2)), start=c(2,1), end=c(3,2))
  expect_error(roverlaps(o1, o2, verbose=TRUE))
})

# test_that("test massive input", {
#   k=1e7
#   o1 <- data.table(seqnames=factor(rep(1,k)), start=1, end=2)
#   o2 <- data.table(seqnames=factor(2), start=1,end=3)
#   o <- roverlaps(o1,o2,verbose=TRUE,cores=4)
# })
