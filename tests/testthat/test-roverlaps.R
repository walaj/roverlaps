context("test-roverlaps.R")
library(data.table)
library(GenomicRanges)

test_that("unit test 1", {
  k=1e1
  verbose=FALSE
  index_only=FALSE
  o1 <- data.table(seqnames=factor(rep(c(1, "X"), each=k)), start=seq(k*2), end=seq(k*2)+2)
  k=1e1
  o2 <- data.table(seqnames=factor(rep(c(1, "X"), each=k)), start=seq(k*2), end=seq(k*2)+2)

  oc <- roverlaps(o1, o2, verbose=verbose, index_only=index_only)

  expect_equal(as.character(oc$seqnames)[1], as.character("1"))
  expect_equal(nrow(oc), 20)
  expect_equal(oc$start[1], 1)
  expect_equal(oc$end[1], 3)
  expect_equal(oc$query.id[1], 1)
  expect_equal(oc$subject.id[1], 1)
})

test_that("GenomicRanges convert", {
  o1 <- GRanges(c(1,1), IRanges(start=c(1,1), end=c(2,4)))
  o2 <- GRanges(c(1,1), IRanges(start=c(1,2), end=c(2,6)))
  o <- roverlaps(o1,o2)
  expect_equal(o$start[1], 1)
  expect_equal(o$end[1], 2)
  expect_equal(o$query.id[1], 1)
  expect_equal(o$subject.id[1], 1)
  expect_equal(nrow(o), 4)
})

test_that("verbose", { ## make sure it doesn't crash
  o1 <- data.table(seqnames=factor(c(1,1)), start=c(1,1), end=c(2,2))
  o2 <- data.table(seqnames=factor(c(1,2)), start=c(1,1), end=c(2,3))
  o <- roverlaps(o1, o2, verbose=TRUE)
})

test_that("index_only", { ## make sure it doesn't crash
  query <- data.table(seqnames=factor(c(1,1)), start=c(1,1), end=c(2,2))
  subject <- data.table(seqnames=factor(c(1,2)), start=c(1,1), end=c(2,3))
  o <- roverlaps(query, subject, index_only=TRUE)
})

test_that("same output whether or not sorted", {
  set.seed(42)
  C=10000
  o1 <- data.table(seqnames=factor(c(1:22, "X")[sample(seq(23),C, replace=TRUE)]), start=sample(seq(100000), C, replace=TRUE))
  o1[, end := start + 100]
  o2 <- data.table(seqnames=factor(c(1:22, "X")[sample(seq(23),C, replace=TRUE)]), start=sample(seq(100000), C, replace=TRUE))
  o2[, end := start + 100]
  o <- roverlaps(o1,o2)
  oi <- roverlaps(o1, o2, index_only=TRUE)
  expect(identical(o$query.id, oi$query.id))

  setkey(o1, seqnames, start)
  setkey(o2, seqnames, start)
  os <- roverlaps(o1,o2)
  ois <- roverlaps(o1, o2, index_only=TRUE)
  expect(all(os$start %in% o$start))
  expect(all(os$seqnames %in% o$seqnames))
  expect(all(ois$start %in% oi$start))
  expect(all(ois$seqnames %in% oi$seqnames))

  expect(all(o$start %in% os$start))
  expect(all(o$seqnames %in% os$seqnames))
  expect(all(oi$start %in% ois$start))
  expect(all(oi$seqnames %in% ois$seqnames))
})

test_that("test rcovered",{
 o1 <- data.table(seqnames=factor(c("1","2")),start=1,end=5)
 o2 <- data.table(seqnames=factor(c("1","1")),start=3,end=5)
 o1 <- o1[rcovered(o1,o2)]
 expect_equal(nrow(o1),1)
})

test_that("ragged diff interval", {
  a <- data.table(seqnames=c(2,2,1,1,1,"X"), start=c(3,7,10,67, 95, 1000))
  b <- data.table(seqnames=c(1,1,1), start=c(3,6,80))
  b[, end := start + 10] ## (3,13) (6,16) (80,90)

  expect_equal(rodiff(a, b),c(NA,NA,0,13, 5, NA))
  expect_equal(rodiff(a, b, sign=-1), c(NA,NA,0,13, NA,NA))   ## query must be <= interval
  expect_equal(rodiff(a, b, sign=1), c(NA,NA,0,51,5,NA)) ## query must be >= interval
  expect_equal(rodiff(data.table(), b), numeric(0))
  expect_equal(rodiff(a, data.table()), rep(NA, nrow(a)))
})

# test_that("test massive input", {
#   k=1e7
#   o1 <- data.table(seqnames=factor(rep(1,k)), start=1, end=2)
#   o2 <- data.table(seqnames=factor(2), start=1,end=3)
#   o <- roverlaps(o1,o2,verbose=TRUE,cores=4)
# })
