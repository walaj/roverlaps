#' roverlaps
#'
#' Range overlaps between different sets of genomic intervals
#' can be memory intensive and slow for huge (1M+) interval queries.
#' This package implements fast, memory-efficient interval tree overlaps in C++
#' and provides flexibility for users to choose only what they need,
#' thereby improving memory performance.
#'
#' @docType package
#' @author Jeremiah Wala <jwala@broadinstitute.org>
#' @import Rcpp
#' @importFrom Rcpp evalCpp
#' @useDynLib roverlaps
#' @name roverlaps

##if(getRversion() >= "2.15.1")  utils::globalVariables(c(".", ":="))

#' @name cpp_gr2dt
#' @title Converts \code{GRanges} to \code{data.table}
#' @description
#'
#' Converts \code{GRanges} to \code{data.table}
#' and a field grl.iix which saves the (local) index that that gr was in its corresponding grl item
#'
#' @param x \code{GRanges} to convert
#' @import GenomicRanges
#' @importFrom S4Vectors elementNROWS
#' @importFrom data.table data.table as.data.table
#' @importFrom methods is
#' @return data.table of GRanges columns ('seqnames', 'start', 'end', 'strand', 'width') and metadata columns
#' @noRd
cpp_gr2dt = function(x)
  {
    ## new approach just directly instantiating data table
    cmd = 'data.frame(';

    if (methods::is(x, 'GRanges'))
      {
        ## as.data.table complains if duplicated row names
        if (any(duplicated(names(x)))){
          names(x) <- NULL
        }

        was.gr = TRUE
        f = c('seqnames', 'start', 'end', 'strand') #, 'width')
        f2 = c('as.character(seqnames', 'c(start', 'c(end', 'as.character(strand') #, 'as.numeric(width')
        cmd = paste(cmd, paste(f, '=', f2, '(x))', sep = '', collapse = ','), sep = '')
        value.f = names(values(x))
      } else {
        was.gr = FALSE
        value.f = names(x)
      }

    if (FALSE && length(value.f)>0)
      {
        if (was.gr){
          cmd = paste(cmd, ',', sep = '')
        }
        class.f = sapply(value.f, function(f) eval(parse(text=sprintf("class(x$'%s')", f))))

        .StringSetListAsList = function(x){
          tmp1 = as.character(unlist(x))
          tmp2 = rep(1:length(x), S4Vectors::elementNROWS(x))
          return(split(tmp1, tmp2))
        }

        ## take care of annoying S4 / DataFrame / data.frame (wish-they-were-non-)issues
        as.statement = ifelse(grepl('Integer', class.f), 'as.integer',
          ifelse(grepl('Character', class.f), 'as.character',
                 ifelse(grepl('((StringSet)|(Compressed.*))List', class.f), '.StringSetListAsList',
                        ifelse(grepl('StringSet$', class.f), 'as.character',
                               ifelse(grepl('factor$', class.f), 'as.character',
                                      ifelse(grepl('List', class.f), 'as.list',
                                             ifelse(grepl('factor', class.f), 'as.character',
                                                    ifelse(grepl('List', class.f), 'as.list', 'c'))))))))
        cmd = paste(cmd, paste(value.f, '=', as.statement, "(x$'", value.f, "')", sep = '', collapse = ','), sep = '')
      }

    cmd = paste(cmd, ')', sep = '')

    out = tryCatch(data.table::as.data.table(eval(parse(text =cmd))), error = function(e) NULL)

    if (is.null(out)){
      out = data.table::as.data.table(x)
    }

    return(out)
}

#' @name roverlaps
#' @title Fast overlaps of \code{data.table} or \code{GRanges}
#' @description
#'
#' Performs interval overlaps between two genomic ranges, returning the
#' intersecting set of ranges and the indicies of the query and subject
#' which created the interval.
#'
#' @note Positions in ranges are inclusive. Example: chr1:2-5 (query) and
#' chr1:5-7 (subject) will create an overlap with value chr1:5-5.
#'
#' @param query Query ranges as a \code{data.table} with mandatory fields \code{seqnames} and \code{start}
#' @param subject Subject ranges as a \code{data.table} with mandatory fields \code{seqnames} and \code{start}
#' @param verbose Increase the verbosity \code{[FALSE]}
#' @param index_only Return only the indicies ('query.id' and 'subject.id') \code{[FALSE]}
#' @importFrom data.table data.table as.data.table setkey set
#' @importFrom utils globalVariables
#' @return data.table ('seqnames', 'start', 'end', 'query.id', 'subject.id') of overlaps
#' @examples
#'
#' library(data.table)
#' set.seed(42)
#' C=10000
#' sn1 <- factor(c(1:22, "X")[sample(seq(23),C, replace=TRUE)])
#' s1  <- sample(seq(100000), C, replace=TRUE)
#' o1 <- data.table(seqnames=sn1, start=s1)
#' o1[, end := start + 100]
#' sn2 <- factor(c(1:22, "X")[sample(seq(23),C, replace=TRUE)])
#' s2  <- sample(seq(100000), C, replace=TRUE)
#' o2 <- data.table(seqnames=sn2, start=s2)
#' o2[, end := start + 100]
#' o <- roverlaps(o1,o2)
#'
#' ## output
#' # seqnames start   end query.id subject.id
#' #        7 83405 83451        3       3315
#' #        7 83405 83463        3       3148
#' #        7 83485 83505        3       1022
#'
#' oi <- roverlaps(o1, o2, index_only=TRUE)
#'
#' ## output
#' #    query.id subject.id
#' #           3       3315
#' #           3       3148
#' #           3       1022
#' @export
roverlaps <- function(query, subject, verbose=FALSE, index_only=FALSE) {

  if (verbose)
    print("roverlaps.R: checking input")

  if (inherits(query, "GRanges"))
    query <- cpp_gr2dt(query)
  if (inherits(subject, "GRanges"))
    subject <- cpp_gr2dt(subject)

  if (!inherits(query,"data.table"))
    stop("required input is a data.table (preferred) or GRanges")
  if (!inherits(subject,"data.table"))
    stop("required input is a data.table (preferred) or GRanges")

  # need full bounds (for now)
  if (!"end" %in% colnames(query))
    query$end = query$start
  if (!"end" %in% colnames(subject))
    subject$end = subject$start

  stopifnot(all(c("seqnames", "start","end") %in% colnames(query)))
  stopifnot(all(c("seqnames", "start","end") %in% colnames(subject)))

  ## fix issue where c++ doesn't like identical structs.
  ## skip ahead if trivial
  if (identical(query, subject)) {
    if (verbose)
      print("skipping cpp because input identical. Trival result")
    o <- data.table::copy(query[,list(seqnames, start, end)])
    o$subject.id <- o$query.id <- seq(nrow(o))
    return(o)
  }

  if (verbose)
    print("roverlaps.R: factorizing")

  # convert char to factor for Rcpp
  charswitch = FALSE;
  if (class(query$seqnames) != "factor") {
    query$seqnames = as.factor(query$seqnames)
    charswitch = TRUE
  }
  if (class(subject$seqnames) != "factor") {
    subject$seqnames = as.factor(subject$seqnames)
    charswitch = TRUE
  }

  ## enforce that seqnames work between the two
  # (should be factor at this point, per above)
  new_levels <- union(levels(query$seqnames),levels(subject$seqnames))
  needs_fix = !identical(levels(query$seqnames), levels(subject$seqnames)) ||
    class(query$seqnames) != "factor" || class(subject$seqnames) != "factor"

  if (needs_fix) { ## if index only, don't need to reset factors
    if (verbose)
      print("roverlaps.R: setting new factor levels")
    data.table::set(query, j="seqnames", value=factor(query$seqnames, levels=new_levels))
    data.table::set(subject, j="seqnames", value=factor(subject$seqnames, levels=new_levels))
    if (!identical(levels(query$seqnames), levels(subject$seqnames)))
      stop("query and subject must have same factor levels.")
  }

  if (verbose)
    print("roverlaps.R: calling cpp")

  ## do the actual overlaps
  stopifnot(all(c("seqnames", "start","end") %in% colnames(query)))
  stopifnot(all(c("seqnames", "start","end") %in% colnames(subject)))
  o <- as.data.table(cppoverlaps(query, subject, verbose, index_only))

  if (index_only)
    return (o)

  ## convert back from int to factor
  stopifnot(identical(levels(query$seqnames),new_levels))
  data.table::set(o, j="seqnames", value=factor(new_levels[o$seqnames], levels=new_levels))
  stopifnot(identical(levels(o$seqnames), new_levels))
  ##o[,seqnames := factor(levels(query$seqnames)[seqnames], levels=levels(query$seqnames))]

  return(o)
}

#' @name rcovered
#' @title Return logical of whether a query region intersects subject
#' @description
#'
#' Determines which members of a query set of regions has any overlap with
#' a set of subject intervals. The query does not have to be fully contained within
#' the subject to be counted as \code{TRUE} overlap.
#'
#' @param query Set of query intervals to evaluate for membership in \code{subject}
#' @param subject Set of subject intervals to check \code{query} against
#' @param verbose Set the verbosity \code{[FALSE]}
#' @return Logical vector of same length as \code{query}, \code{TRUE} if interval overlaps subject
#' @export
#'
#' @examples
#'
#' library(data.table)
#' o1 <- data.table(seqnames=factor(c("1","2")),start=1,end=5)
#' o2 <- data.table(seqnames=factor(c("1","1")),start=3,end=5)
#'
#' ## useful for subsetting in one line
#' o1 <- o1[rcovered(o1,o2)]
#'
#' ## output
#' #    seqnames start end
#' # 1:        1     1   5
rcovered <- function(query, subject, verbose=FALSE) {
  o <- roverlaps(query, subject, index_only=TRUE, verbose=verbose)
  if (!nrow(o))
    return(rep(FALSE,nrow(query)))
  return(seq(nrow(query)) %in% o$query.id)
}

#' @name raggeddiff
#' @title Find the minimum/maximum difference between each element of vector A, compared with vector B of any size 
#' @description
#'
#' For each element of a numeric vector A, this function will find the distance to either the closest or
#' further element in B.
#'
#' @param query Numeric vector to query 
#' @param subject Numeric vector to query against
#' @param max Take the max difference instead of min \code{[FALSE]}
#' @importFrom data.table data.table as.data.table setkey set
#' @importFrom utils globalVariables
#' @return Numeric vector of same length as query, holding the differences
#' @examples
#'
#' library(data.table)
#' set.seed(42)
#' query <- sample(100, 3)
#' subject <- sample(1000, 100)
#' o <- raggediff(query, subject)
#' @export
raggeddiff <- function(query, subject, max=FALSE) {

  cppraggeddiff(query, subject, max);

}
