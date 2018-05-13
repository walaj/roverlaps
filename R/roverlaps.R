#' roverlaps
#'
#' Description of your package
#'
#' @docType package
#' @author you <youremail>
#' @import Rcpp
#' @importFrom Rcpp evalCpp
#' @useDynLib roverlaps
#' @name roverlaps

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
        f = c('seqnames', 'start', 'end', 'strand', 'width')
        f2 = c('as.character(seqnames', 'c(start', 'c(end', 'as.character(strand', 'as.numeric(width')
        cmd = paste(cmd, paste(f, '=', f2, '(x))', sep = '', collapse = ','), sep = '')
        value.f = names(values(x))
      } else {
        was.gr = FALSE
        value.f = names(x)
      }

    if (length(value.f)>0)
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
#' @param o1 Query ranges
#' @param o2 Subject ranges
#' @param robust Perform extra checks (\code{is.unsorted}). \code{[TRUE]}
#' @importFrom data.table data.table as.data.table
#' @return data.table ('seqnames', 'start', 'end', 'strand', 'query.id', 'subject.id') of overlaps
#' @export
roverlaps <- function(o1, o2, robust=TRUE) {

  if (inherits(o1, "GRanges"))
    o1 <- cpp_gr2dt(o1)
  if (inherits(o2, "GRanges"))
    o2 <- cpp_gr2dt(o2)

  if (!inherits(o1,"data.table"))
    stop("required input is a data.table (preferred) or GRanges")
  if (!inherits(o2,"data.table"))
    stop("required input is a data.table (preferred) or GRanges")

  # convert char to factor for Rcpp
  charswitch = FALSE;
  if (class(o1$seqnames) != "factor") {
    o1$seqnames = as.factor(o1$seqnames)
    charswitch = TRUE
  }
  if (class(o2$seqnames) != "factor") {
    o2$seqnames = as.factor(o2$seqnames)
    charswitch = TRUE
  }

  ## check sorted
  if (robust) {
    if (any(o1[, is.unsorted(start), by="seqnames"]$V1))
      stop("o1 is not sorted by position")
    #if (o1[, is.unsorted(as.character(seqnames))])
    #  stop("o1 is not sorted by seqnames")

    if (any(o2[, is.unsorted(start), by="seqnames"]$V1))
      stop("o2 is not sorted by position")
    #if (o2[, is.unsorted(as.character(seqnames))])
    #  stop("o2 is not sorted by seqnames")

    if (any(o1[, end < start]))
      stop("o1 has end < start intervals")
    if (any(o2[, end < start]))
      stop("o1 has end < start intervals")
  }

  ## enforce that seqnames work between the two
  # (should be factor at this point, per above)
  new_levels <- union(levels(o1$seqnames),levels(o2$seqnames))
  if (robust) {
    o1$seqnames <- factor(o1$seqnames, levels=new_levels)
    o2$seqnames <- factor(o2$seqnames, levels=new_levels)
  }

  if (!identical(levels(o1$seqnames), levels(o2$seqnames)))
    stop("requires that o1 and o2 have same factor levels in same order. Or that o1 and o2 have same characters")

  ## do the actual overlaps
  o <- data.table::as.data.table(cppoverlaps(o1, o2))

  ## convert back from int to factor
  o$seqnames = factor(levels(o1$seqnames)[o$seqnames], levels=levels(o1$seqnames))

  return(o)
}
