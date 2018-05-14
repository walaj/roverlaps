[![Build Status](https://travis-ci.org/walaj/bxtools.svg?branch=master)](https://travis-ci.org/walaj/bxtools)

## *rovertools* - Fast and memory efficient genomic range overlaps in R

**License:** [MIT][license]

  * [Installation](#installation)
  * [Description](#description)
  * [Examples](#examples)
  * [Performance](#performance)
  * [Attributions](#attributions)

Installation
------------
```R
devtools::install_github("walaj/rovertools")
```

Description
-----------
Range overlaps between different sets of genomic intervals
	     can be memory intensive and slow for huge (1M+) interval queries.
	     Depending on the application, only part of a range overlap is needed 
	     (e.g. just the overlap mappings, not the actual ranges). This package
	     implements fast, memory-efficient interval tree overlaps in C++
	     and provides flexibility for users to choose only what they need,
	     thereby improving memory performance.

Examples
--------
#### Get overlap intervals

```R
library(data.table)
k=1
d1 <- data.table(seqnames=factor(rep(c(1, "X"), each=k)), start=seq(k*2), end=seq(k*2)+2)
k=2
d2 <- data.table(seqnames=factor(rep(c(1, "X"), each=k)), start=seq(k*2), end=seq(k*2)+2)

o <- roverlaps::roverlaps(d1, d2)
```

Will give the following output:
```R
   seqnames start end query.id subject.id
1:        1     1   3        1          1 
2:        1     2   3        1          2
3:        X     3   4        2          3
4:        X     4   4        2          4
```

Performance
-----------

Test of range overlap of 5 million vs 50 million records, repeated 10 times. Each 
overlap produces 24,999,997 overlaps.

Comparison is between ``roverlaps`` and ``gUtils::gr.findoverlaps`` which
is based on ``GenomicRanges::findOverlaps``.

<img src="https://github.com/walaj/roverlaps/blob/master/memgraph.both.png"
width=600/>

Attributions
------------

This project is developed and maintained by Jeremiah Wala (jwala@broadinstitute.org)

Thanks to Erik Garrison for the excellent [interval tree][tree] implementation.

[license]: https://github.com/walaj/rovertools/blob/master/LICENSE
[tree]: https://github.com/ekg/intervaltree

