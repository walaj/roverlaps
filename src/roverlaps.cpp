#include "IntervalTree.h"
#include <vector>
#include <Rcpp.h>

#define rassert(b, msg)                       \
if ((b) == 0)                                 \
  throw std::runtime_error(msg);

#if __cplusplus > 199711L
#include <memory>
#include <unordered_set>
#include <unordered_map>
#define SeqHashMap std::unordered_map
#define SeqHashSet std::unordered_set
#define SeqPointer std::shared_ptr
#define HAVE_C11 1
#else

#ifdef __APPLE__
#include <memory>
#include <unordered_set>
#include <unordered_map>
#define SeqHashMap std::unordered_map
#define SeqHashSet std::unordered_set
#define SeqPointer std::shared_ptr
#else
#include <tr1/memory>
#include <tr1/unordered_set>
#include <tr1/unordered_map>
#define SeqHashMap std::tr1::unordered_map
#define SeqHashSet std::tr1::unordered_set
#define SeqPointer std::tr1::shared_ptr
#endif
#endif

typedef TInterval<int32_t> GenomicInterval;
typedef SeqHashMap<int32_t, std::vector<GenomicInterval> > GenomicIntervalMap;
typedef TIntervalTree<int32_t> GenomicIntervalTree;
typedef SeqHashMap<int32_t, GenomicIntervalTree> GenomicIntervalTreeMap;
typedef std::vector<GenomicInterval> GenomicIntervalVector;

// shared memory items across threads
GenomicIntervalTreeMap * tree2;

//' Construct the interval tree
//' @param c Vector of chromosomes (as integers)
//' @param s Vector of start positions
//' @param e Vector of end positions
//' @param tree GenomicIntervalTreeMap to fill (see https://github.com/walaj/SeqLib)
//' @noRd
void make_tree(const Rcpp::IntegerVector& c, const Rcpp::IntegerVector& s, const Rcpp::IntegerVector& e, GenomicIntervalTreeMap* tree) {

  GenomicIntervalMap map;
  for (size_t i = 0; i < c.size(); ++i)
    map[c.at(i)].push_back(GenomicInterval(s(i), e(i), i));

  // for each chr, make the tree from the intervals
  for (GenomicIntervalMap::iterator it = map.begin(); it != map.end(); ++it) {
    GenomicIntervalTreeMap::iterator ff = tree->find(it->first);
    if (ff != tree->end())
      ff->second = GenomicIntervalTree(it->second);
    else
      tree->insert(std::pair<int, GenomicIntervalTree>(it->first, GenomicIntervalTree(it->second)));
  }
}

//' Calculate if a range overlaps against an interval tree
//' @param tree Interval tree (map) to query against
//' @param c Query chromosome
//' @param s Query start
//' @param e Query end
//' @param index Position of the query range in the query set
//' @param query Will store the index of the query range
//' @param subject Will store the index of the subject range
//' @param chr Will store the chromosome of the output overlap range
//' @param start Will store the start of the output overlap range
//' @param end Will store the end of the output overlap range
//' @noRd
bool find_overlaps(const GenomicIntervalTreeMap* tree, int32_t c, int32_t s, int32_t e,
                   size_t index,
                   std::vector<int32_t>* query, std::vector<int32_t>* subject,
                   std::vector<int32_t>* chr, std::vector<int32_t>* start,
                   std::vector<int32_t>* end, bool index_only) {

  // which chr (if any) are common between query and subject
  GenomicIntervalTreeMap::const_iterator ff = tree->find(c);
  GenomicIntervalVector giv;

  //must as least share a chromosome
  if (ff == tree->end())
    return false;

  // get the subject hits
  ff->second.findOverlapping(s, e, giv);

  // loop through the hits and define the GenomicRegion
  // giv points to positions on subject
  for (GenomicIntervalVector::const_iterator j = giv.begin(); j != giv.end(); ++j) {
    query->push_back(index + 1); // R is 1 indexed
    subject->push_back(j->value + 1);
    if (!index_only) {
      chr->push_back(c);
      start->push_back(std::max(static_cast<int32_t>(j->start), static_cast<int32_t>(s)));
      end->push_back(std::min(static_cast<int32_t>(j->stop), static_cast<int32_t>(e)));
    }
  }
  return true;
}

//' Check if a set of ranges is sorted
//'
//' Will throw and error if seqnames are not relatively sorted (but does not check
//' if absolutely sorted relative to ordering e.g. alphanumeric)
//' or will throw error if start position not sorted or will throw error if
//' range has negative (end < start) width.
//' SORTING ACTUALLY NOT REQUIRED
//' @param c Seqnames
//' @param s start positions
//' @param e end positions
//' @noRd
/*
void check_sort(const Rcpp::IntegerVector& c, const Rcpp::IntegerVector& s,
                const Rcpp::IntegerVector& e) {
  SeqHashSet<int32_t> chr_table;
  size_t sort_start = 0;
  int32_t curr_chr = 0;
  for (size_t i = 0; i < c.size(); ++i) {
    rassert(e[i] >= s[i], "end < start");
    if (chr_table.find(c[i]) == chr_table.end()) { // switch chr, so check the previous chr one
      chr_table.insert(c[i]);
      rassert(std::is_sorted(s.begin() + sort_start, s.begin() + i), "start pos not sorted");
      sort_start = i;
      curr_chr = c[i];
    } else {
      rassert(curr_chr == c[i], "seqnames not sorted"); // assure we are still on same chr
    }
  }
  // check the last chromsome
  rassert(std::is_sorted(s.begin() + sort_start, s.end()), "start pos not sorted");
}
*/

//' Perform the overlaps using an interval tree
//' @param df1 query data.table / data.frame with fields: seqnames, start, end
//' @param df2 subject data.table / data.frame with fields: seqnames, start, end
//' @param verbose Print more
//' @param index_only Only return the index values (saves memory)
//' @return data.frame with ranges (seqnames, start, end) and query.id and subject.id
//' @noRd
// [[Rcpp::export]]
Rcpp::DataFrame cppoverlaps(const Rcpp::DataFrame& df1, const Rcpp::DataFrame& df2, bool verbose, bool index_only)
{
  if (verbose)
    Rprintf("start roverlaps.cpp");

  // define input structures
  const Rcpp::IntegerVector c1 = df1["seqnames"];
  const Rcpp::IntegerVector c2 = df2["seqnames"];
  const Rcpp::IntegerVector s1 = df1["start"];
  const Rcpp::IntegerVector s2 = df2["start"];
  const Rcpp::IntegerVector e1 = df1["end"];
  const Rcpp::IntegerVector e2 = df2["end"];

  // define output structures
  std::vector<int32_t> co, so, eo, query_id, subject_id;

  // check sorted
  if (verbose)
    Rprintf("roverlaps.cpp: checking sorting of input");
  //check_sort(c1, s1, e1);
  //check_sort(c2, s2, e2);

  // loop through and make the intervals for each chromosome
  tree2 = new GenomicIntervalTreeMap();

  if (verbose)
    Rprintf("roverlaps: Making interval tree for %d interval\n",c2.size() > c1.size() ? c1.size() : c2.size());

  // decide which is bigger, it gets looped. Smaller gets tree'ed
  // this is much faster and uses less memory
  // default is that 2 is the tree
  if (c2.size() > c1.size()) { // c1 should be tree
    make_tree(c1, s1, e1, tree2);
    if (verbose)
      Rprintf("robust: looping %d interval\n", c2.size());
    for (int i = 0; i < c2.size(); ++i)
      find_overlaps(tree2,c2.at(i),s2.at(i),e2.at(i),i,&subject_id, &query_id,&co, &so, &eo,index_only);
  } else {
    make_tree(c2, s2, e2, tree2);
    if (verbose)
      Rprintf("robust: looping %d interval\n", c1.size());
    for (int i = 0; i < c1.size(); ++i)
      find_overlaps(tree2,c1.at(i),s1.at(i),e1.at(i),i,&query_id, &subject_id,&co, &so, &eo,index_only);
  }

  if (verbose)
    Rprintf("roverlaps.cpp: done with c++ call");

  delete tree2;

  if (index_only)
    return Rcpp::DataFrame::create(Rcpp::Named("query.id")=query_id,
                                   Rcpp::Named("subject.id")=subject_id);

  return Rcpp::DataFrame::create(Rcpp::Named("seqnames")=co,
                                 Rcpp::Named("start")=so,
                                 Rcpp::Named("end")=eo,
                                 Rcpp::Named("query.id")=query_id,
                                 Rcpp::Named("subject.id")=subject_id);
}
