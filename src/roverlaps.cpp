#include "IntervalTree.h"
#include <vector>
#include <Rcpp.h>
#include <limits>
#include <cassert>

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
typedef std::pair<int32_t, int32_t> ginterval;
typedef std::pair<std::int32_t, size_t> point;

// shared memory items across threads
GenomicIntervalTreeMap * tree2;

// Method to compare which one is the more close. 
// We find the closest by taking the difference 
// between the target and both values. It assumes 
// that val2 is greater than val1 and target lies 
// between these two. 
point getClosest(const point& val1, const point& val2, 
                 int32_t target) 
{ 
  if (target - val1.first >= val2.first - target) 
    return val2; 
  else
    return val1; 
} 


// Returns element closest to target in arr[] 
// sign = 1: q >= s
point findClosest(const std::vector<point>& arr, int32_t target, int sign) 
{ 
  // Corner cases 
  if (target <= arr[0].first && (sign == 0 || sign == -1)) 
    return arr[0]; 
  else if (target <= arr[0].first) // sign 1, q must be bigger
    return point(std::numeric_limits<int32_t>::min(),0);
  else if (target >= arr[arr.size() - 1].first && (sign == 0 || sign == 1)) 
    return arr[arr.size() - 1];
  else if (target >= arr[arr.size() - 1].first) // sign -1, q must be smaller, but impossible
    return point(std::numeric_limits<int32_t>::min(),0);
  
  // Doing binary search 
  int32_t i = 0, j = arr.size(), mid = 0; 
  while (i < j) { 
    mid = (i + j) / 2; 
    
    if (arr[mid].first == target) 
      return arr[mid]; 
    
    /* If target is less than array element, 
    then search in left */
    if (target < arr[mid].first) {   
      
      // If target is greater than previous 
      // to mid, return closest of two 
      if (mid > 0 && target > arr[mid - 1].first) { // --(m-1)----T----M
        if (sign==0)
          return getClosest(arr[mid - 1], 
                          arr[mid], target); 
        else if (sign == 1) // query is bigger than previous to mid, so return that one
            return(arr[mid - 1]);
        else  // sign: -1. query is the smaller one
            return(arr[mid]);
      }
      
      /* Repeat for left half */
      j = mid; 
    } 
    
    // If target is greater than mid 
    else { 
      if (mid < arr.size() - 1 && target < arr[mid + 1].first) {// --(m)----T----(m+1)
        if (sign==0)  
        return getClosest(arr[mid], 
                          arr[mid + 1], target); 
        else if (sign == 1) // query needs to be bigger, so return m
          return(arr[mid]);
        else // sign -1: query needs to be smaller
          return(arr[mid+1]);
      }
        
      // update i 
      i = mid + 1;  
    } 
  } 
  
  // Only single element left after search 
  return point(std::numeric_limits<int32_t>::min(),0); //arr[mid]; 
} 

// ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
// This code is contributed bu Smitha Dinesh Semwal 

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

struct less_than_key {
  inline bool operator() (const point& struct1, const point& struct2) {
    return (struct1.first < struct2.first);
  }
};

//' Perform the ragged difference between a vector and intervals
//' @param query Numeric vector to query 
//' @param subject_start Subject to query (start coordinates)
//' @param subject_end Subject to query (end coordinates)
//' @param max Return the max difference instead of min
//' @param sign If 0, consider values where q > s and s > q, if 1 only consider values where q >= s, if -1 only consider values where s >= q
//' @return Numeric vector of length same as query, of differences between query and subject
//' @noRd
// [[Rcpp::export]]
Rcpp::NumericVector cpprodiff(const Rcpp::DataFrame& query,
                                          const Rcpp::DataFrame& subject, 
                                          bool max, 
                                          int sign = 0)
{

  // define input structures
  const Rcpp::IntegerVector qseq = query["seqnames"];
  const Rcpp::IntegerVector sseq = subject["seqnames"];
  const Rcpp::IntegerVector qstart = query["start"];
  const Rcpp::IntegerVector sstart = subject["start"];
  const Rcpp::IntegerVector qend = query["end"];
  const Rcpp::IntegerVector send = subject["end"];
  
  Rcpp::NumericVector results(qseq.size());
  
  // if max, default everything to 0 so that unfilled values are discarded during max calc
  // if min, default everything to max, so that unfilled values are discarded during min calc
  float def = max ? 0 : std::numeric_limits<float>::max();
  
  typedef SeqHashMap<int32_t, std::vector<ginterval>> int_map;
  typedef SeqHashMap<int32_t, std::vector<point>> vmap;
  
  int_map map;
  vmap starts;
  vmap ends;
  
  // make the intervals
  for (size_t i = 0; i < sseq.size(); ++i) {
     map[sseq.at(i)].push_back(ginterval(sstart.at(i), send.at(i)));
     starts[sseq.at(i)].push_back(point(sstart.at(i), i));
     ends[sseq.at(i)].push_back(point(send.at(i), i));
  }
  
  // sort them
  for (auto i : starts)
    std::sort(i.second.begin(), i.second.end(), less_than_key());
  for (auto i : ends)
    std::sort(i.second.begin(), i.second.end(), less_than_key());
    
  for (size_t i = 0; i < qseq.size(); ++i) {
    const int32_t q = qstart.at(i);
    
    // if there is no seqnames overlap between query and subject, just continue
    int_map::const_iterator ii = map.find(qseq.at(i));
    if (ii == map.end()) {
      results[i] = -1;
      continue;
    }
    
    // check if it is contained.
    // if hits and is bigger than start and smaller than end its same end, and its same interval, then hit
    point nearest_start2 = findClosest(starts[qseq.at(i)], q, 1); // query must be bigger than start
    //point nearest_end2   = findClosest(ends[qseq.at(i)], q, -1); // query must be smaller than end
    //Rprintf("q: %d, ns2: %d, ne2: %d, %d, %d\n", q, nearest_start2.first, nearest_end2.first, nearest_start2.second, nearest_end2.second);
    
    //Rprintf("qsize: %d, q: %d, i %d, ns2: %d, send: %d\n", 
    //        qseq.size(), q, i, nearest_start2.first, send[i]);
    
    if ( (q >= nearest_start2.first && q <= send[nearest_start2.second]) && 
         nearest_start2.first > std::numeric_limits<int32_t>::min()) {
      results[i] = 0;
      continue;
    }
      
    // we know it's not contained, so binary way
    point nearest_start = findClosest(starts[qseq.at(i)], q, -1); // query must be smaller than start
    point nearest_end   = findClosest(ends[qseq.at(i)], q,  1); // query must be bigger than end
    
    //Rprintf("q: %d, ns: %d, ne: %d sign: %d\n", q, nearest_start, nearest_end, sign);
    
    // nothing hits
    if (nearest_end.first == std::numeric_limits<int32_t>::min() && 
        nearest_start.first == std::numeric_limits<int32_t>::min()) {
      results[i] = -1;
      continue;
    }
      
    // for sign 1, query has to be bigger than end, so if nothing hits on nearest_end, then failed
    if (sign == 1) {
      bool hit = nearest_end.first > std::numeric_limits<int32_t>::min();
      results[i] = hit ? std::abs(nearest_end.first - q) : -1;
      continue;
    } 
    
    // for sign -1, query has to be smaller than end, so if nothing hits on nearest_start, fail
    if (sign == -1) { 
      bool hit = nearest_start.first > std::numeric_limits<int32_t>::min();
      results[i] = hit ? std::abs(nearest_start.first - q) : -1;
      continue;
    }
    
    // just get the nearest
    if (sign == 0) {
      if (nearest_start.first == std::numeric_limits<int32_t>::min())
        results[i] = std::abs(nearest_end.first - q);
      else if (nearest_end.first == std::numeric_limits<int32_t>::min())
        results[i] = std::abs(nearest_start.first - q);
      else
        results[i] = std::min(std::abs(nearest_end.first - q), std::abs(nearest_start.first - q));
    }
    
    std::vector<float> tmp(map[qseq.at(i)].size(), def);
    
    bool any_result = false; // set to true if it found at least one hit
    for (size_t j = 0; j < ii->second.size(); ++j) {
      
      const int32_t subject_start = ii->second.at(j).first;
      const int32_t subject_end   = ii->second.at(j).second;
      
      // query contained in subject
      if (q >= subject_start && q <= subject_end) {
        tmp[j] = 0;
        any_result = true;
        continue;
      }
      
      if (sign == -1 && subject_start < q) // -1 require q <= s   ---q----|---s---|
        continue;
      if (sign == 1 && subject_end > q) // 1: require q >= s  ---|--s--|---q---
        continue;
      
      // query can't be equal to subject start or end, because of above
      // query also can't be inside |--s--| because of above
      tmp[j] = q > subject_end ? q - subject_end : subject_start - qstart[i];
      any_result = true;
    }
    results[i] = max ? *std::max_element(tmp.begin(), tmp.end()) : *std::min_element(tmp.begin(), tmp.end());
    results[i] = any_result ? results[i] : -1; // if no result, return -1 (to signal to convert to NA)
  }
  return results;
}
