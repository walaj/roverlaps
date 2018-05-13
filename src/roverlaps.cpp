#include "IntervalTree.h"
#include <vector>
#include <cassert>
#include <unordered_set>
#include <unordered_map>
#include <Rcpp.h>

#include <iostream>

typedef TInterval<int32_t> GenomicInterval;
typedef std::unordered_map<int32_t, std::vector<GenomicInterval> > GenomicIntervalMap;
typedef TIntervalTree<int32_t> GenomicIntervalTree;
typedef std::unordered_map<int32_t, GenomicIntervalTree> GenomicIntervalTreeMap;
typedef std::vector<GenomicInterval> GenomicIntervalVector;

void make_tree(const Rcpp::IntegerVector& c, const Rcpp::IntegerVector& s, const Rcpp::IntegerVector& e, GenomicIntervalTreeMap& tree) {

  GenomicIntervalMap map;
  for (size_t i = 0; i < c.size(); ++i) {
    assert(i == 0 || (c.at(i-1) <= c.at(i)) || (s(i-1) < s(i))); //assure sorted
    map[c.at(i)].push_back(GenomicInterval(s(i), e(i), i));
  }

  // for each chr, make the tree from the intervals
  for (GenomicIntervalMap::iterator it = map.begin(); it != map.end(); ++it) {
    GenomicIntervalTreeMap::iterator ff = tree.find(it->first);
    if (ff != tree.end())
      ff->second = GenomicIntervalTree(it->second);
    else
      tree.insert(std::pair<int, GenomicIntervalTree>(it->first, GenomicIntervalTree(it->second)));
  }
}

// convert chr string (c) to chr numeric (o) in order they come, with mapping (map)
void chr_map(const Rcpp::IntegerVector& c, std::vector<int32_t> *o, std::unordered_map<int32_t, std::string>& map) {

  int chr_iter = -1; // always should find the first one so iterate up to 0;
  int j = 0;
  std::unordered_set<std::string> tmp_store;
  for (j = 0; j < c.size(); ++j) {
    if (tmp_store.find(std::to_string(c(j))) == tmp_store.end())
      map[++chr_iter] = std::to_string(c(j));
    o->at(j) = chr_iter;
  }
}

// [[Rcpp::export]]
Rcpp::DataFrame cppoverlaps(const Rcpp::DataFrame& df1, const Rcpp::DataFrame& df2)
{
  const Rcpp::IntegerVector c1 = df1["seqnames"];
  const Rcpp::IntegerVector c2 = df2["seqnames"];
  const Rcpp::IntegerVector s1 = df1["start"];
  const Rcpp::IntegerVector s2 = df2["start"];
  const Rcpp::IntegerVector e1 = df1["end"];
  const Rcpp::IntegerVector e2 = df2["end"];

  std::unordered_map<int32_t, std::string> c1_map, c2_map;
  std::vector<int32_t> co, so, eo, query_id, subject_id;

  // loop through and make the intervals for each chromosome
  GenomicIntervalTreeMap tree2;

  // decide which is bigger, it gets looped. Smaller gets tree'ed
  // this is much faster and uses less memory
  // default is that 2 is the tree
  const Rcpp::IntegerVector * cloop = &c1;
  bool query_is_query = true;
  if (c2.size() > c1.size()) { // c1 is tree
    make_tree(c1, s1, e1, tree2);
    cloop = &c2;
    query_is_query = false;
  } else {
    make_tree(c2, s2, e2, tree2);
  }

  // loop through the query and overlap with subject
  for (size_t i = 0; i < cloop->size(); ++i) {

    // which chr (if any) are common between query and subject
    GenomicIntervalTreeMap::const_iterator ff = tree2.find(cloop->at(i));
    GenomicIntervalVector giv;

    //must as least share a chromosome
    if (ff != tree2.end()) {

      // get the subject hits
      if (query_is_query)
        ff->second.findOverlapping(s1(i), e1(i), giv); // 	  ff->second.findOverlapping(m_grv->at(i), m_grv->at(i).pos2, giv);
      else
        ff->second.findOverlapping(s2(i), e2(i), giv);

      // loop through the hits and define the GenomicRegion
      // giv points to positions on subject
      for (GenomicIntervalVector::const_iterator j = giv.begin(); j != giv.end(); ++j) {
        query_id.push_back(i+1); // R is 1 indexed
        subject_id.push_back(j->value + 1);
        //co.push_back(c1_map[c1.at(i)]);
        co.push_back(cloop->at(i));
        if (query_is_query) {
          so.push_back(std::max(static_cast<int32_t>(j->start), static_cast<int32_t>(s1(i))));
          eo.push_back(std::min(static_cast<int32_t>(j->stop), static_cast<int32_t>(e1(i))));
        } else {
          so.push_back(std::max(static_cast<int32_t>(j->start), static_cast<int32_t>(s2(i))));
          eo.push_back(std::min(static_cast<int32_t>(j->stop), static_cast<int32_t>(e2(i))));
        }
      }
    }
  }

  if (!query_is_query) // switched query (1) and subject (2) for memory sake
    std::swap(query_id, subject_id);

  return Rcpp::DataFrame::create(Rcpp::Named("seqnames")=co,
                                 Rcpp::Named("start")=so,
                                 Rcpp::Named("end")=eo,
                                 Rcpp::Named("query.id")=query_id,
                                 Rcpp::Named("subject.id")=subject_id);

}
