#include "IntervalTree.h"
#include "pthread-lite.h"
#include <vector>
#include <cassert>
#include <Rcpp.h>

#include <iostream>

#define rassert(b, msg)                        \
if (b == 0)                                         \
{                                                   \
  throw std::runtime_error(msg);                \
}

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
const Rcpp::IntegerVector *cloop;
const Rcpp::IntegerVector *sloop;
const Rcpp::IntegerVector *eloop;
bool global_verbose;
bool global_index_only;

// forward declare
bool find_overlaps(const GenomicIntervalTreeMap* tree, int32_t c, int32_t s, int32_t e,
                   size_t index,
                   std::vector<int32_t>* query, std::vector<int32_t>* subject,
                   std::vector<int32_t>* chr, std::vector<int32_t>* start,
                   std::vector<int32_t>* end);

/** Define a thread-item class to hold data, etc private
 * to each thread. For instance, this can store output from a thread that
 * can be dumped to a file after all processing is done. This is useful
 * because writing to a file in a multi-threaded program requires a mutex lock,
 * thus halting work on other threads. Alternatively useful for holding a pointer
 * for random access to a file, so multiple threads can randomly access the same file
 */
struct OverlapThreadItem {

  OverlapThreadItem(size_t i) : id(i) {}

  size_t id; // id to identify thread

  // include any number of thread-specific data below
  std::vector<int32_t> co, so, eo, query_id, subject_id;

  size_t NumHits() const {
    return co.size();
  }
};

/** Define a work-item class to hold data for specific task
 * (e.g. some operation on a set of sequences stored in char array)
 */
class OverlapWorkItem {

public:
  OverlapWorkItem(int start, int end) : m_start(start), m_end(end) {}

  // define the actual work to be done
  bool runOverlap(OverlapThreadItem* thread_data) {

    // loop through the query and overlap with subject
    for (size_t i = m_start; i < m_end; ++i) {
      find_overlaps(tree2, cloop->at(i),
                    sloop->at(i),
                    eloop->at(i),
                    i,
                    &thread_data->query_id, &thread_data->subject_id,
                    &thread_data->co, &thread_data->so, &thread_data->eo);
    }

    if (global_verbose)
      Rprintf("finished overlaps for intervals (%d,$d)\n",m_start,m_end);

    return true;
  }

  // always include a run function that takes only
  // a thread-item and returns bool
  bool run(OverlapThreadItem* thread_data) {
    // do the actual work
    return runOverlap(thread_data);
  }

private:

  // some chunk of data to be processed as one unit on one thread
  int m_start, m_end;
};

//' Construct the interval tree
//' @param c Vector of chromosomes (as integers)
//' @param s Vector of start positions
//' @param e Vector of end positions
//' @param tree GenomicIntervalTreeMap to fill (see https://github.com/walaj/SeqLib)
//' @noRd
void make_tree(const Rcpp::IntegerVector& c, const Rcpp::IntegerVector& s, const Rcpp::IntegerVector& e, GenomicIntervalTreeMap* tree) {

  GenomicIntervalMap map;
  for (size_t i = 0; i < c.size(); ++i) {
    //assert(i == 0 || (c.at(i-1) <= c.at(i)) || (s(i-1) < s(i))); //assure sorted
    map[c.at(i)].push_back(GenomicInterval(s(i), e(i), i));
  }

  // for each chr, make the tree from the intervals
  for (GenomicIntervalMap::iterator it = map.begin(); it != map.end(); ++it) {
    GenomicIntervalTreeMap::iterator ff = tree->find(it->first);
    if (ff != tree->end())
      ff->second = GenomicIntervalTree(it->second);
    else
      tree->insert(std::pair<int, GenomicIntervalTree>(it->first, GenomicIntervalTree(it->second)));
  }
}

bool find_overlaps(const GenomicIntervalTreeMap* tree, int32_t c, int32_t s, int32_t e,
                  size_t index,
                  std::vector<int32_t>* query, std::vector<int32_t>* subject,
                  std::vector<int32_t>* chr, std::vector<int32_t>* start,
                  std::vector<int32_t>* end) {

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
    if (!global_index_only) {
      chr->push_back(c);
      start->push_back(std::max(static_cast<int32_t>(j->start), static_cast<int32_t>(s)));
      end->push_back(std::min(static_cast<int32_t>(j->stop), static_cast<int32_t>(e)));
    }
  }
  return true;
}

void check_sort(const Rcpp::IntegerVector& c, const Rcpp::IntegerVector& s,
                const Rcpp::IntegerVector& e) {
  SeqHashSet<int32_t> chr_table;
  size_t sort_start = 0;
  int32_t curr_chr = 0;
  for (size_t i = 0; i < c.size(); ++i) {
    rassert(e[i] >= s[i], "end < start");
    if (chr_table.find(c[i]) == chr_table.end() || (i+1 == c.size())) { // not found or last one
      chr_table.insert(c[i]);
      rassert(std::is_sorted(s.begin() + sort_start, s.begin() + i + 1), "start pos not sorted");
      sort_start = i;
      curr_chr = c[i];
    } else {
      rassert(curr_chr == c[i], "seqnames not sorted"); // assure we are still on same chr
    }
  }
}

//' Perform the overlaps using an interval tree
//' @param df1 query data.table / data.frame with fields: seqnames, start, end
//' @param df2 subject data.table / data.frame with fields: seqnames, start, end
//' @param cores Max number of cores to use (process in 1,000,000 unit chunks)
//' @param verbose Print more
//' @param index_only Only return the index values (saves memory)
//' @return data.frame with ranges (seqnames, start, end) and query.id and subject.id
//' @noRd
// [[Rcpp::export]]
Rcpp::DataFrame cppoverlaps(const Rcpp::DataFrame& df1, const Rcpp::DataFrame& df2, int cores, bool verbose, bool index_only)
{
  if (verbose)
    Rprintf("start roverlaps.cpp");

  global_verbose = verbose;
  global_index_only = index_only;
  const Rcpp::IntegerVector c1 = df1["seqnames"];
  const Rcpp::IntegerVector c2 = df2["seqnames"];
  const Rcpp::IntegerVector s1 = df1["start"];
  const Rcpp::IntegerVector s2 = df2["start"];
  const Rcpp::IntegerVector e1 = df1["end"];
  const Rcpp::IntegerVector e2 = df2["end"];

  // check sorted
  if (verbose)
    Rprintf("roverlaps.cpp: checking sorting of input");

  check_sort(c1, s1, e1);
  check_sort(c2, s2, e2);

  // loop through and make the intervals for each chromosome
  tree2 = new GenomicIntervalTreeMap();

  // decide which is bigger, it gets looped. Smaller gets tree'ed
  // this is much faster and uses less memory
  // default is that 2 is the tree
  cloop = &c1;
  sloop = &s1;
  eloop = &e1;

  if (verbose)
    Rprintf("roverlaps: Making interval tree for %d interval\n",(c2.size() > c1.size() ? c1.size() : c2.size()));

  if (c2.size() > c1.size()) { // c1 should be tree
    make_tree(c1, s1, e1, tree2);
    cloop = &c2;
    sloop = &s2;
    eloop = &e2;
  } else {
    make_tree(c2, s2, e2, tree2);
  }

  if (cores==1) {
    std::vector<int32_t> co, so, eo, query_id, subject_id;
    if (verbose)
      Rprintf("robust: looping %d interval\n", cloop->size());

    for (int i = 0; i < cloop->size(); ++i)
      find_overlaps(tree2, cloop->at(i),
                    sloop->at(i),
                    eloop->at(i),
                    i,
                    &query_id, &subject_id,
                    &co, &so, &eo);
      if (c2.size() > c1.size()) // c1 should be tree
        std::swap(query_id, subject_id);
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

  // create the work item queue and consumer threads
  WorkQueue<OverlapWorkItem*>  queue; // queue of work items to be processed by threads
  std::vector<ConsumerThread<OverlapWorkItem, OverlapThreadItem>* > threadqueue;


  // add jobs to the thread
  int interval = 1000000;
  int end = -1;
  while (end < static_cast<int>(cloop->size())) {
    int start = end + 1;
    end = std::min(start + interval, static_cast<int>(cloop->size()));

    // add to the work queue
    OverlapWorkItem * wu = new OverlapWorkItem(start, end);
    queue.add(wu);
  }

  int num_cores = std::min(cores, (int)queue.size());

  if (verbose)
    Rprintf("roverlaps: sending %d work units with up to 1M intervals each on %d threads\n", queue.size(), num_cores);

  for (int i = 0; i < num_cores; ++i) {
    OverlapThreadItem * tu  = new OverlapThreadItem(i);  // create the thread-specific data, must be on heap.
    ConsumerThread<OverlapWorkItem, OverlapThreadItem>* threadr =        // establish new thread to draw from queue
      new ConsumerThread<OverlapWorkItem, OverlapThreadItem>(queue, tu); // always takes WorkQueue and some thread item

    threadr->start();
    threadqueue.push_back(threadr); // add thread to the threadqueue
  }

  // wait for the threads to finish
  for (int i = 0; i < num_cores; ++i)
    threadqueue[i]->join();

  size_t num_hits = 0;
  // get total num hits
  for (int i = 0; i < threadqueue.size(); ++i)
    num_hits += threadqueue[i]->GetThreadData()->NumHits();

  // merge it out
  if (verbose)
    Rprintf("roverlaps: allocing output memory for %d intervals\n",num_hits);

  std::vector<int32_t> co(num_hits);
  std::vector<int32_t> so(num_hits);
  std::vector<int32_t> eo(num_hits);
  std::vector<int32_t> query_id(num_hits);
  std::vector<int32_t> subject_id(num_hits);

  if (verbose)
    Rprintf("roverlaps: merging the output\n");

  size_t k = 0;
  for (int i = 0; i < threadqueue.size(); ++i) {
    const OverlapThreadItem * td = threadqueue[i]->GetThreadData();
    for (int j = 0; j < td->NumHits(); ++j) {
      co[k] = td->co.at(j);
      so[k] = td->so.at(j);
      eo[k] = td->eo.at(j);
      query_id[k] = td->query_id.at(j);
      subject_id[k] = td->subject_id.at(j);
      ++k;
    }
  }

  if (c2.size() > c1.size()) // c1 should be tree
    std::swap(query_id, subject_id);

  if (verbose)
    Rprintf("roverlaps.cpp: done with c++ call\n");

  return Rcpp::DataFrame::create(Rcpp::Named("seqnames")=co,
                                 Rcpp::Named("start")=so,
                                 Rcpp::Named("end")=eo,
                                 Rcpp::Named("query.id")=query_id,
                                 Rcpp::Named("subject.id")=subject_id);

}
