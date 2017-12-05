#ifndef GA_HPP
#define GA_HPP
#define _GLIBCXX_USE_CXX11_ABI 0
#include <algorithm>
#include <boost/algorithm/string.hpp>
//#include <boost/archive/binary_iarchive.hpp>
//#include <boost/archive/binary_oarchive.hpp>
#include "../rlutil/rlutil.h"
#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/flyweight.hpp>
#include <boost/format.hpp>
#include <boost/lexical_cast.hpp>
#include <cmath>
#include <fstream>
#include <iostream>
#include <memory>
#include <queue>
#include <set>
#include <string.h>
#include <string>
#include <unordered_map>
#include <vector>
namespace ga {
typedef std::unordered_map<std::string, unsigned> map_;
typedef std::unordered_map<std::string, std::vector<unsigned>> map2_;
typedef boost::iterator_range<std::string::const_iterator> string_view;
typedef std::shared_ptr<map_> map;
typedef std::shared_ptr<map2_> map2;
// typedef boost::tokenizer<boost::escaped_list_separator<char>> Tokenizer;

typedef enum { sym = 0, annot, scalar, vector } dType;
typedef enum { ignore = 0, maskOnly, getData } passthroughType;
template <typename T> using Pt = std::shared_ptr<std::vector<T>>;
template <typename T> using Pool = std::queue<std::shared_ptr<std::vector<T>>>;

// basic annotation types
struct vectors {
  unsigned nc;
  std::vector<double> v;
};
struct inv {
  std::string chr;
  unsigned start, end;
  bool strand;
  std::string seq;
  bool operator==(const inv &another) {
    return ((chr == another.chr) && (start == another.start) &&
            (end == another.end) && (strand == another.strand));
  };
  bool operator!=(const inv &another) {
    return !((chr == another.chr) && (start == another.start) &&
             (end == another.end) && (strand == another.strand));
  };
};
struct interval {
  // note that intervals are "unsigned"
  // and only count absolute pos on + strand
  boost::flyweight<std::string> chr;
  boost::flyweight<unsigned long> sz;
  unsigned start = 0, end = 0, idx = 0;
  bool strand;
  inv v() { return inv{chr, start, end, strand}; };
  std::string info() {
    boost::format fmt("[%1%] tx: %2%-%3% s:%4%");
    return (fmt % chr % start % end % strand).str();
  };
  bool operator==(const interval &another) {
    return ((chr == another.chr) && (start == another.start) &&
            (end == another.end) && (strand == another.strand));
  };
  bool operator!=(const interval &another) {
    return !((chr == another.chr) && (start == another.start) &&
             (end == another.end) && (strand == another.strand));
  };
};

struct intervalCmp {
  bool operator()(const interval &l, const interval &r) {
    if (l.chr != r.chr)
      return l.chr < r.chr;
    if (r.start > l.end)
      return true;
    if (r.end < l.start)
      return false;
    return (l.start<r.start);
  };
};

struct gene {
  unsigned cds_start, cds_end, tx_start, tx_end, idx;
  std::string sym, id = "";
  boost::flyweight<std::string> chr;
  boost::flyweight<unsigned long> sz;
  bool strand;
  std::multiset<interval, intervalCmp> exons;
  interval inv() { return interval{chr, sz, tx_start, tx_end, idx, strand}; };
  interval cds() { return interval{chr, sz, cds_start, cds_end, idx, strand}; };
  std::string info() {
    boost::format fmt("[%1%] tx: %2%-%3%, cds: %4%-%5%, s:%6%, e:%7%");
    return (fmt % chr % tx_start % tx_end % cds_start % cds_end % strand %
            exons.size())
        .str();
  };
  std::vector<interval> get_exons() {
    std::vector<interval> out;
    for (auto &it : exons) {
      out.push_back(it);
    }
    return out;
  };
  std::vector<interval> get_introns() {
    std::vector<interval> out;
    if (exons.size() < 2)
      return out; // no intron in this situation
    auto it = exons.begin();
    auto itt = ++it;
    --it;
    while (itt != exons.end()) {
      interval i;
      i.chr = chr;
      i.sz = sz;
      i.idx = idx;
      i.strand = strand;
      i.start = it->end + 1;
      i.end = itt->start - 1;
      out.push_back(i);
      ++it;
      ++itt;
    }
    return out;
  };
  interval utr(const bool s) {
    interval i;
    i.chr = chr;
    i.sz = sz;
    i.idx = idx;
    i.strand = strand;
    if (s) {
      i.start = tx_start;
      i.end = cds_start - 1;
    } else {
      i.start = cds_end + 1;
      i.end = tx_end;
    }
    return i;
  };
  interval utr5() {
    if (noncoding()) {
      //std::cerr << "Gene at " << idx << " is noncoding, no 5'UTR" << std::endl;
      return interval();
    }
    return utr(strand);
  };
  interval utr3() {
    if (noncoding()) {
      //std::cerr << "Gene at " << idx << " is noncoding, no 3'UTR" << std::endl;
      return interval();
    }
    return utr(!strand);
  };
  interval get_p(const int &l, const int &r, const bool &s) {
    interval i;
    i.chr = chr;
    i.sz = sz;
    i.idx = idx;
    i.strand = strand;
    if (s) {
      i.start = tx_start - l;
      i.end = tx_start + r;
    } else {
      i.start = tx_end + l;
      i.end = tx_end - r;
    }
    return i;
  };
  interval get_promoter(const int &l, const int &r) {
    if (noncoding()) {
      std::cerr << "Gene at " << idx << " is noncoding, no promoter"
                << std::endl;
      return interval();
    }
    return get_p(l, r, strand);
  };
  interval get_tail(const int &l, const int &r) {
    if (noncoding()) {
      std::cerr << "Gene at " << idx << " is noncoding, no promoter"
                << std::endl;
      return interval();
    }
    return get_p(l, r, !strand);
  };
  virtual bool noncoding() {
    if (cds_start == cds_end) {
      return true;
    }
    return false;
  };
  bool operator==(const gene &another) { return (idx == another.idx); };
  bool operator!=(const gene &another) { return (idx != another.idx); };
};

// a collection of data and annotation
class Base {
public:
  Base() : coreIDs(std::make_shared<map_>(map_())){};
  /* Constructors */
  template <typename T> Pt<T> const add() {
    auto in = std::make_shared<std::vector<T>>();
    in->reserve(maxEle);
    return in;
  };
  // we use concrete types instead of base class ptrs to
  // increase performance somewhat
  map coreIDs;
  unsigned maxEle;
};

// genome based data source
class BaseGenome : public Base {
public:
  using Base::Base;
  BaseGenome() : GeneSym(std::make_shared<map2_>(map2_())){};
  virtual void add_entry(std::vector<string_view> &) = 0;
  virtual void load(const std::string &fname) = 0;
  std::vector<unsigned> find_sym(const std::string &in) {
    return GeneSym->at(in);
  };
  unsigned find_id(const std::string &in) { return coreIDs->at(in); };
  virtual unsigned long size(const std::string &in) { return sizes[in]; };
  virtual std::vector<std::string> chroms() {
    std::vector<std::string> out;
    for (auto &it : sizes)
      out.push_back(it.first);
    return out;
  };
  virtual std::vector<unsigned> find_ids(const std::vector<std::string> &);
  virtual std::vector<unsigned> find_syms(const std::vector<std::string> &);
  virtual std::vector<unsigned> r(const std::string &chr, const int &l,
                                  const int &r);
  virtual gene closest(const std::string &chr, const int &l, const int &r);
  virtual gene closest_right(const std::string &chr, const int &l, const int &r);
  virtual gene closest_left(const std::string &chr, const int &l, const int &r);
  virtual gene *const next_gene(const unsigned &idx);
  virtual gene *const prev_gene(const unsigned &idx);
  virtual interval intergene_up(const unsigned &idx);
  virtual interval intergene_down(const unsigned &idx);
  std::vector<unsigned> all() {
    std::vector<unsigned> out;
    out.reserve(maxEle);
    for (unsigned int i = 0; i < maxEle; ++i)
      out.push_back(i);
    return out;
  };
  unsigned count() { return maxEle; };
  virtual std::multiset<interval, intervalCmp>::iterator begin() {
    return GeneLoci.begin();
  };
  virtual std::multiset<interval, intervalCmp>::iterator end() {
    return GeneLoci.end();
  };
  virtual std::vector<unsigned> iter_genes() {
    std::vector<unsigned> out;
    out.reserve(maxEle);
    for (auto it = begin(); it != end(); ++it) {
      out.push_back(it->idx);
    }
    return out;
  };
  virtual std::vector<unsigned> iter_genes_unique() {
    std::vector<unsigned> out;
    out.reserve(maxEle);
    auto it = begin();
    while (it != end()) {
      it = GeneLoci.upper_bound(*it);
      if (it!=end())
          out.push_back(it->idx);
    }
    return out;
  };
  virtual gene get_gene(const unsigned &idx) { return Genes[idx]; };
  virtual std::vector<gene> get_genes(const std::vector<unsigned> &in) {
    std::vector<gene> out;
    for (auto &i : in) {
      if (i > maxEle)
        continue;
      out.push_back(Genes[i]);
    }
    return out;
  };
  BaseGenome *const get_this() { return this; };
  map2 GeneSym;
  std::unordered_map<std::string, unsigned long> sizes;
  std::multiset<interval, intervalCmp> GeneLoci;
  std::vector<gene> Genes;
};

class UCSCRefGene : public BaseGenome {
public:
  UCSCRefGene(){};
  void add_entry(std::vector<string_view> &);
  void load(const std::string &);
  void load_chromSizes(const std::string &fname);
};

// immutable IDs based data. IDs cannot be changed after init
class BaseIDsImmutable : public Base {};

// mutable IDs, IDs can be joined or broken
class BaseIDsMutable : public Base {};
};
#endif
