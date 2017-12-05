#ifndef GA_CPP
#define GA_CPP
#include "ga.hpp"
namespace ga {
std::vector<char *> tokenize2(std::string &in, char *sep) {
  std::vector<char *> out;
  char *pch = strtok(&in[0], sep);
  while (pch != NULL) {
    out.push_back(pch);
    pch = strtok(NULL, sep);
  }
  return out;
};
void tokenize(std::string const &s, char const *d,
              std::vector<string_view> &ret) {
  using namespace std;
  typedef string::const_iterator iter;
  std::vector<string_view> output;
  bitset<255> delims;
  while (*d) {
    unsigned char code = *d++;
    delims[code] = true;
  }
  iter beg;
  bool in_token = false;
  for (string::const_iterator it = s.begin(), end = s.end(); it != end; ++it) {
    if (delims[*it]) {
      if (in_token) {
        output.push_back(
            typename std::vector<string_view>::value_type(beg, it));
        in_token = false;
      }
    } else if (!in_token) {
      beg = it;
      in_token = true;
    }
  }
  if (in_token)
    output.push_back(
        typename std::vector<string_view>::value_type(beg, s.end()));
  output.swap(ret);
}

void UCSCRefGene::load(const std::string &fname) {
  std::ifstream ifs(fname);
  std::string tmp;
  Genes.reserve(50000);
  rlutil::setColor(rlutil::GREEN);
  while (ifs.good()) {
    std::getline(ifs, tmp);
    std::vector<string_view> vals;
    tokenize(tmp, "\t", vals);
    if (vals.size() < 3)
      continue;
    add_entry(vals);
    if (Genes.size() % 1000 == 0)
      std::cerr << ".";
  }
  std::cerr << std::endl;
  rlutil::setColor(rlutil::GREY);
  maxEle = Genes.size();
  return;
}

void UCSCRefGene::load_chromSizes(const std::string &fname) {
  std::ifstream ifs(fname);
  std::string tmp;
  while (ifs.good()) {
    std::getline(ifs, tmp);
    std::vector<string_view> vals;
    tokenize(tmp, "\t", vals);
    if (vals.size() < 3)
      continue;
    sizes[std::string(vals[0].begin(), vals[0].end())] =
        boost::lexical_cast<unsigned long>(vals[1]);
  }
  return;
}

void UCSCRefGene::add_entry(std::vector<string_view> &vals) {
  auto idx = Genes.size();
  gene g;
  g.idx = idx;
  auto rsq = std::string(vals[1].begin(), vals[1].end());
  auto gsm = std::string(vals[12].begin(), vals[12].end());
  (*coreIDs)[rsq] = idx;
  if (GeneSym->count(gsm) == 0)
    (*GeneSym)[gsm] = {idx};
  else
    (*GeneSym)[gsm].push_back(idx);
  g.chr = boost::flyweight<std::string>(
      std::string(vals[2].begin(), vals[2].end()));
  auto sz = sizes[g.chr];
  g.sz = boost::flyweight<unsigned long>(sz);
  g.sym = gsm;
  g.id = rsq;
  char str = boost::lexical_cast<char>(vals[3]);
  if (str == '+')
    g.strand = true;
  else
    g.strand = false;
  g.tx_start = boost::lexical_cast<unsigned>(vals[4]);
  g.tx_end = boost::lexical_cast<unsigned>(vals[5]);
  g.cds_start = boost::lexical_cast<unsigned>(vals[6]);
  g.cds_end = boost::lexical_cast<unsigned>(vals[7]);
  auto cnt = boost::lexical_cast<unsigned>(vals[8]);
  auto s1 = std::string(vals[9].begin(), vals[9].end()),
       s2 = std::string(vals[10].begin(), vals[10].end());
  auto T1 = tokenize2(s1, ",");
  auto T2 = tokenize2(s2, ",");
  for (int i = 0; i < T1.size(); ++i) {
    auto b1 = T1[i];
    auto b2 = T2[i];
    try {
      g.exons.insert(
          interval{g.chr, boost::flyweight<unsigned long>(sizes[g.chr]),
                   boost::lexical_cast<unsigned>(b1),
                   boost::lexical_cast<unsigned>(b2), cnt, g.strand});
    } catch (...) {
      continue;
    }
  }
  GeneLoci.insert(interval{g.chr, boost::flyweight<unsigned long>(sizes[g.chr]),
                           boost::lexical_cast<unsigned>(g.tx_start),
                           boost::lexical_cast<unsigned>(g.tx_end), idx,
                           g.strand});
  Genes.push_back(g);
}

std::vector<unsigned>
BaseGenome::find_syms(const std::vector<std::string> &in) {
  std::vector<unsigned> out;
  out.reserve(in.size());
  for (auto &v : in)
    if (GeneSym->find(v) != GeneSym->end())
      for (auto &vv : GeneSym->at(v))
        out.push_back(vv);
  std::sort(out.begin(), out.end());
  return out;
}

std::vector<unsigned> BaseGenome::find_ids(const std::vector<std::string> &in) {
  std::vector<unsigned> out;
  out.reserve(in.size());
  for (auto &v : in)
    if (coreIDs->find(v) != coreIDs->end())
      out.push_back(coreIDs->at(v));
  std::sort(out.begin(), out.end());
  return out;
}

std::vector<unsigned> BaseGenome::r(const std::string &chr, const int &l,
                                    const int &r) {
  interval test{boost::flyweight<std::string>(chr),
                boost::flyweight<unsigned long>(0),
                l,
                r,
                0,
                true};
  std::vector<unsigned> out;
  for (auto it = GeneLoci.lower_bound(test); it != GeneLoci.upper_bound(test);
       ++it) {
    out.push_back(it->idx);
  }
  std::sort(out.begin(), out.end());
  return out;
}

gene BaseGenome::closest_left(const std::string &chr, const int &l, const int &r) {
  interval test{boost::flyweight<std::string>(chr),
                boost::flyweight<unsigned long>(0),
                l,
                r,
                0,
                true};
  //auto it = GeneLoci.lower_bound(test);
  auto itt = GeneLoci.lower_bound(test);
  itt--;
  if (itt==GeneLoci.end()) return gene();
  return Genes[itt->idx];
}

gene BaseGenome::closest_right(const std::string &chr, const int &l, const int &r) {
  interval test{boost::flyweight<std::string>(chr),
                boost::flyweight<unsigned long>(0),
                l,
                r,
                0,
                true};
  //auto it = GeneLoci.lower_bound(test);
  auto itt = GeneLoci.lower_bound(test);
  if (itt==GeneLoci.end()) return gene();
  return Genes[itt->idx];
}

gene BaseGenome::closest(const std::string &chr, const int &l, const int &r) {
  interval test{boost::flyweight<std::string>(chr),
                boost::flyweight<unsigned long>(0),
                l,
                r,
                0,
                true};
  //auto it = GeneLoci.lower_bound(test);
  auto itt = GeneLoci.lower_bound(test);
  auto it=itt;
  it--;
  if (it == GeneLoci.end() && itt == GeneLoci.end())
    return gene();
  else if (it == GeneLoci.end())
    return Genes[itt->idx];
  else if (itt == GeneLoci.end())
    return Genes[it->idx];
  else if ((it->start <= l)&&(it->end >= r))
    return Genes[it->idx];
  else if ((itt->start <= l)&&(itt->end >= r))
    return Genes[itt->idx];
  if (it->idx ==itt->idx) it--;
  if (std::abs(l - (int)it->end) < std::abs((int)itt->start - r))
    return Genes[it->idx];
  else
    return Genes[itt->idx];
}
gene *const BaseGenome::next_gene(const unsigned &idx) {
  if (idx >= maxEle)
    return nullptr;
  auto g = Genes[idx];
  interval test{boost::flyweight<std::string>(g.chr),
                boost::flyweight<unsigned long>(0),
                g.tx_start,
                g.tx_end,
                0,
                true};
  auto it = GeneLoci.upper_bound(test);
  if (it == GeneLoci.end())
    return nullptr;
  return &Genes[it->idx];
};
gene *const BaseGenome::prev_gene(const unsigned &idx) {
  if (idx >= maxEle)
    return nullptr;
  auto g = Genes[idx];
  interval test{boost::flyweight<std::string>(g.chr),
                boost::flyweight<unsigned long>(0),
                g.tx_start,
                g.tx_end,
                0,
                true};
  auto it = GeneLoci.lower_bound(test);
  if (--it == GeneLoci.end())
    return nullptr;
  return &Genes[it->idx];
};
interval BaseGenome::intergene_up(const unsigned &idx) {
  if (idx >= maxEle)
    return interval();
  auto g = Genes[idx];
  interval out{boost::flyweight<std::string>(g.chr),
               boost::flyweight<unsigned long>(g.sz),
               0,
               g.tx_start - 1,
               g.idx,
               g.strand};
  if (g.tx_start == 0)
    out.end = 0;
  auto gg = prev_gene(idx);
  if (gg == nullptr || (gg->chr != g.chr))
    return out;
  else {
    out.start = gg->tx_end + 1;
    return out;
  }
};
interval BaseGenome::intergene_down(const unsigned &idx) {
  if (idx >= maxEle)
    return interval();
  auto g = Genes[idx];
  interval out{boost::flyweight<std::string>(g.chr),
               boost::flyweight<unsigned long>(g.sz),
               g.tx_end + 1,
               g.sz - 1,
               g.idx,
               g.strand};
  auto gg = next_gene(idx);
  if (gg == nullptr || (gg->chr != g.chr))
    return out;
  else {
    out.end = gg->tx_start - 1;
    return out;
  }
};
};
#endif
