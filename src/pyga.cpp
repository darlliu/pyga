#include "ga.hpp"
#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
using namespace boost::python;
using namespace ga;
BOOST_PYTHON_MODULE(pyga) {
  // PyEval_InitThreads();
  class_<std::vector<double>>("DoubleVec")
      .def(vector_indexing_suite<std::vector<double>>());
  class_<std::vector<std::vector<double>>>("DoubleMat")
      .def(vector_indexing_suite<std::vector<std::vector<double>>>());
  class_<std::vector<unsigned>>("IdxVec").def(
      vector_indexing_suite<std::vector<unsigned>>());
  class_<std::vector<std::string>>("StrVec").def(
      vector_indexing_suite<std::vector<std::string>>());
  class_<interval>("IntervalWrapper").def("v", &interval::v);
  class_<gene>("Gene")
      .def("info", &gene::info)
      .def("inv", &gene::inv)
      .def("cds", &gene::cds)
      .def_readonly("id", &gene::id)
      .def_readonly("idx", &gene::idx)
      .def_readonly("sym", &gene::sym)
      .def("utr5", &gene::utr5)
      .def("utr3", &gene::utr3)
      .def("exons", &gene::get_exons)
      .def("introns", &gene::get_introns)
      .def("promoter", &gene::get_promoter)
      .def("tail", &gene::get_tail)
      .def("noncoding", &gene::noncoding);
  class_<std::vector<gene>>("GeneVec").def(
      vector_indexing_suite<std::vector<gene>>());
  class_<std::vector<interval>>("IntervalVec")
      .def(vector_indexing_suite<std::vector<interval>>());
  class_<inv>("Interval")
      .def_readwrite("chr", &inv::chr)
      .def_readwrite("start", &inv::start)
      .def_readwrite("end", &inv::end)
      .def_readwrite("strand", &inv::strand)
      .def_readwrite("seq", &inv::seq);
  class_<UCSCRefGene>("UCSCRefGene")
      .def("iter_genes", &UCSCRefGene::iter_genes)
      .def("iter_genes_unique", &UCSCRefGene::iter_genes_unique)
      .def("get_genes", &UCSCRefGene::get_genes)
      .def("get_gene", &UCSCRefGene::get_gene)
      .def("closest_gene", &UCSCRefGene::closest)
      .def("closest_left", &UCSCRefGene::closest_left)
      .def("closest_right", &UCSCRefGene::closest_right)
      .def("intergene_up", &UCSCRefGene::intergene_up)
      .def("intergene_down", &UCSCRefGene::intergene_down)
      .def("find_sym", &UCSCRefGene::find_sym)
      .def("find_syms", &UCSCRefGene::find_syms)
      .def("find_id", &UCSCRefGene::find_id)
      .def("find_ids", &UCSCRefGene::find_ids)
      .def("size", &UCSCRefGene::size)
      .def("chroms", &UCSCRefGene::chroms)
      .def("r", &UCSCRefGene::r)
      .def("load", &UCSCRefGene::load)
      .def("load_sizes", &UCSCRefGene::load_chromSizes);
};
