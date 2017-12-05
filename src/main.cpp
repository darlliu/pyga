#include "ga.cpp"
#include <google/profiler.h>
using namespace ga;
int main(int argc, char **argv) {
  ProfilerStart("/home/yul13/pprof/libgatest.prof");
  auto gn = UCSCRefGene();
  gn.load_chromSizes("/baldig/biotools/annotations/UCSC/mm10/chromInfo.txt");
  gn.load("/baldig/biotools/annotations/UCSC/mm10/refGene.txt");
};
