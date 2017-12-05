#include "ga.cpp"
using namespace ga;
int main(int argc, char **argv) {
  auto gn = UCSCRefGene();
  gn.load_chromSizes("/baldig/biotools/annotations/UCSC/mm10/chromInfo.txt");
  gn.load("/baldig/biotools/annotations/UCSC/mm10/refGene.txt");
};
