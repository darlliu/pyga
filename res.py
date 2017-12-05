from core import *
from tools import *


mm9 = GenomeUCSCRefGene("/baldig/biotools/annotations/UCSC/mm9/chromInfo.txt",
        "/baldig/biotools/annotations/UCSC/mm9/refGene.txt")

mm10 = GenomeUCSCRefGene("/baldig/biotools/annotations/UCSC/mm10/chromInfo.txt",
        "/baldig/biotools/annotations/UCSC/mm10/refGene.txt")

hg19 = GenomeUCSCRefGene("/baldig/biotools/annotations/UCSC/hg19/chromInfo.txt",
        "/baldig/biotools/annotations/UCSC/hg19/refGene.txt")

hg38 = GenomeUCSCRefGene("/baldig/biotools/annotations/UCSC/hg38/chromInfo.txt",
        "/baldig/biotools/annotations/UCSC/hg38/refGene.txt")

REFS = {
        "mm9":mm9,
        "mm10":mm10,
        "hg19":hg19,
        "hg38":hg38,
        }

mm10seq= UCSCGenomeSequence(mm10, "db/mm10_config_test.kch","mm10")

hg38seq= UCSCGenomeSequence(hg38, "db/hg38_config.kch","hg38")

SEQS = {"mm10":mm10seq, "hg38":hg38seq}
