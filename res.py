from core import *
from tools import *
import os
CWD=os.getcwd()

mm9 = GenomeUCSCRefGene(CWD+"/data/mm9.chromInfo.txt",
        CWD+"/data/mm9.refGene.txt")

mm10 = GenomeUCSCRefGene(CWD+"/data/mm10.chromInfo.txt",
        CWD+"/data/mm10.refGene.txt")

hg19 = GenomeUCSCRefGene(CWD+"/data/hg19.chromInfo.txt",
        CWD+"/data/hg19.refGene.txt")

hg38 = GenomeUCSCRefGene(CWD+"/data/hg38.chromInfo.txt",
        CWD+"/data/hg38.refGene.txt")

REFS = {
        "mm9":mm9,
        "mm10":mm10,
        "hg19":hg19,
        "hg38":hg38,
        }
#the sequence db interface works with mafslice
#mm10seq= UCSCGenomeSequence(mm10, "db/mm10_config_test.kch","mm10")

#hg38seq= UCSCGenomeSequence(hg38, "db/hg38_config.kch","hg38")

#SEQS = {"mm10":mm10seq, "hg38":hg38seq}
