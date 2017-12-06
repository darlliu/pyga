from core import *
from db import *

class UCSCGenomeSequence(object):
    """
    Annotation layer for UCSC genome builds
    Requires a motifmapcore.seqdb resource to be set up for
    each genome to be queried
    Please note that the caller is responsible for making sure
    that the core genome annotations match the genome build to be queried
    This annotation layer is special in that it does not annotate actual layers
    but instead returns iterators to sequences
    """
    def __init__(self,genome, dbp, ref):
        self.g=genome
        self.db = SeqDB(dbp,ref)
        return

    def __iter__(self):
        """
        Iterates through the whole gene sequences from tx_start to tx_finish
        """
        for g in self.g:
            inv = g.inv().v()
            inv.seq=self.db.get_seq(inv.chr, inv.start, inv.end)
            yield inv

    def utrs(self):
        for g in self.g:
            inv = g.utr3().v()
            try:
                inv.seq=self.db.get_seq(inv.chr, inv.start, inv.end)
                yield inv
            except:
                print inv.start, inv.end, inv.chr, "Error"
                continue
        for g in self.g:
            inv = g.utr5().v()
            try:
                inv.seq=self.db.get_seq(inv.chr, inv.start, inv.end)
                yield inv
            except:
                print inv.start, inv.end, inv.chr, "Error"
                continue

    def noexons(self):
        """
        iterates through sequences that exclude gene exons
        this include intergene and intron sequences
        """

        for inv in self.g.intergenes():
            if inv.end == 0 or inv.end<= inv.start: continue
            inv.seq= self.db.get_seq(inv.chr, inv.start, inv.end)
            yield inv
        for inv in self.g.introns():
            inv.seq= self.db.get_seq(inv.chr, inv.start, inv.end)
            yield inv

    def intron_junctions(self, span=100):
        """
        A special generator of introns including junction regions
        goes to up to the next exon or the span limit
        """
        for g in self.g:
            for intron in g.introns():
                intron=intron.v()
                start=intron.start-span
                end=intron.end+span
                for exon in g.exons():
                    exon=exon.v()
                    if exon.end==intron.start-1:
                        if exon.start >= intron.start-span:
                            start=exon.start
                    elif exon.start==intron.start+1:
                        if exon.end <= intron.end+span:
                            end=exon.end
                intron.seq=self.db.get_seq(intron.chr, start, end)
                yield intron,start,end

    def introns(self):
        for inv in self.g.introns():
            inv.seq= self.db.get_seq(inv.chr, inv.start, inv.end)
            yield inv

    def exons(self):
        """
        iterates through sequences of gene exons
        """
        for inv in self.g.exons():
            inv.seq= self.db.get_seq(inv.chr, inv.start, inv.end)
            yield inv

    def cds(self):
        """
        iterates through cds sequences
        """
        pass

    def raw(self, chunk=5000):
        """
        iterates through raw sequences in chunk sizes
        """
        for chrom in self.g.p.chroms():
            size = self.g.p.size(chrom)
            start = 0
            end = chunk
            if end >= size:
                inv = pyga.Interval()
                inv.chr=chrom
                inv.start = start
                inv.end=size-1
                inv.strand=True
                inv.seq= self.db.get_seq(chrom, start, size-1)
                yield inv
            while end < size:
                inv = pyga.Interval()
                inv.chr=chrom
                inv.start = start
                inv.end=end
                inv.strand=True
                inv.seq= self.db.get_seq(chrom, start, end)
                yield inv
                start += chunk
                end+= chunk
                if end >=size:
                    inv = pyga.Interval()
                    inv.chr=chrom
                    inv.start = start
                    inv.end=size-1
                    inv.strand=True
                    inv.seq= self.db.get_seq(chrom, start, size-1)
                    yield inv
                    break

    def random_seqs(self, seq_len = 12 , sample_size = 250000,  seed=None, with_annotations=False):
        import numpy as np
        import numpy.random
        import time
        if seed:
            numpy.random.seed(seed)
        else:
            numpy.random.seed(int(time.time()))
        sizes = []
        total = 0
        for data in self.db.db.sizes():
            sizes.append((data.key(), data.data()))
            total += data.data()
        def find_chrom(idx):
            accu = 0
            for chrom, sz in sizes:
                if accu+sz > idx:
                    return chrom,idx-accu
                accu+=sz
            return chrom,idx-accu
        cnt = 0
        for i in numpy.random.randint(0, total-seq_len, sample_size*2):
            chrom, start = find_chrom(i)
            seq = self.db.get_seq(chrom, int(start), int(start+seq_len))
            if "N" in seq:
                continue
            else:
                cnt += 1
                if cnt<=sample_size:
                    if with_annotations:
                        yield seq,chrom,start,start+seq_len
                    else:
                        yield seq
                else: break

    def promoters(self, ups, dns):
        """
        iterates through promoter sequences given upstream and downstream
        """
        for inv in self.g.upstreams(ups, dns):
            inv.seq = self.db.get_seq(inv.chr, inv.start, inv.end)
            yield inv

    def tails(self, ups, dns):
        for inv in self.g.downstreams(ups, dns):
            inv.seq = self.db.get_seq(inv.chr, inv.start, inv.end)
            yield inv

    def background(self, kind="raw"):
        """Calculate ratio of ATCGs"""
        A=0
        T=0
        C=0
        G=0
        if kind == "raw":
            for seq in self.raw():
                seq=seq.seq
                A+=seq.count("A")
                T+=seq.count("T")
                C+=seq.count("C")
                G+=seq.count("G")
        elif kind == "noexon":
            for seq in self.noexons():
                seq=seq.seq
                A+=seq.count("A")
                T+=seq.count("T")
                C+=seq.count("C")
                G+=seq.count("G")
        elif kind == "exons":
            for seq in self.exons():
                seq=seq.seq
                A+=seq.count("A")
                T+=seq.count("T")
                C+=seq.count("C")
                G+=seq.count("G")
        else: raise NotImplemented("Background summarization method not implemented")
        Total = A+C+G+T
        return (float(A)/Total, float(C)/Total, float(G)/Total, float(T)/Total)


def findClosestGeneFromTSS(chrom, start, end, genome=None):
    gene=genome.g.closest_left(chrom,start,end)
    sym=gene.sym
    gene=gene.inv().v()
    gene2=genome.g.closest_right(chrom,start,end)
    sym2=gene2.sym
    gene2=gene2.inv().v()
    dst = int(start + (end-start)/2)
    ds= dst-gene.start
    ds2=dst-gene2.start
    if abs(ds)<=abs(ds2):
        return sym,ds
    else:
        return sym2,ds2

def findClosestGeneFromTES(chrom, start, end, genome=None):
    gene=genome.g.closest_left(chrom,start,end)
    sym=gene.sym
    gene=gene.inv().v()
    gene2=genome.g.closest_right(chrom,start,end)
    sym2=gene2.sym
    gene2=gene2.inv().v()
    dst = int(start + (end-start)/2)
    ds= dst-gene.end
    ds2=dst-gene2.end
    if abs(ds)<=abs(ds2):
        return sym,ds
    else:
        return sym2,ds2

def findClosestGeneFromMiddle(chrom,start, end,genome=None):
    gene=genome.g.closest_left(chrom,start,end)
    sym=gene.sym
    gene=gene.inv().v()
    gene2=genome.g.closest_right(chrom,start,end)
    sym2=gene2.sym
    gene2=gene2.inv().v()
    dst = int(start + (end-start)/2)
    ds= abs(dst-gene.end)
    ds2=abs(gene2.start - dst)
    if abs(ds)<=abs(ds2):
        return sym,ds
    else:
        return sym2,ds2


def overlap(x,y,x2,y2):
    if int(y)>=int(x2) and int(y)<=int(y2):
        return True
    if int(x)>=int(x2) and int(x)<=int(y2):
        return True
    return False
