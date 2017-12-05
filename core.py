# Core functionality of the library
import bin.pyga as pyga
import json
import cPickle as pkl

class GenomeCore(object):

    def __init__(self, p=None,unique=False):
        """
        Create a genome core object
        p is a pyga.genome object
        Unique argument is used for getting non-overlapping genes
        """
        self.p=p;
        self.order_ = None
        self.unique =unique

    def load(self, fname_sz, fname):
        raise NotImplemented("Please Implement load routine!")

    def __iter__(self):
        """
        This routine iterates through genes in chromosomal order
        """
        if self.unique:
            iter_method= self.p.iter_genes_unique
        else:
            iter_method= self.p.iter_genes
        for idx in iter_method():
            if self.order_ is None:
                yield self.p.get_gene(idx)
            else:
                if idx in self.order_:
                    yield self.p.get_gene(idx)
                else: continue
    @property
    def g(self):
        """for legacy code"""
        return self
    def clear(self):
        self.order_ = None

    def find_sym(self, sym):
        try:
            for idx in self.p.find_sym(sym):
                yield self.p.get_gene(idx)
        except:
            print "Gene not found:", sym
            yield None

    def find_id(self, sym):
        try:
            idx =  self.p.find_id(sym)
            yield self.p.get_gene(idx)
        except:
            print "ID not found",sym
            yield None

    def filter_syms(self, syms):
        self.order_=[]
        SYMS = pyga.StrVec()
        SYMS[:] = syms
        for idx in self.p.find_syms(SYMS):
            self.order_.append(idx)
        return

    def filter_ids(self,ids):
        self.order_=[]
        SYMS = pyga.StrVec()
        SYMS[:] = ids
        for idx in self.p.find_ids(SYMS):
            self.order_.append(idx)
        return

    def find_syms(self,syms):
        self.filter_syms(syms)
        for i in self:
            yield i

    def find_ids(self, ids):
        self.filter_ids(ids)
        for i in self:
            yield i

    def intergenes(self):
        if self.unique:
            idxs_ = self.p.iter_genes_unique()
        else:
            idxs_ = self.p.iter_genes()
        if self.order_ is None:
            idxs = idxs_
        else:
            idxs=[]
            for idx in idxs_:
                if idx in self.order_:
                    idxs.append(idx)
        chrs = []
        for idx in idxs:
            inv = self.p.intergene_down(idx).v()
            if inv.start == inv.end: continue
            if inv.chr not in chrs:
                chrs.append(inv.chr)
                inv2 = self.p.intergene_up(idx).v()
                if inv2.start == inv2.end: continue
                yield inv2
            yield inv

    def utr3s(self):
        for g in self:
            yield g.utr3().v()

    def utr5s(self):
        for g in self:
            yield g.utr5().v()

    def introns(self):
        for g in self:
            for i in g.introns():
                yield i.v()

    def exons(self):
        for g in self:
            for i in g.exons():
                yield i.v()

    def upstreams(self, ups = 10, dns=10):
        for g in self:
            yield g.promoter(abs(ups), abs(dns)).v()

    def downstreams(self, ups = 10, dns=10):
        for g in self:
            yield g.tail(abs(ups), abs(dns)).v()

    def closest(self, chrom, l, r):
        return self.p.closest_gene(chrom,l,r)

    def closest_left(self, chrom, l, r):
        return self.p.closest_left(chrom,l,r)

    def closest_right(self, chrom, l, r):
        return self.p.closest_right(chrom,l,r)

    def r(self, chrom, l, r):
        for idx in self.p.r(chrom, l, r):
            yield self.p.get_gene(idx)

    def __len__(self):
        """
        returns the length
        """
        return self.p.count()

    def __setitem__(self, key, value):
        raise ValueError("Genome object does not allow changing values")

    def __getitem__(self, key):
        if (type(key)!=int):
            raise TypeError("Cannot allow non integer key")
        return self.p.get_gene(key)

    def __delitem__ (self, key):
        raise ValueError("Genome object does not allow deletion")

    def __getslice__ (self, i,j):
        out = []
        for idx in xrange(i,j):
            g = self.p.get_gene(idx)
            if g.tx_start != g.tx_end:
                out.append()
        return out

    def __setslice__ (self, i,j, key):
        raise ValueError("Genome object does not allow changing values")

    def __delslice__ (self, i,j):
        raise ValueError("Genome object does not allow deletion")

class GenomeUCSCRefGene(GenomeCore):

    def __init__(self, sz, rf):
        p=pyga.UCSCRefGene()
        super(GenomeUCSCRefGene, self).__init__(p)
        self.load(sz,rf)
        return

    def load(self, sz, rf):
        print type (self.p)
        self.p.load_sizes(sz)
        self.p.load(rf)
        return

