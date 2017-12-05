# Implements db interfaces
# Mostly contains code that are used internally only

import pymysql as MySQLdb
import pandas as pd
class SeqDB(object):
    """
    Motifmap SeqDB interface
    """
    def __init__(self, dbp, ref):
        import motifmap.motifmapcore as motifmapcore
        import numpy as np
        print """Please note that this module requires a correct \
working directory as MotifMap uses relative paths only.\
If you see C++ exception check your paths."""
        self.db = motifmapcore.seqdb()
        self.db.load(dbp, ref)
        return

    def get_seq(self, chrom, start, end):
        if type (chrom)!= str or type(start)!=int or type(end)!=int:
            raise TypeError("Query value type error")
        return self.db.get(chrom, start,end)

    def moods_score(self,seq,mat, bg, thr):
        if len (mat[0])==4:
            mat = list(np.array(mat).transpose())
        MATRIX = motifmapcore.DoubleMat()
        for row in mat:
            d=motifmapcore.DoubleVec()
            for ele in row:
                d.append(ele)
            MATRIX.append(d)
        BACKGROUND = motifmapcore.DoubleVec()
        for ele in bg:
            BACKGROUND.append(ele)
        scores = self.db.moods_scan(seq, MATRIX, BACKGROUND, thr)
        return scores

    def moods_score_naive(self, seq, mat, thr):
        if len (mat[0])==4:
            mat = list(np.array(mat).transpose())
        MATRIX = motifmapcore.DoubleMat()
        for row in mat:
            d=motifmapcore.DoubleVec()
            for ele in row:
                d.append(ele)
            MATRIX.append(d)
        BACKGROUND = motifmapcore.DoubleVec()
        for ele in bg:
            BACKGROUND.append(ele)
        scores = self.db.moods_naive_scan(seq, MATRIX, thr)
        return scores


class ListDB(object):
    """
    A simple interface to delimiter based lists
    """
    def __init__(self, fname, tname=None, kind="simple", sep="\n", col=None, delim="\t"):
        """
        1, delimiter separated lists. Only contains ids, delimiters and white spaces
        2, multiple columns db with no header, specify which column to take as id
        3, a table with a header, specify header name
        """
        if type(fname)==list or type(fname)==set:
            self.data = sorted(list(fname))
            self.data_unique=set(fname)
            if not tname:
                self.tname="CustomSet"
            else:
                self.tname=tname
            return
        if not tname:
            self.tname=fname
        else:
            self.tname=tname
        if kind not in ["simple","multiple","table"]:
            raise NotImplemented("This kind of list is not supported")
        if kind == "multiple":
            if type(col)!=int:
                raise TypeError("Please provide correct col number")
            lines = [i.strip() for i in open(fname).read().split(sep) if i.strip()!=""]
            self.data = [i.split(delim)[col] if len(i.split(delim))> col else "" for i in lines]
            self.data = [i.strip() for i in self.data if i.strip()!=""]
        elif kind == "table":
            if type(col)!=str:
                raise TypeError("Please provide correct col name")
            import pandas as pd
            df = pd.read_csv(fname, sep=sep)
            self.data = df[col].values
        else:
            self.data = [i.strip() for i in open(fname).read().split(sep) if i.strip()!=""]
        print "Loaded a list with ", len(self.data),"elements"
        self.data_unique=set(self.data)
        return

    def cross(self, entry):
        if entry in self.data_unique:
            return True
        else:
            return False
    def cross_mult (self, list_entry):
        return sorted(list(self.data_unique& set(list_entry)))

class TableDB(object):
    """
    A database to query csv like tables
    """
    def __init__(self, fname, sep="\t", index=None):
        import pandas as pd
        if index:
            self.df = pd.read_csv(fname, sep=sep, index= index)
        else:
            self.df = pd.read_csv(fname, sep=sep)
    def cross_mult(self, col, list_in):
        return self.df[self.df[col].isin(list_in)]

class MySQLDB(object):
    """
    Interface to MySQL databases
    """
    def __init__(self,db,
            MYSQL_HOST="",
            MYSQL_PORT=,
            MYSQL_USER="",
            MYSQL_PASS=""
            ):
        print "TRYING TO CONNECT TO DB",db
        self.dbname=db
        self.host = MYSQL_HOST
        self.port=MYSQL_PORT
        self.user=MYSQL_USER
        self.pwd=MYSQL_PASS
        self.db = MySQLdb.connect(host=self.host,
                                  user=self.user,
                                  passwd=self.pwd,
                                  db=db,
                                  port=self.port)
        self.cur = self.db.cursor()

    def connect(self):
        self.db = MySQLdb.connect(host=self.host,
                                  user=self.user,
                                  passwd=self.pwd,
                                  db=self.dbname,
                                  port=self.port)
    def query(self, sql):
        try:
            cursor = self.db.cursor()
            cursor.execute(sql)
        except (AttributeError, MySQLdb.OperationalError):
            self.connect()
            cursor = self.db.cursor()
            cursor.execute(sql)
        return cursor
    def execute_sql(self,sql):
        """
        Execute SQL with the opened dbd
        """
        self.query(sql)

    def fetch_sql(self,sql):
        """
        Fetchall instead of commit changes
        """
        self.query(sql)
        return self.cur.fetchall()

    def pandas_sql(self,sql):
        print "EXECUTING SQL", sql
        df= pd.read_sql(sql, self.db)
        if len(df)==0:
            print "WARNING: pandas table has 0 rows"
        return df

class MongoDB(object):
    """
    Interface to pymongo databases
    """
    def __init__(self,db,
            MONGO_SERVER = '',
            MONGO_PORT = 0 
            ):
        import pymongo
        self.con = pymongo.MongoClient(MONGO_SERVER, MONGO_PORT)
        self.db = self.con[db]

    def find_one_from_table(self,tbl, key ="_id", val=None):
        if val is None:
            raise ValueError("Please provide query value")
        tbl = self.db[tbl]
        return tbl.find_one({key:str(val)})

    def find_multi_from_table(self, tbl, key=None, val=None):
        if key is None:
            return self.db[tbl].find()
        else:
            if val is None:
                raise ValueError("Please provide query value")
            tbl = self.db[tbl]
            return tbl.find({key:str(val)})
    def drop_table(self,tbl):
        print "DROPPING TABLE {}".format(tbl)
        self.db[tbl].drop()

    def insert(self, tbl, query={}):
        print "INSERTING TO {}".format(tbl)
        self.db[tbl].insert(query)
        return



class IDMapperDB(object):
    """
    A db layer for id conversion via idmapper.py
    """
    def __init__(self,
        PORTNUM=0,
        HOST=''
            ):
        import os, sys, socket
        self.port=PORTNUM
        self.host=HOST

    def convertIDs (self, ids=None, fromType="OFFICIAL_GENE_SYMBOL", toType="UNIPROT_ACCESSION", species="",
            asType="pairs", filterNA=True):
        """
        This client code uses the ID mapper server to convert IDs, currently it is mainly useful for
        converting various IDs to and from Uniprots
        the fromType and toType listings are included above and species only work if the to type is
        uniprot and it filters the converted uniprots to validated species specific prots only.
        """
        HOST=self.host
        PORTNUM=self.port
        FROMTYPES=[
                "OFFICIAL_GENE_SYMBOL",
                "ENTREZ_GENE_ID",
                "UNIPROT_ACCESSION",
                "UCSC_GENE_ID",
                "REFSEQ_MRNA",
                "REFSEQ_PROTEIN",
                "REFSEQ_NT2",
                "REFSEQ_NT",
                ]
        TOTYPES=[
                "OFFICIAL_GENE_SYMBOL",
                "ENTREZ_GENE_ID",
                "UNIPROT_ACCESSION",
                "UNIPROT_ACCESSION_SPECIES",
                "GO_ACCESSION",
                "UCSC_GENE_ID",
                "REFSEQ_MRNA",
                "REFSEQ_PROTEIN",
                "REFSEQ_NT2",
                "REFSEQ_NT",
                ]
        SPECIES=[
                "MOUSE",
                "HUMAN"]
        if ids is None:
            print "Please input list of ids or text of ids"
            return None
        if fromType not in FROMTYPES:
            print "from type not available, available ones are ", FROMTYPES
            return None
        if toType not in TOTYPES:
            print "to type not available, available ones are ", TOTYPES
            return None
        if species is not "":
            if species not in SPECIES or toType is not "UNIPROT_ACCESSION":
                print "Species can only be HUMAN or MOUSE and only when to type is UNIPROT_ACCESSION"
                return None
        if type(ids)==str:
            ids=[i.strip() for i in ids.split("\n") if i.strip()!=""]
            ids=[i.split("\t") for i in ids]
            ids = [i for j in ids for i in j if i!=""]
            #grabs formatted tsv
        if type (ids)==list:
            ids=[str(i) for i in ids if str(i)!=""]
        if species:
            cmd = "from "+ fromType + " to " + toType+"_"+species
        else:
            cmd = "from "+ fromType + " to " + toType
        try:
            results = []
            for _ids in [ids [i:i+100] for i in range(0, len(ids), 100)]:
                msg = cmd+ "\n" + "\n".join(_ids)+"\nFIN"
                sock=socket.socket(socket.AF_INET,socket.SOCK_STREAM)
                sock.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
                sock.connect((HOST,PORTNUM))
                sock.sendall(msg)
                buf=""
                while True:
                    tbuf=sock.recv(4096)
                    if not tbuf:
                        break
                    buf+=tbuf
                _results = [(i.split(",")[0],i.split(",")[1]) for i in buf.split("\n") if "," in i]
                results += sorted(list(set(_results)))
            if filterNA:
                results=filter(lambda x: x[1]!="N/A", results)
            for i in ids:
                if i not in dict(results):
                    results+=[(i,"NA")]
            if asType=="pairs":
                pass
            elif asType == "tabs":
                results = [i.split(",")[0]+"\t"+i.split(",")[1] for i in results]
                results = "\n".join(results)
            return results
        except socket.error as err:
            print "Error occured: " + repr(err)
            return None
