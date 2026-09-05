"""Odczyt obiektu Seurat .rds BEZ R.

rdata parsuje format RDS poprawnie, ale konwersja pada na asercji w galezi CLO
(funkcje R zapisane w historii polecen Seurata jako SeuratCommand). Te funkcje
nie maja atrybutow i nie sa nam do niczego potrzebne - podmieniamy je na None.
"""
import rdata, numpy as np, scipy.sparse as sp, psutil, sys, time
from rdata.conversion import _conversion as C
GiB=1024**3
_orig=C.SimpleConverter._convert_next
def _patched(self,obj):
    info=getattr(obj,'info',None)
    if info is not None and info.type==rdata.parser.RObjectType.CLO:
        return None                      # funkcje R pomijamy
    return _orig(self,obj)
C.SimpleConverter._convert_next=_patched

def mem(): return psutil.Process().memory_info().rss/GiB
def load(path):
    t=time.time()
    parsed=rdata.parser.parse_file(path)
    print(f"  sparsowane {time.time()-t:.0f}s RAM {mem():.2f} GiB",flush=True)
    t=time.time()
    obj=rdata.conversion.convert(parsed)
    print(f"  skonwertowane {time.time()-t:.0f}s RAM {mem():.2f} GiB",flush=True)
    return obj

def walk_keys(o,pref="",d=0,maxd=2):
    if d>maxd: return
    if isinstance(o,dict):
        print(f"{'  '*d}{pref} -> {list(o.keys())[:18]}")
        for k,v in o.items():
            if k in ("assays","RNA","counts","data","meta.data","cell.embeddings","reductions"):
                walk_keys(v,k,d+1,maxd)
    elif hasattr(o,'__dict__'):
        print(f"{'  '*d}{pref} -> obiekt {type(o).__name__}, pola {list(vars(o))[:14]}")

if __name__=="__main__":
    p=sys.argv[1] if len(sys.argv)>1 else "GSM8315551_SCA7.5wk.seurat.rds.gz"
    print(f"=== {p} ===")
    o=load(p)
    print("typ najwyzszy:",type(o).__name__)
    walk_keys(o,"<root>")
