import pandas as pd, numpy as np, os, re, ast, glob
ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", ".."))  # repository root
REPO=ROOT + "/data/transient-dynamics_output-current/example_scenarios/standard_model"
S=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
names=sorted(os.listdir(f"{S}/runs"))
for n in names:
    ref=pd.read_excel(f"{REPO}/{n}/results_standard/array_outputs.xlsx")
    new=pd.read_excel(f"{S}/runs/{n}/generated-data/standard_model/results_standard/array_outputs.xlsx")
    print(f"\n===== {n}  shape ref{ref.shape} new{new.shape}  cols equal={list(ref.columns)==list(new.columns)}")
    worst=[]
    for col in ref.columns:
        a,b=ref[col],new[col]
        if pd.api.types.is_numeric_dtype(a) and pd.api.types.is_numeric_dtype(b):
            A,B=a.to_numpy(float),b.to_numpy(float)
            nanmis=int((np.isnan(A)!=np.isnan(B)).sum())
            m=~(np.isnan(A)|np.isnan(B))
            d=np.abs(A[m]-B[m]); rel=d/np.maximum(np.abs(A[m]),1e-300)
            exact=np.array_equal(A,B,equal_nan=True)
            worst.append((col,'num',exact,nanmis,d.max() if d.size else 0, np.nanmax(np.where(np.abs(A[m])>1e-12,rel,0)) if d.size else 0))
        else:
            ne=int((a.astype(str)!=b.astype(str)).sum())
            worst.append((col,'str',ne==0,ne,None,None))
    for w in worst:
        if not w[2]: print("  DIFF", w)
    print("  identical columns:", sum(w[2] for w in worst), "/", len(worst))
    def parse(path):
        line=open(path).read().split('\n')[1]
        lists=re.findall(r'\[[^\[\]]*\]', line.split('{',1)[0])
        stats=eval(line[line.index('{'):], {'nan':float('nan')})
        return lists, stats
    rl,rs=parse(f"{REPO}/{n}/results_standard/non_array_outputs.txt")
    nl,ns=parse(f"{S}/runs/{n}/generated-data/standard_model/results_standard/non_array_outputs.txt")
    print("  tipping ref:", rl); print("  tipping new:", nl)
    md=0;mism=0
    for reg in rs:
        for k in rs[reg]:
            for s in rs[reg][k]:
                x,y=rs[reg][k][s],ns.get(reg,{}).get(k,{}).get(s,np.nan)
                if np.isnan(x) and np.isnan(y): continue
                if np.isnan(x)!=np.isnan(y): mism+=1; continue
                md=max(md,abs(x-y))
    print(f"  regime stats: regimes ref={len(rs)} new={len(ns)}  nan-mismatch={mism}  max abs diff={md:.3g}")
    rf=sorted(os.path.basename(p) for p in glob.glob(f"{REPO}/{n}/figures_standard/*"))
    nf=sorted(os.path.basename(p) for p in glob.glob(f"{S}/runs/{n}/generated-data/standard_model/figures_standard/*"))
    print("  figures only in ref:", sorted(set(rf)-set(nf)), " only in new:", sorted(set(nf)-set(rf)))
