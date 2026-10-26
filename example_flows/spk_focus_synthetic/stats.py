import sys; sys.path.insert(0,'.')
from rd import rd
import numpy as np, json
from scipy.signal import hilbert
TT=[1922,2522,3105,3605]; ZT=[1500,2100,2800,3500]
def dtau(p):
    d,h,dt=rd(p+'_tau.sgy'); d=d.reshape(9,126,-1); key=np.array(h['key'][:126])
    R=[]
    for T in TT:
        vals=[]
        for q in range(9):
            P=np.abs(d[q]); i0=int(T/2); w=P[:,i0-60:i0+60]; j,k=np.unravel_index(np.argmax(w),w.shape); vals.append((key[j],(i0-60+k)*2))
        R.append((float(np.median([a for a,_ in vals])), float(np.median([b for _,b in vals]))))
    return R
def depths(p):
    d,h,dt=rd(p+'_mig.sgy'); tr=d[180:340].mean(0); e=np.abs(hilbert(tr)); out=[]
    for z in ZT:
        i0=int(z/5); w=e[i0-60:i0+60]; out.append((i0-60+int(np.argmax(w)))*5)
    return out
res={}
for p in ('V90','IT1','IT2','IT3','TRUE'):
    res[p]=dict(dt=dtau(p),z=depths(p)); print(p,'DT',[int(a) for a,_ in res[p]['dt']],'z',res[p]['z'])
json.dump(res,open('stats.json','w'))
