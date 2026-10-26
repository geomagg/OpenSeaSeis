import sys; sys.path.insert(0,'.')
from rd import rd
import numpy as np
from scipy.ndimage import maximum_filter1d, uniform_filter1d
def picks(fn, npt, vmin=1000., vmax=3000., tmin=1500., thr=0.3):
    """picks (T ms, Vavg) por ponto: maximos do painel de coerencia apos AGC em T"""
    d,h,dt=rd(fn); ntr=d.shape[0]//npt; d=d.reshape(npt,ntr,-1); sr=dt/1000.; dv=(vmax-vmin)/(ntr-1)
    out=[]
    for p in range(npt):
        P=np.abs(d[p]); E=P.max(0); E[:int(tmin/sr)]=0
        A=E/(maximum_filter1d(E,int(600/sr))+1e-30)      # AGC 600 ms
        L=[]
        for i in range(2,len(E)-2):
            if A[i]>thr and E[i]==E[max(0,i-int(30/sr)):i+int(30/sr)].max():
                L.append((i*sr, vmin+int(np.argmax(P[:,i]))*dv, A[i]))
        out.append(L)
    return out
def events(P, tol=40.):
    allp=sorted([q for L in P for q in L]); ev=[]
    for q in allp:
        if ev and q[0]-ev[-1][-1][0]<=tol: ev[-1].append(q)
        else: ev.append([q])
    return [(float(np.median([t for t,_,_ in e])), float(np.median([v for _,v,_ in e])), len(e)) for e in ev if len(e)>=len(P)//2]
if __name__=='__main__':
    for t,v,n in events(picks(sys.argv[1],int(sys.argv[2]))): print('T %7.1f ms  Vavg %7.1f m/s  (%d pontos)'%(t,v,n))
