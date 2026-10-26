import sys; sys.path.insert(0,'.')
from rd import rd
import numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from scipy.ndimage import maximum_filter1d
plt.rcParams.update({'font.size':10})
Q=4  # ponto de analise DPN 260 (x = 6,5 km)
def agc(P, w=300):
    E=np.abs(P).max(0); return P/(maximum_filter1d(E,w)+1e-30)
def fig(p, title, fn):
    m,_,_=rd(p+'_mig.sgy'); t,h,_=rd(p+'_tau.sgy'); t=t.reshape(9,126,-1)[Q]; c,_,_=rd(p+'_vfix.sgy'); c=c.reshape(9,126,-1)[Q]; va,_,_=rd(p+'_vavg.sgy'); va=va[Q]
    T=np.arange(t.shape[1])*0.002
    f,ax=plt.subplots(1,3,figsize=(13,5.2),gridspec_kw=dict(width_ratios=[1.6,1,1]))
    cl=np.percentile(np.abs(m),99)
    ax[0].imshow(m.T,aspect='auto',cmap='gray',vmin=-cl,vmax=cl,extent=[0,12.8,4.0,0]); ax[0].set_xlim(2,10.8)
    for z in (1.5,2.1,2.8,3.5): ax[0].axhline(z,color='tab:red',lw=0.6,ls='--')
    ax[0].axvline(6.5,color='tab:blue',lw=0.8); ax[0].set_xlabel('x (km)'); ax[0].set_ylabel('profundidade (km)'); ax[0].set_title('Imagem migrada (vermelho: refletores verdadeiros)',fontsize=10)
    A=agc(t); ax[1].imshow(np.abs(A).T,aspect='auto',cmap='jet',vmin=0,vmax=1,extent=[-500,500,4.0,0]); ax[1].axvline(0,color='w',lw=0.8)
    ax[1].set_ylim(4,1.0); ax[1].set_xlabel('desvio tau (ms)'); ax[1].set_ylabel('T (s)'); ax[1].set_title('Painel tau, DPN 260',fontsize=10)
    B=agc(c); ax[2].imshow(np.abs(B).T,aspect='auto',cmap='jet',vmin=0,vmax=1,extent=[1000,3000,4.0,0]); ax[2].plot(va,T,'w',lw=1.2)
    ax[2].set_ylim(4,1.0); ax[2].set_xlim(1200,2600); ax[2].set_xlabel('Vavg (m/s)'); ax[2].set_ylabel('T (s)'); ax[2].set_title('Coerência (SPK_TAU2VFIX) + VAVG da migração',fontsize=10)
    f.suptitle(title,fontsize=12,fontweight='bold'); plt.tight_layout(); plt.savefig(fn,dpi=110); plt.close()
fig('TRUE','Passo 1 — migração com a velocidade verdadeira','fig_true.png')
fig('V90','Passo 2 — migração com 90% da velocidade nos sedimentos','fig_v90.png')
fig('IT2','Passo 4 — migração com o modelo da segunda atualização','fig_it2.png')
# registro sintetico


d,h,dt=rd('node9.sgy'); c=np.percentile(np.abs(d),99.5)
f,ax=plt.subplots(1,2,figsize=(12,5),gridspec_kw=dict(width_ratios=[1.3,1]))
ax[0].imshow(d.T,aspect='auto',cmap='gray',vmin=-c,vmax=c,extent=[0,12.8,4,0]); ax[0].set_xlabel('posicao do tiro x (km)'); ax[0].set_ylabel('t (s)'); ax[0].set_title('Registro do no 9 (x = 6,4 km)',fontsize=10)
z=np.arange(0,4.5,0.005)
def prof(s):
    L=[(0,1500.),(1.0,1700.),(1.5,2000.),(2.1,2400.),(2.8,2800.),(3.5,3200.)]; v=np.zeros_like(z)
    for zt,vl in L: v[z>=zt]= vl if zt==0 else vl*s
    return v
ax[1].plot(prof(1.0),z,'k',label='verdadeiro'); ax[1].plot(prof(0.9),z,'tab:red',ls='--',label='90% nos sedimentos')
ax[1].set_ylim(4.2,0); ax[1].set_xlabel('V (m/s)'); ax[1].set_ylabel('profundidade (km)'); ax[1].legend(loc='lower left'); ax[1].set_title('Modelos de velocidade',fontsize=10)
plt.tight_layout(); plt.savefig('fig_dados.png',dpi=110)
