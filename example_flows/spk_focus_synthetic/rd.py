import numpy as np, struct
def rd(fn):
    b=open(fn,'rb').read(); ns=struct.unpack('>h',b[3220:3222])[0]; dt=struct.unpack('>h',b[3216:3218])[0]
    rec=240+4*ns; n=(len(b)-3600)//rec
    a=np.frombuffer(b[3600:3600+n*rec],'u1').reshape(n,rec)
    d=a[:,240:].copy().view('>f4').astype('f4')
    h=lambda o:a[:,o:o+4].copy().view('>i4').ravel()
    return d,dict(fldr=h(8),trc=h(12),key=h(36),il=h(188),xl=h(192)),dt
