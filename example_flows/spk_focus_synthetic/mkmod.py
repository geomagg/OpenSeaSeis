import numpy as np, struct, sys
# modelo 2D: nx colunas x nz amostras (dz), escreve ONL SEISPAK (nbyte 4)
def wonl(fn, V, dz):
    nx, nz = V.shape; ifw0=23
    w=bytearray(4*(ifw0-1))
    for k in range(ifw0-1): w[4*k:4*k+4]=struct.pack('<i',999)
    w[0:4]=struct.pack('<i',ifw0); w[4:8]=struct.pack('<i',0); w[8:12]=struct.pack('<i',nz); w[12:16]=struct.pack('<f',dz)
    w[16:20]=struct.pack('<f',float(nx)); w[20:24]=struct.pack('<f',1.0); w[32:36]=struct.pack('<i',4); w[80:84]=b'END '; w[4*(ifw0-2):4*(ifw0-1)]=struct.pack('<i',ifw0)
    with open(fn,'wb') as fh:
        fh.write(b'ICEp'+struct.pack('<ii',1024,0)); fh.write(w); fh.write(V.astype('<f4').tobytes())
# camadas: (topo [m], velocidade) -- verdadeiro
WB=1000.
LAY=[(0,1500.),(WB,1700.),(1500.,2000.),(2100.,2400.),(2800.,2800.),(3500.,3200.)]
def model(nx, nz, dz, sed=1.0):
    z=np.arange(nz)*dz; v=np.zeros(nz)
    for zt,vl in LAY: v[z>=zt]= vl if zt==0 else vl*sed
    return np.tile(v,(nx,1))
if __name__=='__main__':
    nx=int(sys.argv[1]); nz=int(sys.argv[2]); dz=float(sys.argv[3]); sed=float(sys.argv[4]); fn=sys.argv[5]
    wonl(fn, model(nx,nz,dz,sed), dz)
