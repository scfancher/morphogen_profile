import math
import numpy as np
import scipy.special as spsp
import matplotlib
import matplotlib.pyplot as plt

def PSDC(j,lam,SP):
	expfact = math.exp(-float(2*j)/lam)
	if SP==1:
		M = expfact*math.sinh(1.0/lam)
		Sig = expfact*(4.0*math.sinh(1.0/lam)-math.exp(-1.0/lam)*(2.0/lam+math.sinh(2.0/lam)))/4.0
		br = 1.0
	elif SP==2:
		thetaN = 50
		uN = 50
		upN = 50
		M = 0.0
		Sig = 0.0
		dsig = 2.0*math.pi/(lam*float(thetaN*uN*upN))
		for i in range(thetaN):
			theta = math.pi*float(2*i+1)/float(thetaN)
			st = math.sin(theta)
			M += (1.0-math.exp(-st/lam)*(1.0+st/lam))/(2.0*st*st)*2.0*math.pi/float(thetaN)
			for k in range(uN):
				u = float(2*k+1)/(2.0*lam*float(uN))
				for l in range(upN):
					up = math.atanh(float(2*l+1)/float(2*upN))
					Sig += u*spsp.jv(0.0,u*up)*spsp.jv(1.0,up/lam)*math.exp(-abs(float(2*j)/lam+u*st))/(1.0+up*up)*dsig/(1.0-float((2*l+1)*(2*l+1))/float(4*upN*upN))
		M *= expfact
		Sig /= 2.0*lam
		br = 2.0
	elif SP==3:
		M = 2.0*math.pi*expfact*(math.cosh(1.0/lam)/lam-math.sinh(1.0/lam))
		Sig = 2.0*math.pi*expfact*((4.0/(lam*lam)+5.0/lam-1.0)*math.exp(-1.0/lam)/8.0+(1.0+1.0/lam)*math.exp(-3.0/lam)/8.0+3.0*math.cosh(1.0/lam)/(4.0*lam)-5.0*math.sinh(1.0/lam)/4.0)
		br = 2.0*math.sqrt(3.0)
	return pow(lam,float(SP-1))*pow(1.0-math.exp(-2.0/lam),2.0)*M*M/(2.0*br*Sig)

def PDT(j,N,lam,phi):
	gam = np.zeros(N,dtype=float)
	if math.isinf(phi):
		if phi > 0.0:
			kap = 1.0/lam
			for i in range(N):
				gam[i] = math.exp(-kap*float(2*(i+1)))
		elif phi < 0.0:
			for i in range(N):
				gam[i] = 1.0
	elif phi == 0.0:
		for i in range(N):
			gam[i] = 1.0/(1.0+float(2*(i+1))/lam)
	else:
		kap = np.sign(phi)*(math.exp(abs(phi))-1.0)*(abs(phi)-math.log(math.exp(abs(phi))-1.0))/lam
		for i in range(N):
			gam[i] = math.exp(-kap*float(2*(i+1)))*(1.0-math.exp(-phi))/(1.0-math.exp(-phi-kap*float(2*(i+1))))
	if j == N:
		DR = gam[j-2]/gam[j-1]-1.0
	else:
		DR = 1.0-gam[j]/gam[j-1]
	return gam[j-1]*DR*DR/(4.0*np.sum(gam))

colmap = matplotlib.cm.get_cmap('bwr')
fignum = 0

lam = np.linspace(7.0,17.0,num=51)
phi = np.logspace(-2.0,0.0,num=51)
rho50 = np.empty((lam.size,3),dtype=float)
PDT50 = np.empty((lam.size,phi.size),dtype=float)
for i in range(lam.size):
	print(i)
	for j in range(phi.size):
		PDTval = PDT(50,100,lam[i],phi[j])
		if j == 0:
			PDTmax = PDTval
		elif PDTval > PDTmax:
			PDTmax = PDTval
		PDT50[i][j] = PDTval
	for j in range(3):
		rho50[i][j] = PDTmax/PSDC(50,lam[i],j+1)
fignum += 1
plt.figure(fignum,figsize=(8.0,6.0))
plt.plot(lam,rho50)
plt.plot(lam,np.ones(lam.size,dtype=float),'k--')
plt.yscale('log')
plt.legend([r'1D',r'2D',r'3D'],loc='upper right',frameon=False)
plt.xlabel(r'Profile Length, $\hat{\lambda}=\lambda/a$')
plt.ylabel(r'$\rho_{j} = P_{DT}^{2}/P_{SDC}^{2}$')
plt.savefig('Fig2B.jpg',dpi=720)
plt.close(fignum)
fignum += 1
plt.figure(fignum,figsize=(8.0,6.0))
for i in range(lam.size):
	plt.plot(phi,PDT50[i],c=colmap(float(lam.size-1-i)/float(lam.size-1)))
plt.ylim([0.0000001,0.000002])
plt.xscale('log')
plt.yscale('log')
plt.xlabel(r'Shape parameter, $\phi$')
plt.ylabel(r'$P_{DT}^{2}/\beta T$')
#ax = plt.gca()
#ax.set_yticks([0.0000004,0.0000005,0.0000006,0.0000007,0.0000008,0.0000009,0.000001,0.000002])
#ax.set_yticklabels([r'$4\times 10^{-7}$','','','','','',r'$10^{-6}$',r'$2\times 10^{-6}$'])
plt.savefig('Fig2A.jpg',dpi=720)
plt.close(fignum)
