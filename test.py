from milne import milne as milne
import numpy as np
import matplotlib.pyplot as pl
import datetime as dt

lambda0 = 6301.5080
JUp = 2.0
JLow = 2.0
gUp = 1.5
gLow = 1.833
lambdaStart = 6300.8
lambdaStep = 0.03
nLambda = 50

lineInfo = np.asarray([lambda0, JUp, JLow, gUp, gLow, lambdaStart, lambdaStep])

s = milne(nLambda, lineInfo)


stokes = np.zeros((4,nLambda))

BField = 100.0
BTheta = 20.0
BChi = 20.0
VMac = 2.0
damping = 0.0
B0 = 0.8
B1 = 0.2
mu = 1.0
VDop = 0.085
kl = 5.0
model = np.asarray([BField, BTheta, BChi, VMac, damping, B0, B1, VDop, kl])

start = dt.datetime.now()
for i in range(10000):
	stokes = s.synth(model,mu)
end = dt.datetime.now()
print "Computing 1000 models without derivatives took {0} s".format((end-start).microseconds*1e-6)

model = np.tile(model, (10000,1)).T
start = dt.datetime.now()
stokes = s.synthGroup(model,mu)
end = dt.datetime.now()
print "Computing 1000 models without derivatives took {0} s".format((end-start).microseconds*1e-6)


#start = dt.datetime.now()
#for i in range(1000):
	#stokes, stokesDer = s.synthDerivatives(model,mu)
#end = dt.datetime.now()
#print "Computing 1000 models with derivatives took {0} s".format((end-start).microseconds*1e-6)

#fig = pl.figure(num=0)

#for i in range(4):
	#ax = fig.add_subplot(2,2,i+1)
	#ax.plot(wavelength, stokes[i,:])
	
#fig.tight_layout()

#fig2 = pl.figure(num=1, figsize=(16,8))
#loop = 1
#for i in range(4):
	#for j in range(9):
		#ax = fig2.add_subplot(4,9,loop)
		#ax.plot(wavelength, stokesDer[i,j,:])
		#loop += 1
#fig2.tight_layout()
