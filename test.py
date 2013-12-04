from milne import milne as milne
import numpy as np
import matplotlib.pyplot as pl

lambda0 = 6301.5080
JUp = 2.0
JLow = 2.0
gUp = 1.5
gLow = 1.833
lambdaStart = 6300.8
lambdaStep = 0.01
nLambda = 100

lineInfo = [lambda0, JUp, JLow, gUp, gLow, lambdaStart, lambdaStep, nLambda]

s = milne(lineInfo)

stokes = np.zeros((4,nLambda))

BField = 100.0
BTheta = 20.0
BChi = 20.0
VMac = 0.0
damping = 0.0
beta = 3.0
mu = 1.0
VDop = 0.1
kl = 5.0
model = [BField, BTheta, BChi, VMac, damping, beta, mu, VDop, kl]

wavelength, stokes = s.synth(model)

fig = pl.figure(num=0)

for i in range(4):
	ax = fig.add_subplot(2,2,i+1)
	ax.plot(wavelength, stokes[i,:])
	
w, d = s.synthDerivatives(model)