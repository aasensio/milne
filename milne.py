import pymilne as pm
import numpy as np

class milne:
	"""
	Class that synthesizes spectral lines using Milne Eddington.
	To use it:
	from milne import milne as milne
	line = milne(lineInfo)
	lineInfo = [lambda0, JUp, JLow, gUp, gLow, lambdaStart, lambdaStep, nLambda]
	wavelength, stokes = line.synth(model)
	model = [BField, theta, chi, vmac, damping, B0, B1, doppler, kl]
	"""
	
	def __init__(self, nLambda, lineInfo):
		self.lineInfo = lineInfo
		self.nLambda = nLambda
		self.wavelength = pm.init(nLambda, lineInfo)
		
	
	def synth(self, model, mu=1.0):
		"""
		Synthesize a spectral line using the Milne-Eddington model
		"""
#		modelNew = self.varChange(model)
		
		stokes = pm.synth(self.nLambda, model, mu)
		
		return stokes
	
	def synthGroup(self, model, mu=1.0):
		"""
		Synthesize many spectral line using the Milne-Eddington model
		"""
		nModels = model.shape[1]		
		stokes = pm.synthGroup(self.nLambda, nModels, model, mu)
		
		return stokes
	
	def __perturbParameter(self, model, index, relativeChange):
		
		newModel = model[:]
		if (model[index] == 0):
			change = 1.0e-3
		else:
			change = model[index] * relativeChange
		
		newModel[index] = model[index] + change
				
		return newModel, change
	
	def varChange(self, p):
		model = p
		model[6] = p[6] - p[5]
		model[8] = (p[8] / p[7])**2		
		
		return model
	
	def varInvChange(self, p):
		model = p
		model[6] = p[5] + p[6]		
		model[8] = p[7] * np.sqrt(p[8])
		
		return model
			
	def synthDerivatives(self, model, mu=1.0, relativeChange=1e-3):
		"""
		Compute the derivative of the Stokes profiles with respect to all the variables
		"""
		
		stokes = self.synth(model, mu)
		
		stokesDeriv = np.zeros((9,4,self.nLambda))
		
		for i in range(9):			
			newModel, change = self.__perturbParameter(model, i, relativeChange)			
			stokesNew = self.synth(newModel, mu)
		
			stokesDeriv[i,:,:] = (stokesNew - stokes) / change			
		
		return stokes, stokesDeriv	
