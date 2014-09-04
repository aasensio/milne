from numpy cimport ndarray as ar
from numpy import empty, ascontiguousarray

cdef extern from "pymilne.h":
	void c_setline(int *nLambda, double *lineData, double *waveOut)
	void c_milnesynth(int *nLambda, double *modelIn, double *muIn, double *stokesOut)
	void c_milnesynthmany(int *nLambda, int *nModels, double *modelIn, double *muIn, double *stokesOut)

cdef class milne:
	"""
	Class that synthesizes spectral lines using Milne Eddington.
	To use it:
	from milne import milne as milne
	line = milne(lineInfo)
	lineInfo = [lambda0, JUp, JLow, gUp, gLow, lambdaStart, lambdaStep, nLambda]
	wavelength, stokes = line.synth(model)
	model = [BField, theta, chi, vmac, damping, B0, B1, doppler, kl]
	"""
	
	def __init__(self, int nLambda, ar[double,ndim=1] lineInfo):
		
		cdef:
			ar[double,ndim=1] waveOut = empty(nLambda, order='F')
			
		c_setline(&nLambda, &lineInfo[0], <double*> waveOut.data)
		
		self.lineInfo = lineInfo
		self.nLambda = nLambda		
		self.wavelength = waveOut
		
	def synth(self, ar[double,ndim=1] modelIn, double muIn=1.0):
		
		cdef:
			ar[double,ndim=2] stokesOut = empty((4,nLambda), order='F')
			
		c_milnesynth(&self.nLambda, &modelIn[0], &muIn, <double*> stokesOut.data)
		
		return stokesOut