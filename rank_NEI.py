# Estimation of Noise Variance base on Rank,K., algorithm
#
#

import numpy
import scipy
import math

#Function implements Ranks noise estimation index
def rank_NEI(im):
	K = 1 # 1 gives a 3x3 window
	im = numpy.float64(im)
	#Step 1: Suppression
	#Vertical filter
	vert = numpy.zeros(im.shape,dtype=numpy.float64)
	for i in range(1,im.shape[0] - 1):
		for j in range(1,im.shape[1]-1):
			vert[i,j] = (1.0/math.sqrt(2.0)) * (im[i+1,j] - im[i,j])
	#horizontal filter
	horiz = numpy.zeros(im.shape,dtype=numpy.float64)
	for i in range(1,im.shape[0] - 1):
		for j in range(1,im.shape[1]-1):
			horiz[i,j] = (i/math.sqrt(2)) * (vert[i,j+1] - vert[i,j])
	print "Step 1 (Supression) complete.\n"

	#Step 2: Compute the histogram of local standard  deviations
	L = ((2*K)+1)
	numPix = pow(L,2)
	
	mu = numpy.zeros(im.shape,dtype=numpy.float64)
	sigma = numpy.zeros(im.shape,dtype=numpy.float64)

	for i in range(K+1,im.shape[0] - K):
		for j in range(K+1,im.shape[1]-K):
			tempMu = 0;
			tempSigma = 0;
			#Calculate the local mean 
			winL = horiz[i-K:i+K,j-K:j+K]
			tempMu = winL.mean() #calculates the mean of all the elements i
			#Calculate an estimate of the local noise variance
			for wi in range(1,L-1):
				for wj in range(1,L-1):
					tempSigma = tempSigma + pow((winL[wi,wj] - tempMu),2.0)
			mu[i,j] = tempMu
			sigma[i,j] = tempSigma * (1.0/(numPix-1))

	print "Step 2a complete : Local standard deviation\n"
	sigmaSQ = numpy.sqrt(sigma)
	Hmax = 	numpy.uint8(numpy.ceil(sigmaSQ.max()))
	h = numpy.zeros((1,Hmax))
		
	#create a floating point histogram using rande-based bins
	alpha = 1	
	k = 1
	index = 0

	for i in range(K+1,im.shape[0]-K):
		for j in range(K+1, im.shape[1]-K):
			if(sigmaSQ[i,j] >= 0) and (sigmaSQ[i,j] < 0.5):
				index = index + 1

	h[0][0] = index
	
	for k in range(2,Hmax):
		index = 0
		for i in range(K+1,im.shape[0] - K):
			for j in range(K+1, im.shape[1]-K):
				if (sigmaSQ[i,j] >= (k-0.5)) and (sigmaSQ[i,j] < (k+0.5)):
					index = index + 1
		
		h[0][k] = index

	print "Step 2b complete: (Histogram creation)\n"

	#Step 3: Evaluation of the histogram

	#Calculate the mean square value of the histogram

	sumK = 0.0
	for k in range(1,Hmax):
		sumK = sumK + pow(numpy.float(k),2.0)*h[0][k]

	sumh = numpy.sum(h)
	nS = sumK/sumh

	#calculate an initial global estimate for noise variance
	GnS = nS/(pow(alpha,2))
	print "Step 3 complete: Histogram equaluation\n"

	#Step 4 : Fade out

	gl = numpy.zeros((1,Hmax))
	gl = numpy.float64(gl)
	beta = 2.15
	
	#use a soft-fade out (cosine function to reduce the influence of the original image
	for l in range(1,4): #fade out process is called iteratively (lmax = 4)
		for k in range(1,Hmax):
			sl = math.sqrt(nS) #inital s1
			if(k <= sl):
				gl[0][k] = 1.0
			elif (k > sl and k <= (beta*sl)):
				T = ((beta-numpy.float(k)/sl)/(1.0-beta)) * math.pi
				T1 = 1.0 - numpy.cos(T)
				T2 = 0.5 *T1
				gl[0][k] = T2
			elif (k > (beta*sl)):
				gl[0][k] = 0.0

		#calculate an improved value of the mean square

		sumK = 0.0
		sumG = 0.0

		for k in range(1,Hmax):
			sumK = sumK + pow(numpy.float(k),2) * gl[0][k] *h[0][k]
			sumG = sumG + gl[0][k] * h[0][k]
		
		nS = sumK/sumG
		print "sumK = ",sumK , "sumG = ",sumG,"nS = ",nS	
	print "Step 4 complete: Soft fade-out\n"
	rank = nS/pow(alpha,2)
	return rank
