# Description: Image Noise Removal main 

import os
import glob
import time
import sys
import numpy

from SimpleCV import *
from imageIO import *
from rank_NEI import *
from imenh_lib import *
from imfilter_lib import *

#Function calculates Rank's noise estimation index
def calcRankMetric(image):
	#Calculate ranks noise estimation index
	NEI_metric = rank_NEI(image)
	return NEI_metric

#Function displays the histogram of an image
def dispHisto(image):
	plot_IMGhist(image,256)

#Function that applies a noise suppression algorithm "mode" filtering
def mode_filter(image):
	result_img = enh_truncMedian(image,5)
    	
	return result_img   

#Function that applies the hybrid median filter
def hybrid_filter(image):
	result_img = enh_hybridMedian(image,5)
	return result_img

#Function that applies alpha trimmed means
def alpha_filter(image,alpha):
    alpha = float(alpha)
	result_img = enh_alphaTMean(image,alpha,5)
	return result_img

#Function performs Gaussian filtering 
def gaus_filter(image,sigma,width):
	result_img = img_filterGAUSSIAN(image,sigma,width)
	return result_img

#Function to save results to file
def saveFile(image,arg,algo):
	result_path = os.getcwd()+"/result_"+algo+"_"+arg
	imwrite_gray(result_path,image)

#Function reads gray scale image and returns 8bit 2d array
def convertGray(image_path):
	img = imread_gray(image_path)
	return img

#Checks command line arguments to make sure an image was specified
if len(sys.argv) < 2:
	print "Format: >filename pic.ext [filter_type]"
	print "Must specify an image file with .ext"
else:
	img_ext = sys.argv[1].find(".")
	filter_type = sys.argv[2]
	if(filter == "alpha"):
		if(sys.argv[3] == None):
			print "Must enter an integer number\nFormat: pic.ext alpha [1.41]"
			sys.exit()
	if(filter_type == "gaussian"):
		if(sys.argv[3] == None or sys.argv[4] == None):
			print "Must enter an integer number\nFormat: pic.ext gaussian [1.31] [2]"		
			sys.exit()

	if(img_ext > 0 ):
		path = os.getcwd() #gets current directory
		img_path = path + '/'+ sys.argv[1] #appends image file to current directory
		img = convertGray(img_path) #reads gray scale image and returns 8bit 2d array
		before_rank = calcRankMetric(img)
		if(filter_type == "mode"):
			#method performs mode filtering
			print "\nPerforming mode filtering"
			mode_img = mode_filter(img)
			print "Rank before",before_rank
			after_rank = calcRankMetric(mode_img)
			print "Rank after",after_rank
			saveFile(mode_img,sys.argv[1],"mode")
		elif(filter_type == "hybrid"):
			print "\nPerforming hybrid filtering"
			#method performs hybrid filtering best with (Rayleigh images)
			hybrid_img = hybrid_filter(img)
			print "Rank before",before_rank
			after_rank = calcRankMetric(hybrid_img)
			print "Rank after",after_rank
			saveFile(hybrid_img,sys.argv[1],"hybrid")
		elif(filter_type == "alpha"):
			if(sys.argv[3] == None):
				print "Must enter an alpha number"
			else:
				print "\nPerforming alpha trim filtering"
				#method performs alpha trimmed filtering (best with Gaussian images)
				alpha_img = alpha_filter(img,sys.argv[3])
				#print "Rank before",before_rank
				after_rank = calcRankMetric(alpha_img)
				print "Rank after",after_rank
				path = sys.argv[3]+sys.argv[1]
				saveFile(alpha_img,path,"alpha")
		elif(filter_type == "gaussian"):
			print "\nPerforming gaussian filtering"
			if(sys.argv[3] == None and sys.argv[4] == None):
				print "Must enter an integer and a whole number"
			else:		
				#method performs Gaussian filtering (image smoothing)
				gaus_im = gaus_filter(img,sys.argv[3],sys.argv[4])
				print "Rank before",before_rank
				after_rank = calcRankMetric(gaus_im)
				print "Rank after",after_rank				
				path = sys.argv[4]+sys.argv[1]
				saveFile(gaus_im,path,"gaus")
		else:
			print "Incorrect argument for filter_type:\nOptions: alpha\nmode\nhybrid\ngaussian"

		dispHisto(img)

		
    	else:
		print "Invalid extension: example(.tif, .png, .jpeg)"


