# Library to process images using spatial filtering
#
# img_filter         - Image filtering with a 3x3 kernel
# img_filterMEAN     - Image smoothing with a 3x3 neighbourhood
# img_filterGAUSSIAN - Image smoothing with Gaussian filter
# img_sharpUSM       - Image sharpening with unsharp masking
# img_sharpUSMgauss  - Image sharpening using Gaussian
# img_sharpCUSM      - Image sharpening with Cubic UM

import numpy
import math
import scipy.ndimage as nd
import pylab

# Function to filter an image using a 3x3 kernel
#
def img_filter(im, kernel):
    img = numpy.zeros(im.shape,dtype=numpy.int16)
    
    # Process every pixel in the image
    for r in range(1,im.shape[0]-1):
        for c in range(1,im.shape[1]-1):
            sumP = 0
            for i in range(3):
                for j in range(3):
                    pixel = im[r+i-1][c+j-1]
                    sumP = sumP + pixel * kernel[i][j]
            # Suppress over/undershooting pixels
            if sumP > 255:
                img[r][c] = 255
            elif sumP < 0:
                img[r][c] = 0
            else:
                img[r][c] = sumP
        print "."
    return img
    
# Function to perform image smoothing using a 3x3 neighborhood mean
#
def img_filterMEAN(im):
    img = numpy.zeros(im.shape,dtype=numpy.int16)
    
    for r in range(1,im.shape[0]-1):
        for c in range(1,im.shape[1]-1):
            sumP = 0
            for i in range(3):
                for j in range(3):
                    pixel = im[r+i-1][c+j-1]
                    sumP = sumP + pixel
            img[r][c] = sumP / 9
        print "."
    return img
    
    
# Function to perform Gaussian smoothing using a Gaussian filter
# of a specified width and a certain value of sigma. 
# Common values include sigma = 1.0/1.41, width=3/5
#
def img_filterGAUSSIAN(im,sigma,width):

    half = width / 2
    
    # Create a Gaussian kernel
    kernel = numpy.zeros([width,width])
    for i in range(width):
        for j in range(width):
            r = i - half
            s = j - half
            kernel[i][j] = math.exp(-((r**2.0+s**2.0)/(2.0*sigma**2.0)))
    # Normalize the kernel
    kSum = kernel.sum()
    for i in range(width):
        for j in range(width):
            kernel[i,j] = kernel[i,j] / kSum
    
    img = numpy.zeros(im.shape,dtype=numpy.int16)
    
    for r in range(half,im.shape[0]-half):
        for c in range(half,im.shape[1]-half):
            sumP = 0
            for i in range(width):
                for j in range(width):
                    pixel = im[r+i-half][c+j-half]
                    sumP = sumP + pixel * kernel[i][j]
            if sumP > 0:
                img[r][c] = sumP
            else:
                img[r][c] = 0
        #print "."
    return img

# Function to perform unsharp masking
#   Ref(s):
#   Schreiber, W.F., "Wirephoto quality improvement by unsharp masking",
#   Pattern Recognition, Vol.2(2), pp.117-120 (1970)
#
def img_sharpUSM(im):

    kernel = [[0,-1,0],[-1,4,-1],[0,-1,0]]
    kernel = numpy.array(kernel)
    imK = img_filter(im, kernel)
 
    img = numpy.zeros(im.shape,dtype=numpy.int16)
    
    for i in range(0,im.shape[0]):
        for j in range(0,im.shape[1]):
            newPixel = im[i][j] + imK[i][j]
            if newPixel > 255:
                img[i][j] = 255
            else:
                img[i][j] = newPixel
    return img

# Unsharp masking using Gaussian blurring
#
def img_sharpUSMgauss(im,sigma,width):

    imG = img_filterGAUSSIAN(im,sigma,width)
    
    mask = numpy.zeros(im.shape,dtype=numpy.int16)
    img = numpy.zeros(im.shape,dtype=numpy.int16)
    
    for i in range(0,im.shape[0]):
        for j in range(0,im.shape[1]):
            mask[i][j] = im[i][j] - imG[i][j]
            
    for i in range(0,im.shape[0]):
        for j in range(0,im.shape[1]):
            newPixel = im[i][j] + mask[i][j]
            if newPixel <= 255:
                img[i][j] = newPixel
            else:
                img[i][j] = 255
            
    return img  


# Function to perform cubic unsharp masking. A good value for lambda is 
# 0.0015, as per suggested in the reference.
#   Ref(s):
#   Ramponi, G., "A cubic unsharp masking technique for contrast
#   enhancement", Signal Processing, Vol.67, pp.211-222 (1998)
#
def img_sharpCUSM(im,mode,lmbda):

    img = numpy.zeros(im.shape,dtype=numpy.int32)
    imK = numpy.zeros(im.shape)

    # Separable cubic unsharp masking (NS-CUM) (Eq.7)
    if mode == 1:
        for i in range(1,im.shape[0]-1):
             for j in range(1,im.shape[1]-1):        
                 imK[i][j] = (im[i-1][j] - im[i+1][j])**2.0 * \
                             (2.0*im[i][j] - im[i-1][j] - im[i+1][j]) + \
                             (im[i][j-1] - im[i][j+1])**2.0 * \
                             (2.0*im[i][j] - im[i][j-1] - im[i][j+1])

    # Non-separable cubic unsharp masking (NS-CUM) (Eq.8)
    if mode == 2:
       for i in range(1,im.shape[0]-1):
             for j in range(1,im.shape[1]-1):        
                 imK[i][j] = (im[i-1][j] + im[i+1][j] - im[i][j-1] - im[i][j+1])**2.0 * \
                        (4 * im[i][j] - im[i-1][j] - im[i+1][j] - im[i][j-1] - im[i][j+1])
    
    for i in range(0,im.shape[0]):
        for j in range(0,im.shape[1]):
            newPixel = im[i][j] + imK[i][j]*lmbda
            if newPixel < 0:
                img[i][j] = 0
            elif newPixel > 255:
                img[i][j] = 255
            else:
                img[i][j] = newPixel
    
    return img
    
