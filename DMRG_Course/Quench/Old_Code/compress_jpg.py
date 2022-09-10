import numpy as np
from scipy import linalg
from PIL import Image as image_tool

filename = "dogs.bmp"
chi      = 25

image = image_tool.open(filename)
image = np.asarray(image) 

# image is a matrix of size ( length, width, 3 )

U, S, V     = linalg.svd( image[:,:,0] )
red_compr   = np.dot( U[:,:chi] * S[:chi] , V[:chi,:] )

U, S, V     = linalg.svd( image[:,:,1] )
green_compr = np.dot( U[:,:chi] * S[:chi] , V[:chi,:] )

U, S, V     = linalg.svd( image[:,:,2] )
blue_compr  = np.dot( U[:,:chi] * S[:chi] , V[:chi,:] )

compr        = np.zeros( image.shape, np.float16 )
compr[:,:,0] = red_compr
compr[:,:,1] = green_compr
compr[:,:,2] = blue_compr

compr[ compr < 0.   ] = 0.
compr[ compr > 255. ] = 255.

compr = compr.astype(np.uint8)
compr = image_tool.fromarray(compr)

size_orig = 1.0 * 3.0 *   image.shape[0] * image.shape[1]
size_comp = 2.0 * 3.0 * ( image.shape[0]*chi + chi + chi*image.shape[1] )

print ("size before / after compression: {0:.2f} / {1:.2f} MB".format(size_orig/1024**2, size_comp/1024**2))

compr.show()

