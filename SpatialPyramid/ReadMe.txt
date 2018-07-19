========================================================================
Spatial Pyramid Code
Created by Joe Tighe (jtighe@cs.unc.edu) and Svetlana Lazebnik (lazebnik@cs.unc.edu)
1/17/2009

This MATLAB code implements spatial pyramid computation and matching as 
described in the following paper:

S. Lazebnik, C. Schmid, and J. Ponce, ``Beyond Bags of Features: Spatial 
Pyramid Matching for Recognizing Natural Scene Categories,'' CVPR 2006.

========================================================================

2/29/2012:
-Updated sift normalization. Previous update normalized sift across all values causeing 
 sift descriptors in flat regions to effectly be noise.(see sp_gen_sift.m for details)

9/3/2010:
-Upadated sift generation to allow for fast generation of sift descriptor
 at every pixel (see sp_gen_sift.m for details).
-Added progress bars for user feedback
-Slight change to interface to allow for clearner code. See Example.m for usage.


The main function to build the spatial pyramid is BuildPyramid.

BuildPyramid first extracts SIFT descriptors on a regular grid from each 
image. It then runs k-means to find the dictionary. Each sift descriptor 
is given a texton label corresponding to the closest dictionary codeword. 
Finally, the spatial pyramid is generated from these labels.

Each of these steps are split up into individual functions and can be called 
independently, provided the data from the previous step is stored in the 
correct location. The functions are as follows:

GenerateSiftDescriptors
CalculateDictionary
BuildHistograms
CompilePyramid

NOTE: This code does not include functionality for SVM classification, though it
does include functions for computing the histogram intersection kernel matrix
(hist_isect.m and hist_isect_c.c). For classification, we have used the svm_v0.55
MATLAB toolbox: 

http://theoval.sys.uea.ac.uk/~gcc/svm/toolbox

However, any other SVM package (and kernels other than histogram intersection) 
can be adapted.
