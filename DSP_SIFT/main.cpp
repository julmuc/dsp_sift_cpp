/** @internal
 ** @file   main.cpp
 ** @brief   DSP_SIFT based on (c) 2014-2015, Jingming Dong & 
 ** http://synaptic-activity.blogspot.de/2012/02/vlfeat-sift-with-opencv-code.html
 ** @author Julian heuser
 **/

/************************************************** Macros ***********************************************************/
#define NDEBUG


/************************************************** Includes *********************************************************/
#include "vlfeat_helperlib.h"
#include "dsp_sift_helperlib.h"

#include <opencv2\core\core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>

/**************************************************** Main ***********************************************************/

int main(int argc, char** argv)
{	
	//// set vl_sift options 
	////vlfeat_helperlib::vl_sift_options _vlsift_opt = vlfeat_helperlib::vl_sift_options();

	// set dsp options
	dspsift_helperlib::dspOptions dsp_opt = dspsift_helperlib::dspOptions();

	// load template image:
	// needs to be grayscale image
   // IplImage* Timage = cvLoadImage("C:/Users/Julian/Documents/Visual Studio 2010/Projects/dsp_sift_cpp/Debug/Lena.png",0);
	IplImage* Timage = cvLoadImage("C:/Users/Julian/Documents/Visual Studio 2010/Projects/dsp_sift_cpp/Debug/img1.jpg",0);
	// static allocation and zero-initialization of arrays for features and descriptors
    double* TFrames = (double*)calloc(4*50000, sizeof(double));
    float* TDescr  = (float*)calloc(128*50000, sizeof(float));
	
	// stores number of features
    int Tnframes = 0;
	
	//// call sift
	////vlfeat_helperlib::vlsift(Timage, TDescr, TFrames, &Tnframes, _vlsift_opt);
	
	// call dsp_sift
	dspsift_helperlib::dsp_sift(Timage,dsp_opt,TDescr,TFrames,&Tnframes);

	// reallocate memory block (in case to much space allocated before) 
    TFrames = (double*)realloc(TFrames, 4*sizeof(double)*Tnframes); // = Y X Scale Angle
    TDescr = (float*)realloc(TDescr, 128*sizeof(float)*Tnframes);
    
	// draw each feature region as a circle
    for(int i=0; i<Tnframes; i++)
	{
        cvCircle(Timage,													// image
				 cvPoint((int)TFrames[0+i*4], (int)TFrames[1+i*4]),			// center (x,y)
				 (int)TFrames[2+i*4],										// radius
				 cvScalar(255, 0, 0, 0),									// colour
				 1,															// thickness
				 8,															// linetype
				 0);														// shift
    }
	
	// show image
    cvShowImage("Final Output Features", Timage);
    cvWaitKey(0);

	free(TFrames);
	free(TDescr);

	return 0;
}