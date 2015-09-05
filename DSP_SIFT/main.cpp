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
	// set dsp options
	dspsift_helperlib::dspOptions dsp_opt;
	dsp_opt.ns = 10;
	dsp_opt.sc_max = 2.0f;
	dsp_opt.sc_min = 0.5f;

	// load template image:
	// needs to be grayscale image
    IplImage* Timage = cvLoadImage("C:/Users/Julian/Documents/Visual Studio 2010/Projects/DSP_SIFT/Debug/Lena.png",0);
 
	// static allocation and zero-initialization of arrays for features and descriptors
    double* TFrames = (double*)calloc(4*10000, sizeof(double));
    vl_uint8* TDescr  = (vl_uint8*)calloc(128*10000, sizeof(vl_uint8));
	
	// stores number of features
    int Tnframes = 0;
	
	// call sift
    //vlfeat_helperlib::vlsift(Timage, TDescr, TFrames, &Tnframes);
	
	// call dsp_sift
	dspsift_helperlib::dsp_sift(Timage,dsp_opt,TDescr,TFrames,&Tnframes);

	// reallocate memory block (in case to much space allocated before) 
    TFrames = (double*)realloc(TFrames, 4*sizeof(double)*Tnframes); // = Y X Scale Angle
    TDescr = (vl_uint8*)realloc(TDescr, 128*sizeof(vl_uint8)*Tnframes);
    
	// draw each feature + descriptor region as a circle
    for(int i=0; i<Tnframes; i++){
 
        cvCircle(Timage,													// image
				 cvPoint((int)TFrames[0+i*4], (int)TFrames[1+i*4]),			// center (x,y)
				 (int)TFrames[2+i*4],										// radius
				 cvScalar(255, 0, 0, 0),									// colour
				 1,															// thickness
				 8,															// linetype
				 0);														// shift
    }
	
	// show image
    cvShowImage("FrameT", Timage);
    cvWaitKey(0);

	free(TFrames);
	free(TDescr);

	return 0;
}