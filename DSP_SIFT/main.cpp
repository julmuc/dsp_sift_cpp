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

#include <fstream>
#include <platformstl\performance\performance_counter.hpp>		// for profiling
/**************************************************** Main ***********************************************************/

int main(int argc, char** argv)
{	
	platformstl::performance_counter ctimer;
	dspsift_helperlib::dsp_times _dsptimes;
	const char* filename = "C:/Users/Julian/Documents/Visual Studio 2010/Projects/dsp_sift_cpp/times.txt";
	// write times to file
	std::ofstream fout(filename);

	int algorithm = 1; // 1 == dsp_sift; 2 == vl_sift

	//// set vl_sift options 
	////vlfeat_helperlib::vl_sift_options _vlsift_opt = vlfeat_helperlib::vl_sift_options();
	for(int iterations=0; iterations<1; iterations++)
	{	
		std::cout << "Run: " << iterations << std::endl;

		cv::Mat dsp_descr, dsp_features;
		// set dsp options
		dspsift_helperlib::dspOptions dsp_opt = dspsift_helperlib::dspOptions();
		//dsp_opt.vlsift_opt.verbose = 0; // no cout output

		// load template image:
		// needs to be grayscale image
		// IplImage* Timage = cvLoadImage("C:/Users/Julian/Documents/Visual Studio 2010/Projects/dsp_sift_cpp/Debug/Lena.png",0);
		IplImage* Timage = cvLoadImage("C:/Users/Julian/Documents/Visual Studio 2010/Projects/dsp_sift_cpp/Debug/img1.jpg",0);
		
		// static allocation and zero-initialization of arrays for features and descriptors
		double* TFrames = (double*)calloc(4*50000, sizeof(double));
		float* TDescr  = (float*)calloc(128*50000, sizeof(float));


		// stores number of features
		int Tnframes = 0;

		// call sift
		if(algorithm == 2)
		{
			vlfeat_helperlib::vlsift(Timage, TDescr, TFrames, &Tnframes, dsp_opt.vlsift_opt);
		}
		// call dsp_sift
		else if(algorithm == 1)
		{
			ctimer.start();
			dspsift_helperlib::dsp_sift(Timage, dsp_opt, _dsptimes, &Tnframes, dsp_descr, dsp_features);
			ctimer.stop();
			_dsptimes.time_total_dspsift = ctimer.get_microseconds();

			//	///************************ DEBUG ****************/
			//std::cout << "Final Call Descriptors: 1st col OK?" << std::endl;
			//for(int i=0; i<10; i++)
			//{
			//	printf("i: %i \t d_char: %i \n", i, dsp_descr.at<uint8_t>(i,10)); 
			//}
			/////************************ DEBUG ****************/

		}

		// reallocate memory block (in case to much space allocated before) 
		if(algorithm == 2)
		{
			TFrames = (double*)realloc(TFrames, 4*sizeof(double)*Tnframes); // = Y X Scale Angle
			TDescr = (float*)realloc(TDescr, 128*sizeof(float)*Tnframes);

			// draw each feature region as a circle
			for(int i=0; i<Tnframes; i++)
			{
				cvCircle(Timage,													// image
					cvPoint((int)TFrames[0+i*4], (int)TFrames[1+i*4]),				// center (x,y)
					(int)TFrames[2+i*4],											// radius
					cvScalar(255, 0, 0, 0),											// colour
					1,																// thickness
					8,																// linetype
					0);																// shift
			}

			// show image
			cvShowImage("Final Output Features", Timage);
			cvSaveImage("C:/Users/Julian/Pictures/sift_descriptors.png" ,Timage);
		}

		cvWaitKey(0);

		free(TFrames);
		free(TDescr);

		// write times to file
		if(algorithm == 1)
		{
			if(!fout)
			{
				std::cout<<"File Not Opened"<<std::endl;  return -1;
			}

			fout<<_dsptimes.time_total_dspsift<<"\t";
			fout<<_dsptimes.time_vl_sift_normal<<"\t";
			fout<<_dsptimes.time_samplescales<<"\t";
			fout<<_dsptimes.time_getalldescr<<"\t";
			fout<<_dsptimes.time_sort4rows<<"\t";
			fout<<_dsptimes.time_dMat2dArray<<"\t";
			fout<<_dsptimes.time_vl_sift_all<<"\t";
			fout<<_dsptimes.time_sortgenericmat<<"\t";
			fout<<_dsptimes.time_pooldescr<<"\t";
			fout<<_dsptimes.time_normalizehist<<"\t";

			fout<<std::endl;

			std::cout<<"\n total dsp time"<<_dsptimes.time_total_dspsift<<std::endl;
		}
	}

	fout.close();
	return 0;
}