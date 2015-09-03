/** @internal
 ** @file   dsp_sift.cpp
 ** @brief   dsp_sift based on (c) vlfeat.org & J. Dong
 ** http://vision.ucla.edu/~jingming/proj/dsp/
 ** @author Julian heuser
 **/

/************************************************** Includes *********************************************************/
#include "dsp_sift_helperlib.h"
#include <iostream>

/******************************************** Function Definitions ***************************************************/

void dspsift_helperlib::DSP_SIFT(IplImage* i_image,
									dspsift_helperlib::dspOptions i_opt,
									vl_uint8* o_DATAdescr,
									double* o_DATAframes,
									int* o_nframes)
{	

	//----------------------------------------------- Detection -----------------------------------------------------//
	// static allocation and zero-initialization of arrays for features(frames) and descriptors
    double* siftFrames = (double*)calloc(4*10000, sizeof(double));
    vl_uint8* siftDescr  = (vl_uint8*)calloc(128*10000, sizeof(vl_uint8));
	
	o_DATAframes = (double*)calloc(4*10000, sizeof(double));
    o_DATAdescr  = (vl_uint8*)calloc(128*10000, sizeof(vl_uint8));


	// stores number of features(frames)
    int nframes = 0;
	
	// call sift
    vlfeat_helperlib::VLSIFT(i_image, siftDescr, siftFrames, &nframes);
	
	// reallocate memory block (in case to much space allocated before) 
    siftFrames = (double*)realloc(siftFrames, 4*sizeof(double)*nframes); // = Y X Scale Angle
    siftDescr = (vl_uint8*)realloc(siftDescr, 128*sizeof(vl_uint8)*nframes);

	//-------------------------------------- Sample scales around detection -----------------------------------------//
	// Eingabedaten: Anzahl der Merkmalspunkte, Merkmalspunkte(4 Dimensionen Vector)
	
	// Verarbeitungsschritt
	// Glätte die Region um jeden extrahierten Merkmalspunkt für "ns" verschiedene Skalierungen,
	// das ergibt "ns" mal die Anzahl der Eingangsmerkmalspunkte
	// Jeder Merkmalspunkt hat 4 Dimensionen ->
	// Ausgabedaten: Array/Matrix der Größe 4 * num_features * num_scales(== ns)
	// Die 3. Dimension gibt dann jeweils die skalierten Merkmalspunkte an


	// todo
	int dimFeature = 4;

	cv::Mat featureMat;
	featureMat = cv::Mat::zeros(dimFeature, nframes*i_opt.ns, CV_64F);	// 4x(ns*nf) matrix

	sampleScales(siftFrames,&nframes,i_opt,featureMat);

	//--------------------------------- Compute un-normalized SIFT at each scales -----------------------------------//


	// todo


	//------------------------------------------- Aggregate and normalize -------------------------------------------//


	// todo


	//tmp return values
	o_DATAframes = (double*)realloc(o_DATAframes, 4*sizeof(double)*nframes);
    o_DATAdescr  = (vl_uint8*)realloc(o_DATAdescr, 128*sizeof(vl_uint8)*nframes);

	*o_nframes = nframes;

	return;
}

void dspsift_helperlib::sampleScales(double* i_DATAframes, int* i_nframes, dspOptions i_opt, cv::Mat &o_sampledfeatureMat)
{	
	// generate linear spaced scales
	double scale_diff;
	int scale_counter = 0;
	int feature_batch_counter = 0;
	std::vector<double> scales(i_opt.ns);

	scale_diff = (i_opt.sc_max - i_opt.sc_min)/(i_opt.ns - 1);
	for(int scale_iter=0; scale_iter<i_opt.ns; scale_iter++)
	{
		scales.at(scale_iter) = i_opt.sc_min + scale_iter * scale_diff;
	}

	// fill sampledfeatureMat
	for(int row_iter=0; row_iter<o_sampledfeatureMat.rows; row_iter++)
	{
		double* Mat_row = o_sampledfeatureMat.ptr<double>(row_iter);
		for(int col_iter=0; col_iter<o_sampledfeatureMat.cols; col_iter++)
		{	
			if(row_iter != 2)
			{	
				Mat_row[col_iter] = i_DATAframes[4*(col_iter%(*i_nframes))+row_iter];
			}
			else
			{
				Mat_row[col_iter] = scales.at(scale_counter) * i_DATAframes[4*(col_iter%(*i_nframes))+row_iter];
				feature_batch_counter++;
				if(feature_batch_counter >= *i_nframes)
				{
					feature_batch_counter = 0;
					scale_counter++;
				}
			}
		}
	}
	
	/*************************** DEBUG *******************/
	for(int c=0; c<5; c++)
	{
		std::cout << "feature: " << c << std::endl;
		for(int r=0; r<4; r++)
		{	
			if(r==0)
				std::cout << "x: " << o_sampledfeatureMat.at<double>(r,c) << std::endl;
			else if(r==1)
				std::cout << "y: " <<  o_sampledfeatureMat.at<double>(r,c) << std::endl;
			else if(r==2)
				std::cout << "sigma: " <<  o_sampledfeatureMat.at<double>(r,c) << std::endl;			// needs to be different, as its now scaled accordingly!
			else if(r==3)
				std::cout << "angle: " <<  o_sampledfeatureMat.at<double>(r,c) << "\n" << std::endl;
		}
	}
	/*************************** DEBUG *******************/

	return;
}