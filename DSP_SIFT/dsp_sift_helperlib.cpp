/** @internal
 ** @file   dsp_sift.cpp
 ** @brief   dsp_sift based on (c) vlfeat.org & J. Dong
 ** http://vision.ucla.edu/~jingming/proj/dsp/
 ** @author Julian heuser
 **/

/************************************************** Includes *********************************************************/
#include "dsp_sift_helperlib.h"

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

	// stores number of features(frames)
    int nframes = 0;
	
	// call sift
    vlfeat_helperlib::VLSIFT(i_image, siftDescr, siftFrames, &nframes);
	
	// save variables:
	memcpy(o_DATAframes, siftFrames, 4*nframes*sizeof(double));
	memcpy(o_DATAdescr, siftDescr, 128*nframes*sizeof(vl_uint8));

	// reallocate memory block (in case to much space allocated before) 
    siftFrames = (double*)realloc(siftFrames, 4*sizeof(double)*nframes); // = Y X Scale Angle
    siftDescr = (vl_uint8*)realloc(siftDescr, 128*sizeof(vl_uint8)*nframes);

	//-------------------------------------- Sample scales around detection -----------------------------------------//
	int dimFeature = 4;

	cv::Mat featureMat;
	featureMat = cv::Mat::zeros(dimFeature, nframes*i_opt.ns, CV_64F);	// 4x(ns*nf) double matrix

	sampleScales(siftFrames,&nframes,i_opt,featureMat);

	//--------------------------------- Compute un-normalized SIFT at each scales -----------------------------------//


	// todo

	platformstl::performance_counter c;

	/******************** DEBUG ******************/
	//c.start();
	//dspsift_helperlib::sortTest();
	//c.stop();
	//   
	//std::cout << "time (s): " << c.get_seconds() << std::endl;
	//   std::cout << "time (ms): " << c.get_milliseconds() << std::endl;
	//   std::cout << "time (us): " << c.get_microseconds() << std::endl;
	/******************** DEBUG ******************/
	
	c.start();

	cv::Mat sorted_featureMat;
	sorted_featureMat = cv::Mat::zeros(dimFeature, nframes*i_opt.ns, CV_64F);

	cv::Mat sorted_idx, sorted_idx_back, sorted_mat;

	// get the indices of the InputMat second row data sorted in ascending order
	cv::sortIdx(featureMat.row(2), sorted_idx, cv::SORT_EVERY_ROW + cv::SORT_ASCENDING);

	// get the indices of the back assignment (so far just for testing)
	cv::sortIdx(sorted_idx,sorted_idx_back, cv::SORT_EVERY_ROW + cv::SORT_ASCENDING);

	// row pointers of sorted mat
	double *p_r0 = 0;
	double *p_r1 = 0;
	double *p_r2 = 0;
	double *p_r3 = 0; 
	p_r0 = sorted_featureMat.ptr<double>(0);
	p_r1 = sorted_featureMat.ptr<double>(1);
	p_r2 = sorted_featureMat.ptr<double>(2);
	p_r3 = sorted_featureMat.ptr<double>(3);

	// row pointers of unsorted mat
	double *p_r0_unsorted = 0;
	double *p_r1_unsorted = 0;
	double *p_r2_unsorted = 0;
	double *p_r3_unsorted = 0;
	p_r0_unsorted = featureMat.ptr<double>(0);
	p_r1_unsorted = featureMat.ptr<double>(1);
	p_r2_unsorted = featureMat.ptr<double>(2);
	p_r3_unsorted = featureMat.ptr<double>(3);

	//sorting
	for(int col_iter=0; col_iter<sorted_featureMat.cols; col_iter++)
	{	
		p_r0[col_iter] = p_r0_unsorted[sorted_idx.at<int>(0,col_iter)];
		p_r1[col_iter] = p_r1_unsorted[sorted_idx.at<int>(0,col_iter)];
		p_r2[col_iter] = p_r2_unsorted[sorted_idx.at<int>(0,col_iter)];
		p_r3[col_iter] = p_r3_unsorted[sorted_idx.at<int>(0,col_iter)];
	}


	c.stop();
	std::cout << "time (s): " << c.get_seconds() << std::endl;
	std::cout << "time (ms): " << c.get_milliseconds() << std::endl;
	std::cout << "time (us): " << c.get_microseconds() << std::endl;
	
	
	//------------------------------------------- Aggregate and normalize -------------------------------------------//


	// todo





	//--------------------------------------------- clean up and return ---------------------------------------------//
	*o_nframes = nframes;

	free(siftFrames);
	free(siftDescr);

	return;
}

void dspsift_helperlib::sampleScales(double* i_DATAframes, int* i_nframes, dspOptions i_opt, cv::Mat &o_sampledfeatureMat)
{	
	/******
	Eingabedaten: Anzahl der Merkmalspunkte, Merkmalspunkte(4 Dimensionen Vector), Anzahl & Range Skalierungen
	
	Verarbeitungsschritt
	Glätte die Region um jeden extrahierten Merkmalspunkt für "ns" verschiedene Skalierungen,
	das ergibt "ns" mal die Anzahl der Eingangsmerkmalspunkte
	Jeder Merkmalspunkt hat 4 Dimensionen
	
	Ausgabedaten: Array/Matrix der Größe 4 * num_features * num_scales
	******/


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
		debug(c);
		for(int r=0; r<4; r++)
		{	
			debug(o_sampledfeatureMat.at<double>(r,c));
		}
	}
	/*************************** DEBUG *******************/

	return;
}


void dspsift_helperlib::sortTest()
{	
	std::cout << "\n" << " =========== sortTest() start ========== " << std::endl;

	cv::Mat sorted_idx, assign_back, sorted_mat;

	// initialize Test input matrix
	cv::Mat InputMat = (cv::Mat_<double>(4,3) <<   1,    2,   3,
											4,    2,   1,
										  0.5,  1.2, 0.3,
										 -0.2,  0.2, 0.1);
	// initialize Test output sorted matrix
	sorted_mat = cv::Mat::zeros(4,3, CV_64F);


	std::cout << "InputMat = " << std::endl << " " << InputMat << std::endl << std::endl;

	// get the indices of the InputMat second row data sorted in ascending order
	cv::sortIdx(InputMat.row(2), sorted_idx, cv::SORT_EVERY_ROW + cv::SORT_ASCENDING);

	// get the indices of the back assignment (so far just for testing)
	cv::sortIdx(sorted_idx,assign_back, cv::SORT_EVERY_ROW + cv::SORT_ASCENDING);

	std::cout << "InputMat.row(2) = " << std::endl << " " << InputMat.row(2) << std::endl << std::endl;
	std::cout << "sorted_idx = " << std::endl << " " << sorted_idx << std::endl << std::endl;
	std::cout << "assign_back = " << std::endl << " " << assign_back << std::endl << std::endl;

	// row pointers of sorted mat
	double *p_r0 = 0;
	double *p_r1 = 0;
	double *p_r2 = 0;
	double *p_r3 = 0; 
	p_r0 = sorted_mat.ptr<double>(0);
	p_r1 = sorted_mat.ptr<double>(1);
	p_r2 = sorted_mat.ptr<double>(2);
	p_r3 = sorted_mat.ptr<double>(3);

	// row pointers of unsorted mat
	double *p_r0_unsorted = 0;
	double *p_r1_unsorted = 0;
	double *p_r2_unsorted = 0;
	double *p_r3_unsorted = 0;
	p_r0_unsorted = InputMat.ptr<double>(0);
	p_r1_unsorted = InputMat.ptr<double>(1);
	p_r2_unsorted = InputMat.ptr<double>(2);
	p_r3_unsorted = InputMat.ptr<double>(3);


	for(int col_iter=0; col_iter<sorted_mat.cols; col_iter++)
	{	
		p_r0[col_iter] = p_r0_unsorted[sorted_idx.at<int>(0,col_iter)];
		p_r1[col_iter] = p_r1_unsorted[sorted_idx.at<int>(0,col_iter)];
		p_r2[col_iter] = p_r2_unsorted[sorted_idx.at<int>(0,col_iter)];
		p_r3[col_iter] = p_r3_unsorted[sorted_idx.at<int>(0,col_iter)];
	}
	std::cout << "sorted_mat = " << std::endl << " " << sorted_mat << std::endl << std::endl;
	std::cout << "\n" << " =========== sortTest() end ========== " << std::endl;

}