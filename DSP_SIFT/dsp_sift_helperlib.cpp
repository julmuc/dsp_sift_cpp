/** @internal
 ** @file   dsp_sift.cpp
 ** @brief   dsp_sift based on (c) vlfeat.org & J. Dong
 ** http://vision.ucla.edu/~jingming/proj/dsp/
 ** @author Julian heuser
 **/

/************************************************** Includes *********************************************************/
#include "dsp_sift_helperlib.h"

/******************************************** Function Definitions ***************************************************/

platformstl::performance_counter c;


void dspsift_helperlib::dsp_sift(IplImage* i_image,
									dspsift_helperlib::dspOptions i_opt,
									float* o_DATAdescr,
									double* o_DATAframes,
									int* o_nframes)
{	
	dspsift_helperlib::dspOptions _dsp_opt = i_opt;

	//----------------------------------------------- Detection -----------------------------------------------------//
	// static allocation and zero-initialization of arrays for features(frames) and descriptors
    double* siftFrames = (double*)calloc(4*10000, sizeof(double));
    void* siftDescr  = (vl_uint8*)calloc(128*10000, sizeof(vl_uint8));

	// stores number of features(frames)
    int nframes = 0;
	
	// call sift
    vlfeat_helperlib::vlsift(i_image, siftDescr, siftFrames, &nframes, _dsp_opt.vlsift_opt);

	// reallocate memory block (in case to much space allocated before) 
    siftFrames = (double*)realloc(siftFrames, 4*sizeof(double)*nframes); // = Y X Scale Angle
    siftDescr = (vl_uint8*)realloc(siftDescr, 128*sizeof(vl_uint8)*nframes);

	//-------------------------------------- Sample scales around detection -----------------------------------------//
	int dimFeature = 4;

	cv::Mat featureMat;
	featureMat = cv::Mat::zeros(dimFeature, nframes*_dsp_opt.ns, CV_64F);	// 4x(ns*nf) double matrix

	dspsift_helperlib::sample_scales(siftFrames,&nframes,_dsp_opt,featureMat);

	//------------------------- Compute (un-normalized) SIFT descriptors at each scales -----------------------------//
	
	cv::Mat allfeatureMat, alldescriptorMat;
	
	dspsift_helperlib::get_all_descriptors(	i_image,
											featureMat,
											_dsp_opt,
											allfeatureMat,
											alldescriptorMat,
											o_DATAdescr,
											o_DATAframes);
	
	std::cout << "Output frames in original order" << std::endl;
	for(int row_iter=0; row_iter<allfeatureMat.rows; row_iter++)
	{
		std::cout << "row: " << row_iter << " value: " << allfeatureMat.at<double>(row_iter,19) << std::endl;
	}
	std::cout << "Output descriptors in original order" << std::endl;
	for(int row_iter=0; row_iter<alldescriptorMat.rows; row_iter++)
	{
		std::cout << "row: " << row_iter << " value: " << alldescriptorMat.at<float>(row_iter,19) << std::endl;
	}


	//------------------------------------------- Aggregate and normalize -------------------------------------------//


	// todo





	//--------------------------------------------- clean up and return ---------------------------------------------//
	*o_nframes = allfeatureMat.cols;

	free(siftFrames);
	free(siftDescr);
	
	debug("dsp_sift successfully ended");
	return;
}

void dspsift_helperlib::sample_scales(double* i_DATAframes, int* i_nframes, dspOptions i_opt, cv::Mat &o_sampledfeatureMat)
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
	//for(int c=0; c<5; c++)
	//{
	//	debug(c);
	//	for(int r=0; r<4; r++)
	//	{	
	//		debug(o_sampledfeatureMat.at<double>(r,c));
	//	}
	//}
	/*************************** DEBUG *******************/
	
	debug("sample_scales successfully ended");
	return;
}


void dspsift_helperlib::get_all_descriptors(IplImage* i_image, 
											cv::Mat& i_featureMat, 
											dspOptions i_opt, 
											cv::Mat &o_featureMat,
											cv::Mat &o_descriptorMat,
											float* o_DATAdescr, 
											double* o_DATAframes)
{
	cv::Mat sorted_idx, sorted_idx_back, sorted_featureMat, descriptorMat;
	int num_sampledframes = 0;
	int dimDescriptor = 128;

	// get the indices of the InputMat second row sorted in ascending order
	cv::sortIdx(i_featureMat.row(2), sorted_idx, cv::SORT_EVERY_ROW + cv::SORT_ASCENDING);

	// get the indices of the back assignment
	cv::sortIdx(sorted_idx,sorted_idx_back, cv::SORT_EVERY_ROW + cv::SORT_ASCENDING);

	// fast sorting of matrix with only 4 rows
	dspsift_helperlib::sort_4Row_matrixcolsbyindices(i_featureMat,sorted_idx,sorted_featureMat);
	
	// allocate memory
	double* all_output_frames = (double*)calloc(sorted_featureMat.rows*sorted_featureMat.cols, sizeof(double));
    double* sorted_input_frames = (double*)calloc(sorted_featureMat.rows*sorted_featureMat.cols, sizeof(double));
	void* all_output_desc  = (float*)calloc(128*sorted_featureMat.cols, sizeof(float));
	float* all_descr = (float*)calloc(128*50000, sizeof(float));

	//  transform scaled keypoints mat to double array
	dspsift_helperlib::transformKeypointMat_to_Array(sorted_featureMat,sorted_input_frames);

	// set new options for computation of sift descriptors only (no sift feature detection!)
	i_opt.vlsift_opt.ikeys_provided = true;
	i_opt.vlsift_opt.ikeys = sorted_input_frames;
	i_opt.vlsift_opt.nikeys = sorted_featureMat.cols;
	i_opt.vlsift_opt.floatDescriptors = 1;

	//compute sift descriptors of all scaled features
	vlfeat_helperlib::vlsift(i_image, all_output_desc, all_output_frames, &num_sampledframes, i_opt.vlsift_opt);

	// show image
	// draw each feature region as a circle
    for(int i=0; i<num_sampledframes; i++)
	{	
        cvCircle(i_image,																// image
				 cvPoint((int)all_output_frames[0+i*4], (int)all_output_frames[1+i*4]),			// center (x,y)
				 (int)all_output_frames[2+i*4],												// radius
				 cvScalar(255, 0, 0, 0),												// colour
				 1,																		// thickness
				 8,																		// linetype
				 0);																	// shift
    }
	// draw input sift features in black -- not possible in this function anymore
	//for(int i=0; i<nframes; i++)
	//{
 //       cvCircle(i_image,															// image
	//			 cvPoint((int)siftFrames[0+i*4], (int)siftFrames[1+i*4]),			// center (x,y)
	//			 (int)siftFrames[2+i*4],											// radius
	//			 cvScalar(0, 255, 255),												// colour
	//			 1,																	// thickness
	//			 8,																	// linetype
	//			 0);																// shift
 //   }
    cvShowImage("Sampled Features", i_image);
	
	// save variables// temporarily:
	memcpy(o_DATAframes, all_output_frames, 4*num_sampledframes*sizeof(double));
	memcpy(o_DATAdescr, all_output_desc, 128*num_sampledframes*sizeof(float));
	
	
	memcpy(all_descr, all_output_desc, 128*num_sampledframes*sizeof(float));
	all_descr = (float*)realloc(all_descr, 128*sizeof(float)*num_sampledframes);

	// transform back to matrix 

	// fill sorted_featureMat
	for(int row_iter=0; row_iter<sorted_featureMat.rows; row_iter++)
	{
		double* Mat_row = sorted_featureMat.ptr<double>(row_iter);
		for(int col_iter=0; col_iter<sorted_featureMat.cols; col_iter++)
		{		
			Mat_row[col_iter] = all_output_frames[4*col_iter+row_iter];
		}
	}
	// fill descriptorMat
	descriptorMat = cv::Mat::zeros(dimDescriptor, num_sampledframes, CV_32F);	// 128x(ns*nf) float matrix
	for(int row_iter=0; row_iter<descriptorMat.rows; row_iter++)
	{
		float* Mat_row = descriptorMat.ptr<float>(row_iter);
		for(int col_iter=0; col_iter<descriptorMat.cols; col_iter++)
		{		
			Mat_row[col_iter] = all_descr[128*col_iter+row_iter];
		}
	}

	// sort back to original order
	c.start();
	cv::Mat out_featureMat;
	dspsift_helperlib::sort_4Row_matrixcolsbyindices(sorted_featureMat,sorted_idx_back,o_featureMat);	//only 4 rows
	c.stop();
	std::cout << "time (s): " << c.get_seconds() << std::endl;
	std::cout << "time (ms): " << c.get_milliseconds() << std::endl;
	std::cout << "time (us): " << c.get_microseconds() << std::endl;
	
	c.start();
	cv::Mat out_descriptorMat;
	dspsift_helperlib::sort_genericf32_matrixcolsbyindices(descriptorMat,sorted_idx_back,o_descriptorMat);	//need 128 rows
	c.stop();
	std::cout << "time (s): " << c.get_seconds() << std::endl;
	std::cout << "time (ms): " << c.get_milliseconds() << std::endl;
	std::cout << "time (us): " << c.get_microseconds() << std::endl;

	free(all_output_frames);
	free(all_output_desc);
	free(sorted_input_frames);
	free(all_descr);

	debug("get_all_descriptors successfully ended");
	return;
}





void dspsift_helperlib::sort_4Row_matrixcolsbyindices(cv::Mat &i_mat, cv::Mat &i_indices, cv::Mat &o_mat)
{
	//cv::Mat sorted_featureMat;
	o_mat = cv::Mat::zeros(i_mat.rows,i_mat.cols,i_mat.type());

	// row pointers of sorted output mat
	double *p_r0 = 0;
	double *p_r1 = 0;
	double *p_r2 = 0;
	double *p_r3 = 0; 
	p_r0 = o_mat.ptr<double>(0);
	p_r1 = o_mat.ptr<double>(1);
	p_r2 = o_mat.ptr<double>(2);
	p_r3 = o_mat.ptr<double>(3);

	// row pointers of unsorted input mat
	double *p_r0_unsorted = 0;
	double *p_r1_unsorted = 0;
	double *p_r2_unsorted = 0;
	double *p_r3_unsorted = 0;
	p_r0_unsorted = i_mat.ptr<double>(0);
	p_r1_unsorted = i_mat.ptr<double>(1);
	p_r2_unsorted = i_mat.ptr<double>(2);
	p_r3_unsorted = i_mat.ptr<double>(3);

	//sorting
	for(int col_iter=0; col_iter<o_mat.cols; col_iter++)
	{	
		p_r0[col_iter] = p_r0_unsorted[i_indices.at<int>(0,col_iter)];
		p_r1[col_iter] = p_r1_unsorted[i_indices.at<int>(0,col_iter)];
		p_r2[col_iter] = p_r2_unsorted[i_indices.at<int>(0,col_iter)];
		p_r3[col_iter] = p_r3_unsorted[i_indices.at<int>(0,col_iter)];
	}
	return;
}

void dspsift_helperlib::sort_genericf32_matrixcolsbyindices(cv::Mat &i_mat, cv::Mat &i_indices, cv::Mat &o_mat)
{	
	o_mat = cv::Mat::zeros(i_mat.rows,i_mat.cols,i_mat.type());

	float *row_out_ptr = 0;
	float *row_input_ptr = 0;

	for(int row_iter=0; row_iter<o_mat.rows; row_iter++)
	{
		row_out_ptr = o_mat.ptr<float>(row_iter);
		row_input_ptr =  i_mat.ptr<float>(row_iter);

		for(int col_iter=0; col_iter<o_mat.cols; col_iter++)
		{	
			row_out_ptr[col_iter] = row_input_ptr[i_indices.at<int>(0,col_iter)];
		}
	}
	return;
}



void dspsift_helperlib::transformKeypointMat_to_Array(cv::Mat &i_Dmat, double* o_pDarray)
{
	for(int nframe=0; nframe<i_Dmat.cols; nframe++)
	{
		o_pDarray[4*nframe + 0] = i_Dmat.at<double>(0,nframe);
		o_pDarray[4*nframe + 1] = i_Dmat.at<double>(1,nframe);
		o_pDarray[4*nframe + 2] = i_Dmat.at<double>(2,nframe);
		o_pDarray[4*nframe + 3] = i_Dmat.at<double>(3,nframe);

		//if(nframe<5)
		//{	
		//	std::cout << "col: " << nframe << std::endl;
		//	std::cout << "x: " << i_Dmat.at<double>(0,nframe) << std::endl;
		//	std::cout << "x via ptr: " << o_pDarray[4*nframe + 0] << std::endl;
		//	std::cout << "y: " << i_Dmat.at<double>(1,nframe) << std::endl;
		//	std::cout << "scale: " << i_Dmat.at<double>(2,nframe) << std::endl;
		//	std::cout << "angle: " << i_Dmat.at<double>(3,nframe) << std::endl;
		//}
	}

}

void dspsift_helperlib::sorttest()
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

	double *row_ptr = 0;
	double *row_input_ptr = 0;

	for(int row_iter=0; row_iter<sorted_mat.rows; row_iter++)
	{
		row_ptr = sorted_mat.ptr<double>(row_iter);
		row_input_ptr =  InputMat.ptr<double>(row_iter);

		for(int col_iter=0; col_iter<sorted_mat.cols; col_iter++)
		{	
			row_ptr[col_iter] = row_input_ptr[sorted_idx.at<int>(0,col_iter)];
			//p_r0[col_iter] = p_r0_unsorted[sorted_idx.at<int>(0,col_iter)];
			//p_r1[col_iter] = p_r1_unsorted[sorted_idx.at<int>(0,col_iter)];
			//p_r2[col_iter] = p_r2_unsorted[sorted_idx.at<int>(0,col_iter)];
			//p_r3[col_iter] = p_r3_unsorted[sorted_idx.at<int>(0,col_iter)];
		}
	}
	std::cout << "sorted_mat = " << std::endl << " " << sorted_mat << std::endl << std::endl;
	std::cout << "\n" << " =========== sortTest() end ========== " << std::endl;

}