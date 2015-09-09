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
									int* o_nframes,
									cv::Mat &o_descr,
									cv::Mat &o_features)
{	
	dspsift_helperlib::dspOptions _dsp_opt = i_opt;

	int dim_descriptor = 128;

	//----------------------------------------------- Detection -----------------------------------------------------//
	// static allocation and zero-initialization of arrays for features(frames) and descriptors
    double* siftFrames = (double*)calloc(4*10000, sizeof(double));
    void* siftDescr  = (vl_uint8*)calloc(dim_descriptor*10000, sizeof(vl_uint8));

	// stores number of features(frames)
    int nframes = 0;
	
	// call sift
    vlfeat_helperlib::vlsift(i_image, siftDescr, siftFrames, &nframes, _dsp_opt.vlsift_opt);

	// reallocate memory block (in case to much space allocated before) 
    siftFrames = (double*)realloc(siftFrames, 4*sizeof(double)*nframes); // = Y X Scale Angle
    siftDescr = (vl_uint8*)realloc(siftDescr, dim_descriptor*sizeof(vl_uint8)*nframes);

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
											alldescriptorMat);
	
	/*********************************** DEBUG *****************************/
	//std::cout << "Output frames in original order" << std::endl;
	//for(int row_iter=0; row_iter<allfeatureMat.rows; row_iter++)
	//{
	//	std::cout << "row: " << row_iter << " value: " << allfeatureMat.at<double>(row_iter,19) << std::endl;
	//}
	//std::cout << "Output descriptors in original order" << std::endl;
	//for(int row_iter=0; row_iter<alldescriptorMat.rows; row_iter++)
	//{
	//	std::cout << "row: " << row_iter << " value: " << alldescriptorMat.at<float>(row_iter,19) << std::endl;
	//}
	/*********************************** DEBUG *****************************/


	//------------------------------------------- pool and normalize -------------------------------------------//

	int batchsize = allfeatureMat.cols/_dsp_opt.ns;
	cv::Mat out_featureMat, firstbatch;

	firstbatch = allfeatureMat.colRange(0,batchsize).rowRange(0,allfeatureMat.rows);	// 4 x batchsize
	firstbatch.copyTo(out_featureMat);

	std::cout << "out_featureMat.cols: " << out_featureMat.cols << " .rows: " << out_featureMat.rows << std::endl;

	int row_nr_scale = 2;
	double* Mat_row = out_featureMat.ptr<double>(row_nr_scale);
	for(int col_iter=0; col_iter<batchsize; col_iter++)
	{
		Mat_row[col_iter] = siftFrames[4*(col_iter)+row_nr_scale];
	}
	


	const char* filename = "C:/Users/Julian/Pictures/featrs.txt";
    writeMatToFile(out_featureMat,filename);
	
	//std::cout << "frames after overwriting " << std::endl;
	//for(int row_iter=0; row_iter<out_featureMat.rows; row_iter++)
	//{
	//	std::cout << "row: " << row_iter << " value: " << out_featureMat.at<double>(row_iter,0) << std::endl;
	//}
	//for(int row_iter=0; row_iter<out_featureMat.rows; row_iter++)
	//{
	//	std::cout << "row: " << row_iter << " value: " << out_featureMat.at<double>(row_iter,out_featureMat.cols-1) << std::endl;
	//}

	cv::Mat out_pooledDescriptors;
	// pooling --- working !
	c.start();
	dspsift_helperlib::pool_descriptors(alldescriptorMat,batchsize,_dsp_opt.ns,out_pooledDescriptors);
	c.stop();
	std::cout << "time (s): " << c.get_seconds() << std::endl;
	std::cout << "time (ms): " << c.get_milliseconds() << std::endl;
	std::cout << "time (us): " << c.get_microseconds() << std::endl;
	

	/*********************************** DEBUG *****************************/
	//int nf = 3;
	//int ns = 3;
	//cv::Mat InputMat = (cv::Mat_<float>(4,9) <<		1,		2,		3,		2,		3,		4,		3,		4,		5,
	//												1,		2,		3,		2,		3,		4,		3,		4,		5,
	//												1,		2,		3,		2,		3,		4,		3,		4,		5,
	//												1,		2,		3,		2,		3,		4,		3,		4,		5);
	//cv::Mat outMat;
	//std::cout << "InputMat = " << std::endl << " " << InputMat << std::endl << std::endl;
	//dspsift_helperlib::pool_descriptors(InputMat,nf,ns,outMat);
	/*********************************** DEBUG *****************************/

	// normalizing

	//reshape pooledDescriptormat to 1d array
	cv::Mat outdescrArray;
	outdescrArray = cv::Mat(out_pooledDescriptors.t()).reshape(1,1);

	//std::cout << "rows: " << outdescrArray.rows << " cols: " << outdescrArray.cols << std::endl;
	//for(int col_iter=0; col_iter<130; col_iter++)
	//{
	//	std::cout << "col: " << col_iter << " value: " << outdescrArray.at<float>(0,col_iter) << std::endl;
	//}

	// TODO

	float *d = (float*)outdescrArray.data;
	// Pointer to the 0-th row

	float *d_out =  (float*)calloc(outdescrArray.cols, sizeof(float));
	uint8_t *d_char =  (uint8_t*)calloc(outdescrArray.cols, sizeof(uint8_t));

	//d_out = (float*)outdescrArray.data;

	// TODO: use defines here
	const int NBO = 8;
	const int NBP = 4;

	int nel = dim_descriptor*batchsize;
            
	for (int i = 0; i < nel; ++i)
	{
		d_out[i] = d[i];
	}

	/************************ DEBUG ****************/
	std::cout << "D_out Input: OK" << std::endl;
	for(int i=128; i<138; i++)
	{
		std::cout << i << " " << "value: " << d_out[i] << std::endl;
	}
	/************************ DEBUG ****************/
	int DEBUGCOL = 1;
	for (int i = 0; i < batchsize; ++i)
	{	
		float norm = dspsift_helperlib::normalize_histogram(d_out, d_out + NBO*NBP*NBP);
		
		/************************ DEBUG ****************/
		if(i==DEBUGCOL)
		{	
			printf("After first normalization: \n");
			for(int j = 0; j< 10; j++)
				printf("j: %i value: %f \n", j, d_out[j]); 
			std::cout << "Norm: " << norm << std::endl;
		}
		/************************ DEBUG ****************/

		for (int bin = 0; bin < NBO*NBP*NBP; ++bin)
		{
			if (d_out[bin] > 0.0667f)
				d_out[bin] = 0.0667f;
			
			/************************ DEBUG ****************/
			if(i==DEBUGCOL && bin<10)
				printf("bin: %i value: %f \n", bin, d_out[bin]);
			/************************ DEBUG ****************/
		}
		
		float norm2 = dspsift_helperlib::normalize_histogram(d_out, d_out + NBO*NBP*NBP);
		
		/************************ DEBUG ****************/
		if(i==DEBUGCOL)
		{	
			printf("After 2nd normalization: \n");
			for(int j = 0; j< 10; j++)
				printf("j: %i value: %f \n", j, d_out[j]); 
			printf("Norm: %f\n",norm2);
		}
		/************************ DEBUG ****************/

		for (int bin = 0; bin < NBO*NBP*NBP; ++bin)
		{
			float x = 512.0F * d_out[bin];
			x = (x < 255.0F) ? x : 255.0F;
			d_char[bin] = (uint8_t)x;
			
			/************************ DEBUG ****************/
			if(i==DEBUGCOL && bin<10)
				printf("x: %f \t d_out: %f \t d_char: %f \n", x, d_out[bin], (float)d_char[bin]);
			/************************ DEBUG ****************/
		}
		d_out += dim_descriptor;
		d_char += dim_descriptor;
	}

	// set pointer back to first element
	d_out = d_out -  batchsize*dim_descriptor;
	d_char = d_char - batchsize*dim_descriptor;

	/************************ DEBUG ****************/
	std::cout << "D_char Output 1st col: OK?" << std::endl;
	for(int i=0; i<10; i++)
	{
		printf("i: %i \t d_out: %f \t d_char: %f \n", i, d_out[i], (float)d_char[i]); 
	}
	d_out += dim_descriptor;
	d_char += dim_descriptor;
	printf("\n 2nd col: OK?\n");
	for(int i=0; i<10; i++)
	{
		printf("i: %i \t d_out: %f \t d_char: %f \n", i, d_out[i], (float)d_char[i]); 
	}
	d_out -= dim_descriptor;
	d_char -= dim_descriptor;
	/************************ DEBUG ****************/


	// todo
	// fill out_pooledandnormalizedDescriptorMat (8 bit unsigned char)
	cv::Mat out_pn_DescriptorMat;
	out_pn_DescriptorMat = cv::Mat::zeros(dim_descriptor, batchsize, CV_8UC1);	// 128x(ns*nf) float matrix
	for(int row_iter=0; row_iter<out_pn_DescriptorMat.rows; row_iter++)
	{
		uint8_t* Mat_row = out_pn_DescriptorMat.ptr<uint8_t>(row_iter);
		for(int col_iter=0; col_iter<out_pn_DescriptorMat.cols; col_iter++)
		{		
			Mat_row[col_iter] = d_char[dim_descriptor*col_iter+row_iter];
		}
	}
	cv::imwrite( "C:/Users/Julian/Pictures/descritors.png", out_pn_DescriptorMat );
	
	
	///************************ DEBUG ****************/
	std::cout << "Final Output Descriptors: 1st col OK?" << std::endl;
	for(int i=0; i<10; i++)
	{
		printf("i: %i \t d_out: %f \t d_char: %i \n", i, d_out[i], out_pn_DescriptorMat.at<uint8_t>(i,10)); 
	}
	///************************ DEBUG ****************/

	//for(int i=0; i<out_featureMat.cols; i++)
	//{	
 //       cvCircle(i_image,																// image
	//			 cvPoint((int)out_featureMat.at<double>(0,i), (int)out_featureMat.at<double>(1,i)),			// center (x,y)
	//			 (int)out_featureMat.at<double>(2,i),												// radius
	//			 cvScalar(255, 0, 0, 0),												// colour
	//			 1,																		// thickness
	//			 8,																		// linetype
	//			 0);																	// shift
 //   }
	//cvShowImage("Final DSP-Sift Features", i_image);
	//cvSaveImage("C:/Users/Julian/Pictures/dsp_sift_final_descriptors.png" ,i_image);

	//--------------------------------------------- clean up and return ---------------------------------------------//
	*o_nframes = out_featureMat.cols;
	o_descr = out_pn_DescriptorMat;
	o_features = out_featureMat;



	free(siftFrames);
	free(siftDescr);
	//free(d_out);

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
											cv::Mat &o_descriptorMat)
{
	cv::Mat sorted_idx, sorted_idx_back, sorted_featureMat, descriptorMat;
	int num_sampledframes = 0;
	int dimDescriptor = 128;

	// get the indices of the InputMat second row sorted in ascending order
	cv::sortIdx(i_featureMat.row(2), sorted_idx, cv::SORT_EVERY_ROW + cv::SORT_ASCENDING);

	// get the indices of the back assignment
	cv::sortIdx(sorted_idx,sorted_idx_back, cv::SORT_EVERY_ROW + cv::SORT_ASCENDING);

	// fast sorting of matrix with only 4 rows
	dspsift_helperlib::s_sort_4rowf64_matrixcolsbyindices(i_featureMat,sorted_idx,sorted_featureMat);
	
	// allocate memory
	double* all_output_frames = (double*)calloc(sorted_featureMat.rows*sorted_featureMat.cols, sizeof(double));
    double* sorted_input_frames = (double*)calloc(sorted_featureMat.rows*sorted_featureMat.cols, sizeof(double));
	void* all_output_desc  = (float*)calloc(128*sorted_featureMat.cols, sizeof(float));
	float* all_descr = (float*)calloc(128*50000, sizeof(float));

	//  transform scaled keypoints mat to double array
	dspsift_helperlib::dmat_to_darray(sorted_featureMat,sorted_input_frames);

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
	//cvSaveImage("C:/Users/Julian/Pictures/dsp_sift_descriptors.png" ,i_image);
	
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
	dspsift_helperlib::s_sort_4rowf64_matrixcolsbyindices(sorted_featureMat,sorted_idx_back,o_featureMat);	//only 4 rows
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

void dspsift_helperlib::pool_descriptors(cv::Mat& i_descriptorMat, int batchsize, int ns, cv::Mat& o_pooleddescriptorMat)
{
	int dim_d = i_descriptorMat.rows;
	o_pooleddescriptorMat = cv::Mat::zeros(i_descriptorMat.rows,batchsize,i_descriptorMat.type());
	
	for(int i=0; i<batchsize; i++)
	{
		cv::Mat tmp = cv::Mat::zeros(dim_d,ns,CV_32F);
		for(int j=0; j<ns; j++)
		{
			i_descriptorMat.col(batchsize*j+i).copyTo(tmp.col(j));
		}
		//std::cout << "tmp = " << std::endl << " " << tmp << std::endl << std::endl;
		cv::reduce(tmp,o_pooleddescriptorMat.col(i),1,CV_REDUCE_AVG,-1);
	}
	//std::cout << "pooledMat = " << std::endl << " " << o_pooleddescriptorMat << std::endl << std::endl;
	debug("pool_descriptors successfully ended");
	return;
}



void dspsift_helperlib::s_sort_4rowf64_matrixcolsbyindices(cv::Mat &i_mat, cv::Mat &i_indices, cv::Mat &o_mat)
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



void dspsift_helperlib::dmat_to_darray(cv::Mat &i_Dmat, double* o_pDarray)
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

static inline float dspsift_helperlib::normalize_histogram(float* begin, float* end)
{
	float* iter;
	float norm = 0.0f;

	for (iter = begin; iter != end; ++iter)
		norm += (*iter)*(*iter);

	norm = vl_fast_sqrt_f(norm) + VL_EPSILON_F;
	//norm = sqrt(norm) + VL_EPSILON_F;

	for (iter = begin; iter != end; ++iter)
		*iter /= norm;

	return norm;
}

/*********************************************************************************************************************/

static inline float dspsift_helperlib::vl_fast_resqrt_f(float x)
{
	/* 32-bit version */
	union 
	{
		float x;
		int  i;
	}u;

	float xhalf = (float)0.5*x;

	/* convert floating point value in RAW integer */
	u.x = x;

	/* gives initial guess y0 */
	u.i = 0x5f3759df - (u.i >> 1);

	/* two Newton steps */
	u.x = u.x * ( (float) 1.5  - xhalf*u.x*u.x);
	u.x = u.x * ( (float) 1.5  - xhalf*u.x*u.x);
	return u.x;
}

/*********************************************************************************************************************/

static inline float dspsift_helperlib::vl_fast_sqrt_f(float x)
{
	return (x < 1e-8) ? 0 : x * vl_fast_resqrt_f(x);
}

void dspsift_helperlib::writeMatToFile(cv::Mat& m, const char* filename)
{
    std::ofstream fout(filename);

    if(!fout)
    {
        std::cout<<"File Not Opened"<<std::endl;  return;
    }

    for(int i=0; i<m.rows; i++)
    {
        for(int j=0; j<m.cols; j++)
        {
            fout<<m.at<double>(i,j)<<"\t";
        }
        fout<<std::endl;
    }

    fout.close();
}