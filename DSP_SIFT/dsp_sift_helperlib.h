/** @internal
 ** @file   dsp_sift.h
 ** @brief   dsp_sift based on (c) vlfeat.org & J. Dong
 ** http://vision.ucla.edu/~jingming/proj/dsp/
 ** @author Julian heuser
 **/

#ifndef DSPSIFTHELPERLIB_H
#define DSPSIFTHELPERLIB_H

//#define NDEBUG	// define for no debug messages

/************************************************** Includes *********************************************************/
#include <opencv2\core\core.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/highgui/highgui.hpp>
#include "vlfeat_helperlib.h"
#include "debug_helper.h"
#include <platformstl\performance\performance_counter.hpp>		// for profiling
#include <stdint.h>
#include <fstream>

/************************************************** Structures *******************************************************/


/**************************************** Begin Namespace dspsift_helperlib ******************************************/
namespace dspsift_helperlib
{

	typedef struct dspOptions
	{	
		// default constructor
		dspOptions(): sc_min(0.5), sc_max(2), ns(10), vlsift_opt(vlfeat_helperlib::vl_sift_options()){}

		double sc_min;			// scale sampling lower limit
		double sc_max;			// scale sampling upper limit. Scales are sampled from (sc_min * s, sc_max * s) where s is the detected scale
		int ns;					// number of scales
		vlfeat_helperlib::vl_sift_options vlsift_opt;
	} dspOptions ;

	typedef struct dsp_times
	{	
		long long time_total_dspsift;
		long long time_vl_sift_normal;
		long long time_vl_sift_all;
		long long time_samplescales;
		long long time_getalldescr;
		long long time_pooldescr;
		long long time_sort4rows;
		long long time_sortgenericmat;
		long long time_dMat2dArray;
		long long time_normalizehist;

	} dsp_times ;

	/******************************************** Function Declarations **********************************************/

	/** ------------------------------------------------------------------
	** @internal
	** @brief based on dsp_sift by J.Dong
	**
	** @param i_image input greyscale image
	** @param i_opt input option for sampling different scales
	** @param o_dsptimes timing struct
	** @param o_descr cv:.Mat (float32) of output descriptors
	** @param o_features cv::Mat (float64) of output features(frames)
	** @param o_nframes int* to the number of frames
	** 
	**/
	void dsp_sift(IplImage* i_image, dspOptions i_opt, dsp_times &o_dsptimes, int* o_nframes, cv::Mat &o_descr, cv::Mat &o_features);

	/** ------------------------------------------------------------------
	** @internal
	** @brief Sample scales around detection
	**
	** @param i_DATAframes input features 
	** @param i_nframes total number of features
	** @param i_opt scale options
	** @param o_sampledfeatureMat sampled output matrix
	** 
	**/
	void sample_scales(double* i_DATAframes, int* i_nframes, dspOptions i_opt, cv::Mat &o_sampledfeatureMat);

	/** ------------------------------------------------------------------
	** @internal
	** @brief Computes all sift descriptors at all scales
	**
	** @param i_image input greyscale image
	** @param i_featureMat input matrix containing all sift features (as columns) 
	** @param i_opt scale options
	** @param o_featureMat output feature matrix (features as columns)
	** @param o_descriptorMat output descriptor matrix (descriptors as columns)
	** 
	**/
	void get_all_descriptors(IplImage* i_image,
							cv::Mat& i_featureMat,
							dspOptions i_opt,
							cv::Mat &o_featureMat,
							cv::Mat &o_descriptorMat,
							dsp_times &o_dsptimes);

	/** ------------------------------------------------------------------
	** @internal
	** @brief Pools all sift descriptors
	**
	** @param i_descriptorMat input mat containing all descriptors as columns [dim_descr x (batchsize * ns)]
	** @param batchsize size of one batch (original size of input features to dsp-sift)
	** @param ns number of scales
	** @param o_pooleddescriptorMat pooled descriptors output matrix [dim_descr x batchsize]
	** 
	**/
	void pool_descriptors(cv::Mat& i_descriptorMat, int batchsize, int ns, cv::Mat& o_pooleddescriptorMat);

	/** ------------------------------------------------------------------
	** @internal
	** @brief Sorts the input matrix (4xN) columns according to the given indices 
	**
	** @param i_mat unsorted input matrix 
	** @param i_indices input indices for sort order of columns
	** @param o_mat sorted output matrix
	** 
	**/
	void s_sort_4rowf64_matrixcolsbyindices(cv::Mat &i_mat, cv::Mat &i_indices, cv::Mat &o_mat);

	/** ------------------------------------------------------------------
	** @internal
	** @brief Sorts the input float32 matrix (MxN) columns according to the given indices 
	**
	** @param i_mat unsorted input float32 matrix 
	** @param i_indices input indices for sort order of columns
	** @param o_mat sorted output matrix
	** 
	**/
	void sort_genericf32_matrixcolsbyindices(cv::Mat &i_mat, cv::Mat &i_indices, cv::Mat &o_mat);

	/** ------------------------------------------------------------------
	** @internal
	** @brief Transform scaled keypoints matrix to double array
	**
	** @param i_Dmat input matrix (type double)
	** @param o_pDarray output pointer to double array
	** 
	**/
	void dmat_to_darray(cv::Mat &i_Dmat, double* o_pDarray);

	// Debug methods
	void sorttest();




		/** ------------------------------------------------------------------
	** @internal
	** @brief Normalizes the histogram 
	**
	** @param begin pointer to begin of the array
	** @param end pointer to end of the array
	**
	** @return norm normalization factor (integral value)
	**/
	static inline float normalize_histogram(float* begin, float* end);

	/** ------------------------------------------------------------------
	** @internal
	** @brief Fast real sqrt(x) approximation
	**
	** @param x float value
	**
	** @return sqrt(x) square root of x
	**/
	static inline float vl_fast_resqrt_f(float x);

	/** ------------------------------------------------------------------
	** @internal
	** @brief The function computes a fast approximation of sqrt(x)
	**
	** @param x float value
	**
	** @return sqrt(x) square root of x
	**/
	static inline float vl_fast_sqrt_f(float x);

	void writeMatToFile(cv::Mat& m, const char* filename);

} 
/**************************************** End Namespace dspsift_helperlib ********************************************/
#endif /* DSPSIFTHELPERLIB_H */