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
#include "vlfeat_helperlib.h"
#include "debug_helper.h"
#include <platformstl\performance\performance_counter.hpp>

/************************************************** Structures *******************************************************/


/**************************************** Begin Namespace dspsift_helperlib ******************************************/
namespace dspsift_helperlib
{

	typedef struct dspOptions
	{
		double sc_min;			// scale sampling lower limit
		double sc_max;			// scale sampling upper limit. Scales are sampled from (sc_min * s, sc_max * s) where s is the detected scale
		int ns;					// number of scales
		vlfeat_helperlib::vl_sift_options vlsift_opt;
	} dspOptions ;

	/******************************************** Function Declarations **********************************************/

	/** ------------------------------------------------------------------
	** @internal
	** @brief based on dsp_sift by J.Dong
	**
	** @param i_image input greyscale image
	** @param i_opt input option for sampling different scales
	** @param o_DATAdescr vl_uint8* pointer to array of output descriptors
	** @param o_DATAframes double* pointer to array of output frames(sift features)
	** @param o_nframes int* to the number of frames
	** 
	**/
	void dsp_sift(IplImage* i_image, dspOptions i_opt, vl_uint8* o_DATAdescr, double* o_DATAframes, int* o_nframes);

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
	void samplescales(double* i_DATAframes, int* i_nframes, dspOptions i_opt, cv::Mat &o_sampledfeatureMat);

	/** ------------------------------------------------------------------
	** @internal
	** @brief Sorts the input matrix columns according to the given indices 
	**
	** @param i_mat unsorted input matrix 
	** @param i_indices input indices for sort order of columns
	** @param o_mat sorted output matrix
	** 
	**/
	void sortmatrixcolsbyindices(cv::Mat &i_mat, cv::Mat &i_indices, cv::Mat &o_mat);

	// Debug methods
	void sorttest();


} 
/**************************************** End Namespace dspsift_helperlib ********************************************/
#endif /* DSPSIFTHELPERLIB_H */