/** @internal
 ** @file   dsp_sift.h
 ** @brief   dsp_sift based on (c) vlfeat.org & J. Dong
 ** http://vision.ucla.edu/~jingming/proj/dsp/
 ** @author Julian heuser
 **/

#ifndef DSPSIFTHELPERLIB_H
#define DSPSIFTHELPERLIB_H

/************************************************** Includes *********************************************************/
#include <opencv2\core\core.hpp>
#include "vlfeat_helperlib.h"

/************************************************** Structures *******************************************************/


/**************************************** Begin Namespace dspsift_helperlib ******************************************/
namespace dspsift_helperlib
{

	typedef struct dspOptions
	{
		double sc_min;			// scale sampling lower limit
		double sc_max;			// scale sampling upper limit. Scales are sampled from (sc_min * s, sc_max * s) where s is the detected scale
		int ns;					// number of scales
	} dspOptions ;

	/******************************************** Function Declarations **********************************************/

	/** ------------------------------------------------------------------
	** @internal
	** @brief based on dsp_sift by J.Dong
	**
	** @param todo
	** @param todo
	** 
	** @return todo
	**/
	void DSP_SIFT(IplImage* i_image, dspOptions i_opt, vl_uint8* o_DATAdescr, double* o_DATAframes, int* o_nframes);

	/** ------------------------------------------------------------------
	** @internal
	** @brief Sample scales around detection
	**
	** @param i_DATAframes input features 
	** @param i_nframes total number of features
	** @param i_opt scale options
	** @param const cv::Mat &o_sampledfeatureMat const as matrix dimension do not change, (data values surely do!)
	** 
	** @return todo
	**/
	void sampleScales(double* i_DATAframes, int* i_nframes, dspOptions i_opt, const cv::Mat &o_sampledfeatureMat);


} 
/**************************************** End Namespace dspsift_helperlib ********************************************/
#endif /* DSPSIFTHELPERLIB_H */