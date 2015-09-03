/** @internal
 ** @file   vlfeat_helperlib.h
 ** @brief   vl_sift based on (c) vlfeat.org &
 ** http://synaptic-activity.blogspot.de/2012/02/vlfeat-sift-with-opencv-code.html
 ** @author Julian heuser
 **/

#ifndef VLFEATHELPERLIB_H
#define VLFEATHELPERLIB_H

/************************************************** Includes *********************************************************/
extern "C" {
#include <vl/generic.h>
#include <vl/sift.h>
#include <vl/mathop.h>
}

#include <opencv2\core\core.hpp>

/************************************************** Structures *******************************************************/

typedef struct Pair
{
  int k1;
  int k2;
  double score;
} Pair;

/**************************************** Begin Namespace vlfeat_helperlib *******************************************/
namespace vlfeat_helperlib
{

	/******************************************** Function Declarations **********************************************/

	/** ------------------------------------------------------------------
	** @internal
	** @brief Transpose desriptor
	**
	** @param dst destination buffer.
	** @param src source buffer.
	**
	** The function writes to @a dst the transpose of the SIFT descriptor
	** @a src. The tranpsose is defined as the descriptor that one
	** obtains from computing the normal descriptor on the transposed
	** image.
	**/
	void transpose_descriptor(vl_sift_pix* dst, vl_sift_pix* src);

	/** ------------------------------------------------------------------
	** @internal
	** @brief based on vl_sift.c
	**
	** @param todo
	** @param todo
	** 
	** @return todo
	**/
	void VLSIFT(IplImage* i_image, vl_uint8* o_DATAdescr, double* o_DATAframes, int* o_nframes);

	/** ------------------------------------------------------------------
	** @internal
	** @brief based on vl_siftmatch.c/vl_ubcmatch.c
	**
	** @param L1_pt is a descriptor array of type vl_uint8
	** @param L2_pt is a descriptor array of type vl_uint8
	** @param K1 is the size if L1_pt
	** @param K2 is the size if L2_pt
	** 
	** @return todo
	**/
	void VLMATCH(vl_uint8* L1_pt,vl_uint8* L2_pt, int K1, int K2, double thresh, int* nMatches, double* MATCHES);

	/** ------------------------------------------------------------------
	** @internal
	** @brief Ordering of tuples by increasing scale
	**
	** @param a tuple.
	** @param b tuple.
	**
	** @return @c a[2] < b[2]
	**/
	static int korder(void const* a, void const* b);

	/** ------------------------------------------------------------------
	** @internal
	** @brief Check for sorted keypoints
	**
	** @param keys keypoint list to check
	** @param nkeys size of the list.
	**
	** @return 1 if the keypoints are storted.
	**/
	vl_bool check_sorted(double const * keys, vl_size nkeys);


} 
/***************************************** End Namespace vlfeat_helperlib ********************************************/
#endif /* VLFEATHELPERLIB_H */