// Copyright (C) 2007-12 Andrea Vedaldi and Brian Fulkerson.
// All rights reserved.

/** @internal
 ** @file   vlfeat_helperlib.h
 ** @brief   vl_sift based on (c) vlfeat.org &
 ** http://synaptic-activity.blogspot.de/2012/02/vlfeat-sift-with-opencv-code.html
 ** @author Julian heuser
 **/

#ifndef VLFEATHELPERLIB_H
#define VLFEATHELPERLIB_H

//#define NDEBUG	// define for no debug messages

/************************************************** Includes *********************************************************/
extern "C" {
#include <vl/generic.h>
#include <vl/sift.h>
#include <vl/mathop.h>
}

#include <opencv2\core\core.hpp>
#include "debug_helper.h"

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
	/*
	VL_SIFT() accepts the following options:

	Octaves:: maximum possible
	Set the number of octave of the DoG scale space.

	Levels:: 3
	Set the number of levels per octave of the DoG scale space.

	FirstOctave:: 0
	Set the index of the first octave of the DoG scale space.

	PeakThresh:: 0
	Set the peak selection threshold.

	EdgeThresh:: 10
	Set the non-edge selection threshold.

	NormThresh:: -inf
	Set the minimum l2-norm of the descriptors before
	normalization. Descriptors below the threshold are set to zero.

	Magnif:: 3
	Set the descriptor magnification factor. The scale of the
	keypoint is multiplied by this factor to obtain the width (in
	pixels) of the spatial bins. For instance, if there are there
	are 4 spatial bins along each spatial direction, the
	``side'' of the descriptor is approximatively 4 * MAGNIF.

	WindowSize:: 2
	Set the variance of the Gaussian window that determines the
	descriptor support. It is expressend in units of spatial
	bins.

	Frames::
	If specified, set the frames to use (bypass the detector). If
	frames are not passed in order of increasing scale, they are
	re-orderded.

	Orientations::
	If specified, compute the orientations of the frames overriding
	the orientation specified by the 'Frames' option.

	Verbose::
	If specfified, be verbose (may be repeated to increase the
	verbosity level).
	*/
	typedef struct vl_sift_options
	{	
		vl_sift_options(): verbose(1),
							O(-1), 
							S(3),
							o_min(0),
							edge_thresh(-1),
							peak_thresh(-1),
							norm_thresh(-1),
							magnif(-1),
							window_size(-1),
							ikeys(0),
							nikeys(-1),
							force_orientations(0),
							floatDescriptors(0),
							ikeys_provided(0){ }   // default Constructor

		int                verbose; 
		int                O; //Octaves
		int                S; //Levels
		int                o_min;
		double             edge_thresh;  //-1 will use the default (as in matlab)
		double             peak_thresh;
		double             norm_thresh;
		double             magnif;
		double             window_size;
		double            *ikeys;			//
		int                nikeys; 
		vl_bool            force_orientations;
		vl_bool            floatDescriptors;
		bool			   ikeys_provided;
	} vl_sift_options;


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
	void vlsift(IplImage* i_image, void* o_DATAdescr, double* o_DATAframes, int* o_nframes, vl_sift_options opts);

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
	void vlmatch(vl_uint8* L1_pt,vl_uint8* L2_pt, int K1, int K2, double thresh, int* nMatches, double* MATCHES);

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

	/** ------------------------------------------------------------------
	** @internal
	** @brief Normalizes the histogram 
	**
	** @param begin pointer to begin of the array
	** @param end pointer to end of the array
	**
	** @return norm normalization factor (integral value)
	**/
	float normalize_histogram(float* begin, float* end);

	/** ------------------------------------------------------------------
	** @internal
	** @brief Fast real sqrt(x) approximation
	**
	** @param x float value
	**
	** @return sqrt(x) square root of x
	**/
	float vl_fast_resqrt_f(float x);

	/** ------------------------------------------------------------------
	** @internal
	** @brief The function computes a fast approximation of sqrt(x)
	**
	** @param x float value
	**
	** @return sqrt(x) square root of x
	**/
	float vl_fast_sqrt_f (float x);

} 
/***************************************** End Namespace vlfeat_helperlib ********************************************/
#endif /* VLFEATHELPERLIB_H */