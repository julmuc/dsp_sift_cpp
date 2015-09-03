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


	//--------------------------------- Compute un-normalized SIFT at each scales -----------------------------------//


	// todo


	//------------------------------------------- Aggregate and normalize -------------------------------------------//


	// todo



	return;
}