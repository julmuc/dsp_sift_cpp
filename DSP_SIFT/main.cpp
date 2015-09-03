/** @internal
 ** @file   main.cpp
 ** @brief   DSP_SIFT based on (c) 2014-2015, Jingming Dong & 
 ** http://synaptic-activity.blogspot.de/2012/02/vlfeat-sift-with-opencv-code.html
 ** @author Julian heuser
 **/

#include <iostream>

extern "C" {
#include <vl/generic.h>
#include <vl/sift.h>
#include <vl/mathop.h>
}
#include <opencv2\core\core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>


#include <math.h>
#include <assert.h>

typedef struct Pair
{
  int k1;
  int k2;
  double score;
} Pair;

typedef struct dspOptions
{
  double sc_min;			// scale sampling lower limit
  double sc_max;			// scale sampling upper limit. Scales are sampled from (sc_min * s, sc_max * s) where s is the detected scale
  int ns;					// number of scales
} dspOptions ;


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
VL_INLINE void
	transpose_descriptor (vl_sift_pix* dst, vl_sift_pix* src)
{
	int const BO = 8;  /* number of orientation bins */
	int const BP = 4;  /* number of spatial bins     */
	int i, j, t;

	for (j=0; j<BP; ++j)
	{
		int jp = BP-1-j;
		for (i=0; i<BP; ++i) 
		{
			int o  = BO*i + BP*BO*j;
			int op = BO*i + BP*BO*jp;
			dst[op] = src[o];
			for (t=1 ; t<BO; ++t)
			{
				dst[BO-t+op] = src[t+o];
			}
		}
	}
}


/** ------------------------------------------------------------------
 ** @internal
 ** @brief Ordering of tuples by increasing scale
 **
 ** @param a tuple.
 ** @param b tuple.
 **
 ** @return @c a[2] < b[2]
 **/

static int korder(void const* a, void const* b)
{
  double x = ((double*)a)[2] - ((double*)b)[2];
  if (x < 0)
	  return -1;
  if (x > 0) 
	  return +1;
  return 0;
}

/** ------------------------------------------------------------------
 ** @internal
 ** @brief Check for sorted keypoints
 **
 ** @param keys keypoint list to check
 ** @param nkeys size of the list.
 **
 ** @return 1 if the keypoints are storted.
 **/

vl_bool check_sorted (double const * keys, vl_size nkeys)
{
  vl_uindex k;
  for (k=0; k + 1<nkeys; ++k)
  {
    if (korder(keys, keys + 4) > 0)
	{
      return VL_FALSE ;
    }
    keys += 4;
  }
  return VL_TRUE;
}



void VLSIFT(IplImage* i_image, vl_uint8* o_DATAdescr, double* o_DATAframes, int* o_nframes);

void VLMATCH(vl_uint8* L1_pt,vl_uint8* L2_pt, int K1, int K2, double thresh, int* nMatches, double* MATCHES);

void VL_DSP_SIFT(IplImage* i_image, dspOptions i_opt, vl_uint8* o_DATAdescr, double* o_DATAframes, int* o_nframes);

int main(int argc, char** argv)
{
	/*
	
	//IplImage* im_rgb = cvLoadImage("C:/Users/Julian/Documents/Visual Studio 2010/Projects/DSP_SIFT/Debug/Lenna.png");
	//Mat i_mat_image = Mat(im_rgb);
	//Mat i_mat_g_image;
	//cv::cvtColor(i_mat_image, i_mat_g_image, CV_RGB2GRAY);

	//if(! i_mat_image.data ) // Check for invalid input
	//{
	//	cout << "Could not open or find the image" << std::endl ;
	//	cv::waitKey(5000);
	//	return -1;
	//}

	//namedWindow( "Display window", WINDOW_AUTOSIZE ); // Create a window for display.
	//imshow( "Display window", i_mat_g_image ); // Show our image inside it.

	////// Convert image to one-dimensional array.
	////double* i_1d_image = new double[i_mat_image.rows*i_mat_image.cols*i_mat_image.channels()];
	////for (int i = 0; i < i_mat_image.rows; ++i) {
	////	for (int j = 0; j < i_mat_image.cols; ++j) {
	////		i_1d_image[j + i_mat_image.cols*i + i_mat_image.cols*i_mat_image.rows*0] = i_mat_image.at<cv::Vec3b>(i, j)[0];
	////		i_1d_image[j + i_mat_image.cols*i + i_mat_image.cols*i_mat_image.rows*1] = i_mat_image.at<cv::Vec3b>(i, j)[1];
	////		i_1d_image[j + i_mat_image.cols*i + i_mat_image.cols*i_mat_image.rows*2] = i_mat_image.at<cv::Vec3b>(i, j)[2];
	////	}
	////}

	//// transform image in cv::Mat to float vector
	//std::vector<float> imgvec;
	//for (int i = 0; i < i_mat_g_image.rows; ++i){
	//	for (int j = 0; j < i_mat_g_image.cols; ++j){
	//		imgvec.push_back(i_mat_g_image.at<unsigned char>(i,j) / 255.0f);                                                                                                                                                                                                        
	//	}
	//}
	//
	*/

	// load template image:
	// needs to be grayscale image
    IplImage* Timage = cvLoadImage("C:/Users/Julian/Documents/Visual Studio 2010/Projects/DSP_SIFT/Debug/Lena.png",0);
 
	// static allocation and zero-initialization of arrays for features and descriptors
    double* TFrames = (double*)calloc(4*10000, sizeof(double));
    vl_uint8* TDescr  = (vl_uint8*)calloc(128*10000, sizeof(vl_uint8));
	
	// stores number of features
    int Tnframes = 0;
	
	// call sift
    VLSIFT(Timage, TDescr, TFrames, &Tnframes);
	
	// reallocate memory block (in case to much space allocated before) 
    TFrames = (double*)realloc(TFrames, 4*sizeof(double)*Tnframes); // = Y X Scale Angle
    TDescr = (vl_uint8*)realloc(TDescr, 128*sizeof(vl_uint8)*Tnframes);
    
	// draw each feature + descriptor region as a circle
    for(int i=0; i<Tnframes; i++){
 
        cvCircle(Timage,											// image
				 cvPoint(TFrames[0+i*4], TFrames[1+i*4]),			// center (x,y)
				 TFrames[2+i*4],									// radius
				 cvScalar(255, 0, 0, 0),							// colour
				 1,													// thickness
				 8,													// linetype
				 0);												// shift
    }
	
	// show image
    cvShowImage("FrameT", Timage);
    cvWaitKey(0);
	return 0;
}

/** ------------------------------------------------------------------
 ** @internal
 ** @brief based on vl_sift.c
 **
 ** @param todo
 ** @param todo
 ** 
 ** @return todo
 **/
void VLSIFT(IplImage* i_image, vl_uint8* o_DATAdescr, double* o_DATAframes, int* o_nframes)
{ 
    //Take IplImage -> convert to SINGLE (float):
    float* frame = (float*)malloc(i_image->height*i_image->width*sizeof(float));
    uchar* Ldata = (uchar*)i_image->imageData;
 
    for(int i=0; i<i_image->height; i++)
	{
        for(int j=0; j<i_image->width; j++)
		{
            frame[j*i_image->height+i*i_image->nChannels] = (float)Ldata[i*i_image->widthStep+j*i_image->nChannels];
		}
	}
 
    // VL SIFT computation:
    vl_sift_pix const *data;
    int M, N;
    data = (vl_sift_pix*)frame;
    M = i_image->height;
    N = i_image->width;
                                                                                                                       
    // VL SIFT PARAMETERS
    int                verbose				= 1; // change to 2 for more verbose..
    int                O					= -1; //Octaves
    int                S					= 3; //Levels
    int                o_min				= 0;
    double             edge_thresh			= -1;  //-1 will use the default (as in matlab)
    double             peak_thresh			= -1;
    double             norm_thresh			= -1;
    double             magnif				= -1;
    double             window_size			= -1;
    double            *ikeys				= 0; //?
    int                nikeys				= -1; //?
    vl_bool            force_orientations	= 0;
    vl_bool            floatDescriptors		= 0;
   
	/* -----------------------------------------------------------------
    *                                                            Do job
    * -------------------------------------------------------------- */
 
	VlSiftFilt	*filt;
	vl_bool		first;
	double		*frames = 0;
	vl_uint8	*descr = 0;
	int			reserved = 0;
	int			i,j,q;

	/* create a filter to process the image */
	filt = vl_sift_new (M, N, O, S, o_min) ;

	if (peak_thresh >= 0) 
		vl_sift_set_peak_thresh(filt, peak_thresh);
	if (edge_thresh >= 0) 
		vl_sift_set_edge_thresh(filt, edge_thresh);
	if (norm_thresh >= 0) 
		vl_sift_set_norm_thresh(filt, norm_thresh);
	if (magnif      >= 0) 
		vl_sift_set_magnif(filt, magnif);
	if (window_size >= 0) 
		vl_sift_set_window_size(filt, window_size);

	if (verbose)
	{
		printf("vl_sift: filter settings:\n") ;
		printf("vl_sift:   octaves      (O)      = %d\n", vl_sift_get_noctaves(filt));
		printf("vl_sift:   levels       (S)      = %d\n", vl_sift_get_nlevels(filt));
		printf("vl_sift:   first octave (o_min)  = %d\n", vl_sift_get_octave_first(filt));
		printf("vl_sift:   edge thresh           = %g\n", vl_sift_get_edge_thresh(filt));
		printf("vl_sift:   peak thresh           = %g\n", vl_sift_get_peak_thresh(filt));
		printf("vl_sift:   norm thresh           = %g\n", vl_sift_get_norm_thresh(filt));
		printf("vl_sift:   window size           = %g\n", vl_sift_get_window_size(filt));
		printf("vl_sift:   float descriptor      = %d\n", floatDescriptors);
		printf((nikeys >= 0) ? "vl_sift: will source frames? yes (%d read)\n" : "vl_sift: will source frames? no\n", nikeys);
		printf("vl_sift: will force orientations? %s\n", force_orientations ? "yes" : "no") ;
	} 
	/* ...............................................................
	*                                             Process each octave
	* ............................................................ */

	i     = 0;
	first = 1;
	while (1) 
	{
		int                   err;
		VlSiftKeypoint const *keys  = 0;
		int                   nkeys = 0;

		if (verbose)
		{
			printf ("vl_sift: processing octave %d\n", vl_sift_get_octave_index(filt));
		}

		/* Calculate the GSS for the next octave .................... */
		if (first) 
		{
			err   = vl_sift_process_first_octave(filt, data);
			first = 0;
		} 
		else 
		{
			err   = vl_sift_process_next_octave(filt);
		}

		if (err) 
			break;

		if (verbose > 1)
		{
			printf("vl_sift: GSS octave %d computed\n", vl_sift_get_octave_index(filt));
		}

		/* Run detector ............................................. */
		if (nikeys < 0)
		{
			vl_sift_detect(filt);
			keys  = vl_sift_get_keypoints(filt);
			nkeys = vl_sift_get_nkeypoints(filt);
			i     = 0;

			if (verbose > 1) 
			{
				printf ("vl_sift: detected %d (unoriented) keypoints\n", nkeys);
			}

		}
		else
		{
			nkeys = nikeys;
		}
		
		/* For each keypoint ........................................ */
		for (; i<nkeys; ++i)
		{
			double                angles[4];
			int                   nangles;
			VlSiftKeypoint        ik;
			VlSiftKeypoint const *k;

			/* Obtain keypoint orientations ........................... */
			if (nikeys >= 0)
			{
				vl_sift_keypoint_init(filt, &ik, ikeys [4 * i + 1] - 1, ikeys [4 * i + 0] - 1, ikeys [4 * i + 2]);

				if (ik.o != vl_sift_get_octave_index(filt))
					break;

				k = &ik;

				/* optionally compute orientations too */
				if (force_orientations) 
				{
					nangles = vl_sift_calc_keypoint_orientations(filt, angles, k);
				} 
				else
				{
					angles[0] = VL_PI/2 - ikeys[4*i + 3];
					nangles   = 1;
				}
			} 
			else
			{
				k = keys + i;
				nangles = vl_sift_calc_keypoint_orientations(filt, angles, k);
			}

			/* For each orientation ................................... */
			for(q=0; q<nangles; ++q) 
			{
				vl_sift_pix  buf[128];
				vl_sift_pix rbuf[128];

				/* compute descriptor (if necessary) */
				vl_sift_calc_keypoint_descriptor(filt, buf, k, angles [q]);
				transpose_descriptor(rbuf, buf);

				/* make enough room for all these keypoints and more */
				if (reserved < (*o_nframes) + 1)
				{
					reserved  += 2*nkeys;
					frames = (double*)realloc(frames, 4*sizeof(double)*reserved);
					descr  = (vl_uint8*)realloc(descr,  128*sizeof(vl_uint8)*reserved);
				}

				/* Save back with MATLAB conventions. Notice that the input
				* image was the transpose of the actual image. */
				frames[4*(*o_nframes) + 0] = k->y + 1;
				frames[4*(*o_nframes) + 1] = k->x + 1;
				frames[4*(*o_nframes) + 2] = k->sigma;
				frames[4*(*o_nframes) + 3] = VL_PI/2 - angles[q];

				for(j=0; j<128; ++j)
				{
					float x = 512.0F * rbuf[j];
					x = (x < 255.0F) ? x : 255.0F;
					descr[128*(*o_nframes) + j] = (vl_uint8)x;
				}
				++(*o_nframes);
			} /* next orientation */
		} /* next keypoint */
	} /* next octave */
 


	if (verbose)
		printf("vl_sift: found %d keypoints\n", (*o_nframes));

	// save variables:
	memcpy(o_DATAframes, frames, 4*(*o_nframes)*sizeof(double));
	memcpy(o_DATAdescr, descr, 128*(*o_nframes)*sizeof(vl_uint8));

	/* cleanup */
	vl_sift_delete(filt);
	
	/* end: do job */
	return;
}

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
void VLMATCH(vl_uint8* L1_pt,vl_uint8* L2_pt, int K1, int K2, double thresh, int* nMatches, double* MATCHES )
{ 
    //Match descriptors!
    int ND = 128;
 
    Pair* pairs_begin = (Pair*) malloc(sizeof(Pair) * (K1+K2)) ;
    Pair* pairs_iterator = pairs_begin ; 
 
    int k1, k2;                                                       
 
    const int maxval = 0x7fffffff ;                        
    for(k1=0; k1<K1; ++k1, L1_pt += ND)
	{                    
        int best = maxval;                                    
        int second_best = maxval;                             
        int bestk = -1 ;                                                 
 
        /* For each point P2[k2] in the second image... */              
        for(k2=0; k2<K2; ++k2, L2_pt  += ND)
		{                      
            int bin;
            int acc = 0;
            for(bin=0; bin<ND; ++bin)
			{                              
                int delta = ((int) L1_pt[bin]) - ((int) L2_pt[bin]);                             
                acc  += delta*delta;                                         
            }
            
			/* Filter the best and second best matching point. */
            if(acc < best)
			{          
                second_best = best;                                         
                best = acc; 
                bestk = k2;
            }
			else if(acc < second_best)
			{                                 
                second_best = acc;
            }
        }
        L2_pt -= ND*K2;

		/* Lowe's method: accept the match only if unique. */
		if(thresh*(float)best < (float)second_best && bestk!=-1)
		{
			pairs_iterator->k1 = k1;
			pairs_iterator->k2 = bestk;
			pairs_iterator->score = best;
			pairs_iterator++;
			(*nMatches)++; //todo: ?rly
		}
	}
    Pair* pairs_end = pairs_iterator;
 
    //double* M_pt = (double*)calloc((pairs_end-pairs_begin)*2,sizeof(double));
    double* M_pt = (double*)calloc((*nMatches)*2,sizeof(double));
    //double* M_start = M_pt;
 
    for(pairs_iterator = pairs_begin; pairs_iterator<pairs_end; ++pairs_iterator) {
            *M_pt++ = pairs_iterator->k1 + 1;
            *M_pt++ = pairs_iterator->k2 + 1;
    }
 
    M_pt -= (*nMatches)*2;
    memcpy(MATCHES, M_pt, (*nMatches)*2*sizeof(double));
    free(pairs_begin);
    free(M_pt);
 
    return; 
}

void VL_DSP_SIFT(IplImage* i_image, dspOptions i_opt, vl_uint8* o_DATAdescr, double* o_DATAframes, int* o_nframes)
{	

	//----------------------------------------------- Detection -----------------------------------------------------//
	// static allocation and zero-initialization of arrays for features(frames) and descriptors
    double* siftFrames = (double*)calloc(4*10000, sizeof(double));
    vl_uint8* siftDescr  = (vl_uint8*)calloc(128*10000, sizeof(vl_uint8));
	
	// stores number of features(frames)
    int nframes = 0;
	
	// call sift
    VLSIFT(i_image, siftDescr, siftFrames, &nframes);
	
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