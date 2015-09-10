# dsp_sift_cpp

The project contains an implementation in c++ of the domain-size pooling sift algorithm by J. Dong.

###How-to Setup for Building

[OpenCV 2.4.11 Setup VS2010](http://docs.opencv.org/doc/tutorials/introduction/windows_visual_studio_Opencv/windows_visual_studio_Opencv.html)

[VLFeat Setup VS2010](http://www.vlfeat.org/vsexpress.html)

[STLSoft](http://www.stlsoft.org/index.html) [temporarily used for timing]



###Content
![Flowchart](https://cdn.pbrd.co/images/sNj45Vx.png)

* **VLFeat Helperlib**
Contains functions to implement the [vl_sift](http://www.vlfeat.org/api/sift.html) algorithm

 `vlsift`   
 `vlmatch`   
 `transpose_descriptor`   
 `korder`   
 `check_sorted`   

* **DSPSift Helperlib**
Contains functions to implement the [dsp sift](http://vision.ucla.edu/~jingming/proj/dsp/) algorithm

 `dsp_sift`   
 `sample_scales`   
 `get_all_descriptors`   
 `pool_descriptors`   
 `s_sort_4rowf64_matrixcolsbyindices`   
 `sort_genericf32_matrixcolsbyindices`   
 `dmat_to_darray`   
 `get_final_output_features`   
 `get_normalized_descriptors`   
 `normalize_histogram`

###References
* [J. Dong and S. Soatto. Domain-Size Pooling in Local Descriptors: DSP-SIFT. In _Proceedings of the IEEE Conference on Computer Vision and Pattern Recognition_ (CVPR), 2015](http://vision.ucla.edu/~jingming/proj/dsp/)

* [VLFeat](http://www.vlfeat.org/index.html)

* [OpenCV](http://opencv.org/)

* [STLSoft](http://www.stlsoft.org/index.html) [temporarily used for timing]