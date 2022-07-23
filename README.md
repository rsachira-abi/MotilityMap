# MotilityMap

A full description of the MotilityMap method is available in the journal article: [High-Resolution Spatiotemporal Quantification of Intestinal Motility With Free-Form Deformation](https://ieeexplore.ieee.org/document/9652118). If you are using this method, please site the original paper:

```S. Kuruppu et al., "High-Resolution Spatiotemporal Quantification of Intestinal Motility With Free-Form Deformation," in IEEE Transactions on Biomedical Engineering, vol. 69, no. 6, pp. 2077-2086, June 2022, doi: 10.1109/TBME.2021.3135855.```

MotilityMap implementation in MATLAB. Tested on MATLAB version 2021b. The polygon region of interest selection tool (used to specify the boundary) may not be availabe in MATLAB versions R2019b and below.

## Downloading sample data

Test data can be downloaded through the following link:
https://drive.google.com/drive/folders/11J3FOrGdclwJXfVqu09DGRxzqxEI3u_M?usp=sharing

* Test Data RAW 
    - A raw video sequence. Each \*.raw file contains a video frame in RAW format.
    - Frames can be read using ImageRegistration.ImportRaw() function.
* sample_data.avi
    - An AVI video sequence.
* sample_timevec.mat
    - Time vector for sample_data.avi
    - Contains the time (in seconds) of each frame.

## Running the unit tests
Unit tests are found in the "Unit Tests" folder. Test scripts can be run with MATLAB.

* test_ImageRegistration = Tests to validate the integer and subpixel shift detection accuracy of the image registration algorithm.
* test_ImageRegistrationCUDA = Same as above, but for the CUDA implementation of the image registration algorithm.
* test_StrainMeasurement = Tests to validate the free-form deformation strain measurement algorithm with and without optimization.

## Apply the method to sample AVI video sequence

1. Run `GenDisplacementFieldsAVI.m` to compute the displacement fields for the video sequence.
2. Run `AnalyzeVideoSegment.m`, select the project directory created by step 1, and follow the instructions to specify the boundary and calculate strain fields.
3. Run `GenStrainMap.m`, select the project directory (as in step 2), and select `sample_timevec.mat` as the timevec, to generate the strain map.
4. The strain map can be differentiated along time with the MATLAB `gradient()` function to obtain the strain-rate map.

If NVIDIA CUDA is not available, set `USE_CUDA` to `false` in the above MATLAB scripts.

## Apply the method to sample RAW video sequence

Same as previous section, but use `GenDisplacementFieldsRAW.m` to compute the displacement fields from RAW frames.
 
 
