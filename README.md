# VANUC: Voxel-wise Anatomical-region-based Non-Uniformity Correction for partial volume effects

<img src="/fig1b.jpg" width="100%">  

[Japanese Page](/README_JP.md)  

# Overview
<img src="/VANUClogoM.jpg" width="15%">  

**VANUC** is a program for performing partial volume effect correction (PVEC) of PET and SPECT images using VANUC method[^1].  
[^1]:[Akira Arai, Yuriko Kagaya, Kentaro Takanami, Kei Takase. A novel partial volume effect correction method which allows for heterogeneous distribution: The potential advantages in the white matter activity estimation on FDG-PET. J Nucl Med. 2016;57(supplement 2):1924](https://jnm.snmjournals.org/content/57/supplement_2/1924)  

* With a simple push of a button, VANUC reads image data, aligns it, corrects partial volume effects, and outputs standard brain coordinate images for statistical analysis.  
* By inputting Hoffmann brain phantom PET images, the program can automatically calculate the parameters necessary for performing PVEC.  

This program runs on MATLAB and requires SPM to be installed.    
(Developer: Akira Arai)  

# Table of Contents
1. [Preparation before Use](#1-preparation-before-use)
    1. [Usage Environment](#1-1-usage-environment)
    1. [Download VANUC](#1-2-download-vanuc)
    1. [Path Settings for VANUC](#1-3-path-settings-for-vanuc)
    1. [Image Preparation](#1-4-image-preparation)
1. [How to Use](#2-how-to-use)
    1. [Execution of VANUC Program](#2-1-execution-of-vanuc-program)
    1. [Output Data Lists](#2-2-output-data-lists)
1. [Principle](#3-principle)
    1. [Program Processing Flow](#3-1-program-processing-flow)
    1. [Theory of Partial Volume Effect Correction Methods](#3-2-theory-of-partial-volume-effect-correction-methods)
    1. [Results of Previous Studies](#3-3-results-of-previous-studies)
1. [Licenses](#4-licenses)
1. [References](#5-references)

# 1. Preparation before Use
## 1-1. Usage Environment
This software runs on the numerical calculation software MATLAB.  
* Click [here](https://jp.mathworks.com/products/matlab.html) to purchase MATLAB (external site)  

To use SPM12 functions, SPM12 download and path setting are required.  
* Click [here](https://www.fil.ion.ucl.ac.uk/spm/software/spm12/) to download SPM12 (external site)  
*  After downloading, please set the path settings for SPM12 in MATLAB.  
   1. Launch MATLAB.  
   1. Open "ENVIRONMENT" > "Set Path" on the toolbar.  
   1. Click "Add Folder..." and select "spm" folder, then click "Save".  
   1. Close the Path Setting window, type "spm" in the MATLAB command window, and confirm that SPM12 starts up.  

## 1-2. Download VANUC
1. Click "Code" on [this page](../../) and click Download ZIP to download all the files in the repository. Unzip the downloaded zip file.  
1. Save the entire "vanuc" folder. To avoid errors, it is preferable to save the files in a location such as directly under the C: drive.  

## 1-3. Path Settings for VANUC
1. Launch MATLAB.  
1. Open "ENVIRONMENT" > "Set Path" on the toolbar.  
1. Click "Add Folder...". and select "vanuc" folder, then click "Save".  
1. Close the Path Setting window, type "vanuc" in the MATLAB command window, and confirm that vanuc starts up.  

## 1-4. Image Preparation
* For partial volume effect correction process, etc.  
   * Brain PET or SPECT and brain 3D-T1-weighted images are required.  
   * DICOM or NIfTI formats are supported.  
      * Conversion to SUV on PET images is only possible with DICOM format images.  
* For phantom images
   * Only DICOM format is supported at this time.  
* For other formats or images that cannot be read by this software  
   * Please use other software to convert to NIfTI format.   

# 2. How to Use
## 2-1. Execution of VANUC Program
* Enter "vanuc" in the MATLAB command window to open the GUI of VANUC.  
* Processed images and other data will be output in the "current folder" of MATLAB, so change the output destination accordingly.

### Partial Volume Effect Correction Tool
* Corrects partial volume effects using the [VANUC](#vanuc), [Muller-Gartner](#muller-gartner) and [RBV](#rbv) methods, and outputs the corrected data and images.   
1. Click **"Start"** button.  
1. **Select PET or SPECT image.** If you want to convert PET to SUV values, load "DICOM" data.  
   * If "DICOM" is selected:  
      1. Select the folder under that the DICOM files are stored.  
      1. Select the target image from the list of DICOM data in the folder and subfolders.  

   * If "NIfTI" is selected:
      1. Select the NIfTI file.
      1. Select whether or not cropping (cutting off the parts other than the head) should be performed. Normally, select "OK", but if cropping has already been performed, select "Skip".  

1. **Align the orientation of the image.** If the image directions and the H (up), F (down), A (front), P (back), R (right), and L (left) are matched, click "Confirm directions". If not, click "Rotate/Flip" and follow the directions in the order shown: up, forward, left and right.  
1. **Cropping** The cropped image will be displayed in color.  
   * If cropping does not work:   
      1. Stop the program.  
      1. If necessary, use other software to crop the NIfTI file (in the case of PET, the file begins with SUV_) that was output in step 3.  
      1. Run the process again, select the created image file as PET or SPECT image, and Skip Cropping.  

1. **Set PET parameters.** Select whether to preset values (Preset) or to estimate from the image (Image-based). Normally, "Preset" is recommended.  
   * If "Preset" is selected:  
      * The parameter values that were used the last time (default values for the first time) are first displayed. If you need to change the parameters, you can either select from the protocols saved in advance (Other Protocol) or directly input them (Direct Imput). The value entered by "Direct Imput" can be saved as a new protocol and selected the next time. If you want to determine parameters from phantom data, [see below](#phantom-data-analysis-tool).  
      * Explanation of setting values:  
         <details>
         <summary>FWHM of PSF: </summary>

         Full width at half maximum (mm) of the point spread function of the PET image (x: left to right, y: front to back, z: up and down).  
         </details>
         <details>
         <summary>Distance from center of local region: </summary>
         Distance (in voxels) from the center voxel to the edge of the kernel used in the VANUC method. Normally, the default value is sufficient, but if the resolution is low and there is a lot of noise, a larger value should be set. Note that larger values increase processing time.  
         </details>

         <details>
         <summary>Regularization parameter: </summary>
         Regularization parameter used in the VANUC method. For clinical PET and SPECT images, the value is between 1 and 100. If the noise is large, a larger value should be used. Strictly determined by L-curve analysis, etc.  
         </details>

1. Select **a 3D T1-weighted image** , orient the image, and cropping is performed. （same as above for PET or SPECT).  
1. Partial volume effect correction begins and takes 10 to 30 minutes. (During the process, if the TPM file (tissue probability map) is not found in the designated location in the SPM folder, you will be asked for the file location.)  

### Phantom Data Analysis Tool
* Tool to determine PET parameters needed for partial volume effect correction  
* [See below](#parameter-determination-in-vanuc) for the principle.  
1. Click on "Phantom."  
1. Select the folder where the phantom images (DICOM) are stored and select the target phantom image from the list.  
1. Orient the phantom image. (Same as above)  
1. Select a digital phantom; the data for the Hoffman 3D brain phantom is pre-populated.  
1. Analysis takes a few minutes.  
1. The FWHM of PSF value is obtained, and the other parameters and the protocol name are entered and saved as a new protocol. If you do not want to save it, click "Cancel" to exit. Once saved, the parameter settings of the saved protocol can be used for the next partial volume effect correction.  

### Other Processing Tools
* dcm>nii: DICOM to NIfTI conversion (including conversion to SUV values for PET)  
* Orientation: Transformation of image orientation (e.g., rotation, left/right flip, etc.)  
* Crop: Cutting the image except for the head (Orientation is processed first)  

## 2-2. Output Data Lists
### Partial volume effect correction output data
<details>
<summary>List in the current folder (main output images such as images in process and partial volume effect corrected images):</summary>

|Data Name|Description|
|:---|:---|
|PT~.nii or NM~.nii|Original PET or SPECT image|
|SUV_PT~.nii|Converted image to SUV values (only for PET)|
|trim_PET.nii|PET (or SPECT) image after cropping|
|MR~.nii|Original MRI image|
|T1.nii|MRI image after cropping|
|coT1.nii|isovoxel MRI image|
|coT1_seg8.mat|Segmentation data by SPM|
|c1~6coT1.nii|Maps of gray matter, white matter, cerebrospinal fluid, skull, soft tissue, and air|
|y_coT1.nii|transformation field from MRI coordinates to MNI standard brain coordinates|
|c1acoT1.nii|cerebellar gray matter map|
|c1bcoT1.nii|map of gray matter other than cerebellum|
|rcoT1.nii|MRI image coregistered to PET (or SPECT)|
|rtrim_PET.nii|PET (or SPECT) image coregistered to MRI |
|PVCvanuc.nii|Partial volume effect corrected image (7 regions) using VANUC method|
|PVCmg.nii|Partial volume effect corrected image (7 areas) using Muller-Gartner method|
|PVCrbv.nii|Partial volume effect corrected image (7 areas) using RBV method|
|rPVC~.nii (3 types)|Partial volume effect corrected images (7 areas) coregistered to MRI|
|PVCC~.nii (3 types)|Partial volume effect-corrected images composited from all 7 regions|
|PVRC~.nii (3 types)|Partial volume effect "reduced" images (FWHM 2.5mm) composited from all 7 regions|
</details>

<details>
<summary>List in MNIspace folder (images converted to MNI standard brain coordinates, used for statistical image analysis by SPM, etc.):</summary>

|Data Name|Description|
|:---|:---|
|wcoT1.nii|Anatomical standardized MRI image|
|wrtrim_PET.nii|Anatomical standardized PET (or SPECT) image (SUV values for PET, original values for SPECT)|
|wrPVC~.nii (3 methods x 3 regions)|Anatomically standardized partial volume effect corrected image (cerebellar gray matter, other gray matter, white matter)|
|wrPVC~\_cbl.nii (3 methods x 3 regions)|Anatomical standardized partial volume effect corrected images normalized by cerebellar gray matter values|
|wrPVC~\_gm.nii (3 methods x 3 regions)|Anatomical standardized partial volume effect corrected images normalized by gray matter values|
|wrPVC~\_wb.nii (3 methods x 3 regions)|Anatomical standardized partial volume effect corrected images normalized by whole brain values|
</details>

<details>
<summary>List in temp folder (raw data of partial volume effect correction):</summary>

|Data Name|Description|
|:---|:---|
|Ppsf.mat|Parameters for point spread function|
|Ploc.mat|Parameters for kernel size in VANUC method|
|Plam.mat|Parameters for regularization in VANUC method|
|G.mat|PET (or SPECT) voxel values|
|M.mat|Tissue map data of 7 regions obtained by MRI segmentation|
|V.mat|Region spread function (RSF) of 7 regions|
|Rgtm.mat|Correction values for partial volume effect in 7 regions by GTM method|
|Rols.mat|Expected value component of 7 regions by VANUC method|
|Rvanuc.mat|Partial volume effect corrected image data of 7 regions by VANUC method|
|Rmg.mat|Partial volume effect corrected image data of 7 regions by Muller-Gartner method|
|Rrbv.mat|Partial volume effect corrected image data of 7 regions by RBV method|
|Voxel.mat|Voxel size parameters|
</details>

### Phantom data analysis output data
(in preparation)  

# 3. Principle
## 3-1. Program Processing Flow
The following is the processing flow of the partial volume effect correction tool.
1. Input of original PET or SPECT image, 3D T1-weighted image, and several parameters.  
1. Segmentation of the 3D T1-weighted image  
1. Anatomical standardization of the 3D T1-weighted image to MNI standard brain coordinates  
1. Inverse transformation of the cerebellar mask to individual brain coordinates  
1. Coregistration of each obtained tissue probability maps (cerebellar gray matter, other gray matter, white matter, cerebrospinal fluid, skull, soft tissue, air) to the PET or SPECT image  
1. Partial volume effect correction  
1. Output of various images ([see above](#partial-volume-effect-correction-output-data) for output images)  

## 3-2. Theory of Partial Volume Effect Correction Methods
Most partial volume effect correction methods are based on a convolutional model as follows. That is, the PET or SPECT image $G$ is represented by the following equation, where the true radioactivity concentration distribution is $\rho$ and the PET or SPECT point spread function is $psf$:  

$$G = \rho \otimes psf$$  

Although we wish to obtain the true radioactivity concentration distribution $\rho$ from the image $G$ based on this model, in reality, $G$ is degraded by discrepancies from the model and statistical noise, so it is difficult to simply estimate $\rho$ by deconvolution. Therefore, the method introduced below compensates for the information on the boundaries of tissue regions lost due to image degradation by using tissue region maps $M_k, (i=1,2,...,n)$ obtained from morphological images such as MRI. Each tissue map $M_k$ represents the volume fraction of the tissue component in a voxel, such as that obtained by SPM segmentation. In addition, as is often the case, the VANUC program models the $psf$ by a 3D Gaussian function.  

### GTM
The GTM (geometry transfer matrix) method is a technique reported by Rousset et al[^2]. The GTM method assumes that the radioactivity concentration within each tissue region is constant.  
Based on this assumption, a transformation matrix $\mathbf{GTM}$ is obtained from the true radioactivity concentration within each region to the average pixel value within each region on the image. Using this transformation matrix, the true radioactivity concentration $\mathbf{\rho}$ is estimated from the average pixel value $\mathbf{g}$ obtained from the actual image by inverse transformation:  
[^2]:[O G Rousset, Y Ma, A C Evans. Correction for partial volume effects in PET: principle and validation. J Nucl Med. 1998 May;39(5):904-11](https://pubmed.ncbi.nlm.nih.gov/9591599/)  

$$\mathbf{\hat{\rho_{GTM}}} = \mathbf{GTM^{-1}} \mathbf{g}, GTM_{i,j} = \dfrac{\displaystyle \sum M_i RSF_j}{\displaystyle \sum M_i}$$  

Where $RSF$ is called the region spreading function:  

$$RSF_k = M_k \otimes psf$$  

### Muller-Gartner
The Muller-Gartner method is a method reported by Muller-Gartner[^3]. The Muller-Gartner method also assumes that the radioactivity concentration in the tissue is constant.  
Based on this assumption, the spill-in from other surrounding tissues is first removed, and the remaining pixel values are divided by the region spread function to obtain the true radioactivity distribution of the tissue.  

$$\hat{\rho_{MG}} (k) = \dfrac{G - \displaystyle \sum_{i \neq k} C_i RSF_i}{RSF_k}$$  

Where the radioactivity concentration $C_i$ of other surrounding tissues must be given, and the VANUC program uses the value $\hat{\rho_{GTM}}(i) obtained by the GTM method.  
[^3]:[H W Muller-Gartner, J M Links, J L Prince, R N Bryan, E McVeigh, J P Leal, C Davatzikos, J J Frost. Measurement of radiotracer concentration in brain gray matter using positron emission tomography: MRI-based correction for partial volume effects. J Cereb Blood Flow Metab. 1992 Jul;12(4):571-83](https://pubmed.ncbi.nlm.nih.gov/1618936/)  

### RBV
The RBV (region-based voxel-wise correction) method is a method reported by Thomas[^4]. The RBV method also assumes that the radioactivity concentration in the tissue is constant.  
First, voxel-wise correction coefficients are obtained by dividing the radioactivity distribution estimated by the GTM method by the image obtained by convolution of the $psf$ . Then, the correction coefficients for each voxel are multiplied by the original image $G$ to obtain the true radioactivity distribution of the tissue.  
[^4]:[Benjamin A Thomas, Kjell Erlandsson, Marc Modat, Lennart Thurfjell, Rik Vandenberghe, Sebastien Ourselin, Brian F Hutton. The importance of appropriate partial volume correction for PET quantification in Alzheimer's disease. Eur J Nucl Med Mol Imaging. 2011 Jun;38(6):1104-19](https://pubmed.ncbi.nlm.nih.gov/21336694/)  

$$\hat{\rho_{RBV}} = \dfrac{\displaystyle \sum_i^n \hat{\rho_{GTM}}(i) M_i}{\displaystyle \sum_i^n \hat{\rho_{GTM}}(i) RSF_i} G$$  
 
### VANUC
Unlike the methods described above, the VANUC method[^1] does not assume uniformity of the radioactivity distribution in the tissue. Instead, it assumes that the radioactivity concentration of the tissue in each voxel follows a certain fixed probability distribution in the neighborhood of the voxel of interest. The true radioactivity distribution of the tissue is obtained by distributing the pixel value $G$ which is a mixture of multiple tissue components due to the partial volume effect, to each tissue component in the most likely proportion based on their probability distributions. If the probability distribution is regarded as following a normal distribution $N(\mu, \sigma^2)$ with a constant expected value and variance, it can be estimated as the sum of the expected value component $\mu$ and the deviation $r$ as follows:  

$$\hat{\rho_{VANUC}}(k) = \hat{\mu}(k) + \hat{r}(k), $$  

$$\hat{r}(k) = \dfrac{\hat{{\sigma_k}^2} B_k}{\displaystyle \sum_t^n {\sigma_t}^2 B_t} \left(\dfrac{G}{RSF_k} - \hat{\mu}(k) \right), $$  

$$B_k = {M_k}^2 \otimes psf^2$$  

Where the expected value $\mu$ is estimated by the generalized least squares method based on multiple voxel data in a neighborhood, and the VANUC program considers the variance-covariance matrix $\mathbf{\Omega}$ as $\sigma^2 \mathbf{I}$ for simplicity:  

$$\begin{align}\mathbf{\hat{\mu}} &= (\mathbf{RSF}^T \mathbf{\Omega}^{-1} \mathbf{RSF} + \lambda^2 \mathbf{I})^{-1} \mathbf{RSF}^T \mathbf{\Omega}^{-1} \mathbf{G}\\
&= (\mathbf{RSF}^T \mathbf{RSF} + \lambda^2 \sigma^2 \mathbf{I})^{-1} \mathbf{RSF}^T \mathbf{G} \end{align}$$  

The variance of the probability distribution $\sigma^2$ is considered to be constant within the tissue region, and is estimated from the variance $s^2$ of the pixel values within each region in the image $G$ , using a principle similar to the GTM method described above: 

$$\mathbf{\hat{V}} = \mathbf{A^{-1}} \mathbf{S}, A_{i,j} = \dfrac{\displaystyle \sum M_i B_j}{\displaystyle \sum M_i}, $$  

$$\mathbf{V} = \begin{pmatrix} {\sigma_1}^2\\
\vdots\\
{\sigma_n}^2 \end{pmatrix}, \mathbf{S} = \begin{pmatrix} {s_1}^2\\
\vdots\\
{s_n}^2 \end{pmatrix}$$  

### Partial volume effect reduced image
In addition to the normal partial volume effect corrected images, the VANUC program outputs partial volume effect "reduced" images. This is an image with a certain resolution (FWHM = 2.5mm in this software) that is processed by applying the principle of partial volume effect correction. This not only reduces the partial volume effect of the original PET or SPECT image with lower resolution, but also produces an image with virtually equal resolution from various images with different conditions. Partial volume effect "reduction" is not simply smoothing of the partial volume effect corrected image.  
Specifically, a partial volume effect "reduced" image $H$ based on the Muller-Gartner or RBV method is obtained using a region spread function of a set resolution $P$ as follows:  

$$H(k) = \hat{\rho}(k) P_k$$  

In the case of the VANUC method, the same operation is performed for the expected value component $\mu$ :  

$$H(k) = \hat{\mu}(k) P_k + \hat{r}(k)$$  

### Parameter Determination in VANUC
The VANUC software is recommended for resolution measurements using the Hoffman brain phantom based on the results of previous studies, but can be applied to other phantoms as well.  
The process flow is as follows.  
1. Coregister the true radioactivity concentration distribution image of the phantom with the actual PET or SPECT image obtained by imaging the phantom.  
1. A model image is created by convolving the true distribution image with a 3-dimensional Gaussian function with a certain full width at half maximum (FWHM).  
1. The difference between the model image and the actual PET or SPECT image is evaluated as the mean squared error (MSE).  
1. Iterate 2. through 3. to search for a FWHM that minimizes the MSE.  

## 3-3. Results of Previous Studies
(in preparation)

# 4. Licenses
When presenting the results obtained by using this software in scientific journals or at conferences, etc., please clearly state the name of the software, etc.  
The developer is not responsible for any damage caused by the use or malfunction of this software.  
Please click [here](/LICENSE) for license terms.

# 5. References
