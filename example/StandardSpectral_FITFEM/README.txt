Example of reconstruction using a tissue components reconstruction by means of fit with fem spectral forward model.

This example computes a fem forward model for a ground truth given by the mask contained in ORIGINALkwaved. Pseudoheterogeneities can also be added. Reconstruction is operated by a fit fem on the components where the two regions are defined in DTkwaved which has been obtained by segmenting and extruding an Ultrasound simulation (by means of the sofware k-wave) of the ground truth in ORIGINALkwaved.

To operate a segmentation and substitute the files ORIGINALkwave and DTkwaved, cd to src/UCLutils/, run Gen3Dshape and select a DICOM image from those in the folder ./US_simulated and operate the segmentation. Once the operation is concluded in the folder a file DICOM*3D.mat is generated while the ground truth ORIGINAL*.mat is already present. You can then rename the files and copy them in the example folders or set the paths of the two priors in Init_DOT and RecSetting. 
