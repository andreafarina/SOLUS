RECONSTRUCTION FROM EXPERIMENTAL DATA

0th order Tikhonov regularisation:  
 - Install SOLUS 
 - Set your directory to a folder in free/ . 
   Each folder contains initialisation files for the reconstruction of the experimental data coming from different combinations of inclusion/bulk.
 - The initialisation file RECsettings.m can be	modified to allow for different parameters of reconstruction 
 - Run run_DOT.m

Reconstruction with US prior
 - Install SOLUS    
 - Set your directory to a folder in /dt_prior_0.01/ .
   Each folder contains initialisation files for the reconstruction of the experimental data coming from different combinations of inclusion/bulk and a file containing the extrapolated 3D shape.
   The extrapolated shape from US is already set in the in the initialisation files
 - The initialisation file RECsettings.m can be modified to allow for different parameters of reconstruction
 - The path to the extrapolated 3D shape is set in REC.solver.prior.path =''  
 - Run run_DOT.m

Reconstruction with CVylindrical prior
 - Install SOLUS
 - Set your directory to a folder in /cylprior_0.01/ .
   Each folder contains initialisation files for the reconstruction of the experimental data coming from different combinations of inclusion/bulk and a file containing the extrapolated 3D shape.
   The use of a cylindrical prior is already set in the initialisation files
 - The initialisation file RECsettings.m can be modified to allow for different parameters of reconstruction
 - The path to the extrapolated 3D shape is set in REC.solver.prior.path =''
 - Run run_DOT.m



Generation of a 3D prior from US images
 - set your current directory to ../../src/UCLutils
 - run Gen3DShape.m by setting the correct input image. It will run the snake and the DT based extrapolation to generate an output *.mat file containing information about the 3D extrapolated shape.
   In order to generate a 3D shape from Experimental data gathered on phantom set loadname to be DICOMimages/US_image_dicom_cylinder.dcm  
 - set the correct path in the initialisation file to use the generated prior in reconstruction with US priors
 


 

