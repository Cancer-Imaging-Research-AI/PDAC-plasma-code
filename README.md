# PDAC-plasma-code
This repository contains **supplementary code** for the article:
### Artificial Neural Network Detection of Pancreatic Cancer from 1H MR Spectral Patterns of Plasma Metabolites

Meiyappan Solaiyappan, Santosh Kumar Bharti, Raj Kumar Sharma, Mohamad Dbouk, Wasay Nizam, Malcolm V. Brock, Michael G. Goggins, Zaver M. Bhujwalla

Submitted to **Nature Communications Medicine**

The MATLAB code for generating the ANN classification results reported in this study using the MR spectral data (included under Data availability), and the system of the trained ANN, are provided here under **JHU Academic Software License Agreement** (included in the repository)

There are three MATLAB code files (pdacplasmaclassifier.m, spectrasort.m, spectranet.m) and two ANN classifier structure files (classifier_cbvsm.mat, classifier_cvsb.mat). Additionally, the copy of the training data (trainingspectra.mat) and the blinded test data (blindedtestspectra.mat) from the Supplementary Data repository are also included.

pdacplasmaclassifier.m is the main file with appropriate comments to run the code.

The code may be potentially tested on any new MRS data with the following provisions. All details related to the plasma sample collection, including sample collection tubes, and the high resolution 1H MR spectra acquisition parameters for the ZGPR sequence, including field strength and spectral processing, will be required to accurately match the details provided in Methods. 

