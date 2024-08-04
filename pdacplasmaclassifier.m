%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                       %
%                       PDAC plasma classifier                          %
%                                                                       %
%                   Developed by Meiyappan Solaiyappan                  %
%                Johns Hopkins University, Baltimore, MD 21287          %
%                                                                       %
%           The following MATLAB script (code) is provided under        %
%               JHU Academic Software License Agreement                 %                                     %
%                 Please refer to the license document                  %
%                      included in the repository                       %
%                                                                       %
%  The code & data shared as supplementary material to the manuscript:  %
%  Artificial Neural Network Detection of Pancreatic Cancer             %
%  from 1H MR Spectral Patterns of Plasma Metabolites                   %
%  Solaiyappan, M., et al. Submitted to Nature Communications Medicine  %
%                                                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% the following code shows the classification results that can be obtained 
% from the blinded-test samples and the traning samples, as reported in the 
% manuscript, using the system of ANN based classifiers described therein.

clear
close all
fprintf("\nLoading the training spectral data ...")
load ('trainingspectra.mat');
fprintf("\nLoading the blinded-test spectral data ...")
load ('blindedtestspectra.mat');

% check the mean spectra of each class
% by displaying them together

figure;
ppmaxis = linspace(10,0.5,30142);
set(gca, 'XDir','reverse');
hold on;
plot(ppmaxis,mean(trainingspectra.normal),'g');
plot(ppmaxis,mean(trainingspectra.disease),'b');
plot(ppmaxis,mean(trainingspectra.malignant),'r');

% the classification results are obtained using a set of two classifiers,
% the first one (cbvsm) classifier is optimized for the combined 
% normal & disease versus malignant and the second classifier (cvsb)
% optimized for normal versus disease.
%
fprintf("\nLoading the clasifiers ...")
load ('classifier_cbvsm.mat');
load ('classifier_cvsb.mat');
fprintf ("done.\n")

% construct the vectors of spectral data for classification
trainingspectravec =  [trainingspectra.normal; trainingspectra.disease; trainingspectra.malignant];
trainingspectravec = spectrasort.stdizedata(trainingspectravec);
blindedtestspectravec = spectrasort.stdizedata(blindedtestspectra.classxyz);

% classify the blinded test samples using the network
[blindedtestclassids,blindedtestclassprobs]=spectrasort.resolvejointclasses(blindedtestspectravec,classifier_cbvsm,classifier_cvsb);
% display the confusion matrix of the results 
spectrasort.classifierplotconfusion(blindedtestclassids,blindedtestspectra.classids);

% classfiy the training spectra similar to blinded-test spectra
% and diplay its confusion matrix
% Turn on the flag to test that
test_training_same_as_blinded=false;
if (test_training_same_as_blinded)
  [trainingclassids,trainingclassprobs]=spectrasort.resolvejointclasses(trainingspectravec,classifier_cbvsm,classifier_cvsb);
  spectrasort.classifierplotconfusion(trainingclassids,trainingspectra.classids);
end

% classify the training spectra using classifier optimized for 
% combined normal & disease vs. malignant (cbvsm) and display 
% its confusion matrix. Uncomment the line below to test the code.
% classifier_cbvsm.trainconfusion(trainingspectra);

% classify the training spectra using classifier optimized for
% combined normal & disease vs. malignant (cbvsm) and 
% the normal versus disease (cvsb)
classifier_cbvsm.trainjointconfusion(trainingspectra,classifier_cvsb);

% when testing spectra of samples belonging to different cancer type,
% the bestperfcutoff limit can be reset to 1 to obtain a more relaxed
% condition for the classfiers when generating the final result.
% Here it is applied on the same blinded test data as an example.
% Turn on the flag below to test that.
test_new_spectra = false;
if (test_new_spectra)
    classifier_cbvsm.cmperfcutoff = 1;
    classifier_cbvsm.bmperfcutoff = 1;
    classifier_cbvsm.cbperfcutoff = 1;
    [blindedtestclassids,blindedtestclassprobs]=spectrasort.resolvejointclasses(blindedtestspectravec,classifier_cbvsm,classifier_cvsb);
    spectrasort.classifierplotconfusion(blindedtestclassids,blindedtestspectra.classids);
end
