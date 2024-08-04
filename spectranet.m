%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                       %
%                           spectranet.m class                          %
%                                                                       %
%                   Developed by Meiyappan Solaiyappan                  %
%                Johns Hopkins University, Baltimore, MD 21287          %
%                                                                       %
%           The following MATLAB script (code) is provided under        %
%               JHU Academic Software License Agreement                 %                                     %
%                 Please refer to the license document                  %
%                      included in the repository                       %
%                                                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


classdef spectranet
    %
    % the class members represents the spectral region of interest
    % in the given instance of the class and the ANN prediction function
    % based on that spectral regions and also the corresponding best 
    % performance obtained in the training
    %
   properties
        spectraROIs;
        nupROIs;
        ndownROIs;
        baseclassmean = [];
        bestperf;
        predictnn;
    end
    
    methods (Static)   

        function [featurevec] = capturefeaturevec(spectrastackstdzd, spnet)
        %
        % the member function extracts the feature vector from the input spectra
        % based on the spectral region of interests (ROIs) provided in the 
        % input structure
        %
            nsamples = size(spectrastackstdzd,1);

            nregions = size(spnet.spectraROIs,2);
            featurevec = zeros(nregions, nsamples);

            for iregion = 1:nregions
                iregionstart = spnet.spectraROIs (1, iregion);
                iregionend = iregionstart +  spnet.spectraROIs(2, iregion) - 1;


                regionslice=spectrastackstdzd(:,iregionstart:iregionend);
                baseslice = spnet.baseclassmean(iregionstart:iregionend);
                for isample = 1:nsamples
                    featurevec(iregion, isample) = sum (regionslice(isample,:)-baseslice);
                end
            end
        end        

    end
    
end
   

    

