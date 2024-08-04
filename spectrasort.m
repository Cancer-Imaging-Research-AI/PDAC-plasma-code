%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                       %
%                           spectrasort.m class                         %
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

classdef spectrasort
%
% class members represent the cross validation (cv) based distribution of
% the classifier such as the associated cv training samples and cv test samples
% and the corresponding classifiers: normal vs. malignant (denoted as cm), 
% disease vs. malignant (denoted as bm), and normal vs disease (denoted as cb).
% Additionally, includes information on pivot samples that are persistent 
% members of the training set.
%
    properties
        version;

        ctrdim;
        bendim;
        maldim;
        cvrunidx;        
        ctrpivotsz;
        benpivotsz;
        malpivotsz;
        ctrpivot=[]; 
        benpivot=[];
        malpivot=[];

        ctrstocksz; 
        benstocksz;
        malstocksz;
        ctrstock=[];
        benstock=[];
        malstock=[];

        ctrsplitsz; 
        bensplitsz; 
        malsplitsz;
                                              
        trainvec = {};
        testvec =  {};

        ncvruns;
        ncvrunactual;

        netarray_cm = struct([]);
        netarray_bm = struct([]);
        netarray_cb = struct([]);

        cmnetbestperf=[];
        cmnetperfsortndx = [];

        bmnetbestperf=[];
        bmnetperfsortndx = [];

        cbnetbestperf=[];
        cbnetperfsortndx = [];

        cutoffcv;
        cmperfcutoff; 
        bmperfcutoff;
        cbperfcutoff;
       
        cmperfminoccurcvs=[]; 
        bmperfminoccurcvs=[];
        cbperfminoccurcvs=[]; 

    end
    
    methods (Static)

        function [classids, classprobs] = resolvejointclasses(testspectravector, cmclassifier, cbclassifier)
        %
        % the member function takes the input test spectra array and resolves 
        % the class of the individual spectrum by processing them through two
        % different classifiers each one fine tuned for two different optimal
        % classification. The first classifier normal versus malignant (cmclassifier)
        % is optimized primarily for normal versus malignant (while it can still
        % process normal versus disease) and the second classifier normal versus
        % disease (cbclassifier) that is optimized for the normal versus disease
        % classification. The combined clasifiers results are provided in the output
        % results (classids), along with the class probabilities (classprobs)
        %
            nt2samples = size(testspectravector,1);

            t2cmflags = zeros(cmclassifier.ncvrunactual,nt2samples);
            t2bmflags = zeros(cmclassifier.ncvrunactual,nt2samples);
            t2cbflags = zeros(cmclassifier.ncvrunactual,nt2samples);

            gentblmsgname = "Generating classification table ...";
            classifiername = "normal vs. malignant classifier:";

            cmncvbest = 0;
            bmncvbest = 0;
            cbncvbest = 0;

            tstart = tic;
            wb = waitbar(0, classifiername, 'Name', gentblmsgname);
            for j=1:cmclassifier.ncvrunactual
                if (mod(j,10)==0)
                  waitbar(j/cmclassifier.ncvrunactual,wb);
                end

                t2Y = cmclassifier.testspectra(j,testspectravector);
                if (cmclassifier.netarray_cm(j).bestperf <= cmclassifier.cmperfcutoff)
                    t2cmflags(j,:) = t2Y(1,:);
                    cmncvbest = cmncvbest+1;
                end
                if (cmclassifier.netarray_bm(j).bestperf <= cmclassifier.bmperfcutoff)
                    t2bmflags(j,:) = t2Y(2,:);
                    bmncvbest = bmncvbest+1;
                end
            end

            delete(wb);
            classifiername = "normal vs. disease classifier:";
            wb = waitbar(0, classifiername, 'Name', gentblmsgname);
            
            for j=1:cbclassifier.ncvrunactual
                if (mod(j,10)==0)                
                   waitbar(j/cbclassifier.ncvrunactual,wb);
                end

                t2Y = cbclassifier.testspectra(j,testspectravector);

                if (cbclassifier.netarray_cb(j).bestperf <= cbclassifier.cbperfcutoff)
                    t2cbflags(j,:) = t2Y(3,:);
                    cbncvbest = cbncvbest + 1;
                end

            end

            t2Ycm1sum = sum(t2cmflags==1);
            t2Ycm2sum = sum(t2cmflags==2);
            t2Ybm1sum = sum(t2bmflags==1);
            t2Ybm2sum = sum(t2bmflags==2);
            t2Ycb1sum = sum(t2cbflags==1);
            t2Ycb2sum = sum(t2cbflags==2);

            t2Ycmfmax = ones(1,nt2samples);
            t2Ycmfmax(t2Ycm2sum >= t2Ycm1sum) = 2;

            NM_Dfrac = t2Ycm2sum/cmncvbest;

            t2Ybmfmax = ones(1,nt2samples);
            t2Ybmfmax(t2Ybm2sum >= t2Ybm1sum) = 2;
            DM_Mfrac = t2Ybm2sum/bmncvbest;

            t2Ycbfmax = ones(1,nt2samples);
            t2Ycbfmax(t2Ycb2sum >= t2Ycb1sum) = 2;
            ND_Dfrac = t2Ycb2sum/cbncvbest;

            classids = zeros(1,nt2samples);
            classprobs = zeros(3,nt2samples);

            classids(1,(t2Ycmfmax ==2) & (t2Ybmfmax==2)) = 3;
            classids(1,(classids(1,:)~=3) & (t2Ycbfmax==1)) = 1;
            classids(1,(classids(1,:)~=3) & (t2Ycbfmax==2)) = 2;

            minNMDMfrac = min(NM_Dfrac,DM_Mfrac);
            classprobs(3,:) = round(minNMDMfrac*100,2);

            classprobs(1,:) = (1-classprobs(3,:)/100).*round((1-ND_Dfrac)*100,2);
            classprobs(2,:) = (1-classprobs(3,:)/100).*round(ND_Dfrac*100,2);

            delete(wb);
            fprintf ("\nCompleted. Time taken = %4.2f sec.\n",toc(tstart));
            
        end

        function classifierplotconfusion (testvector,targetvector)
        %
        % this member function takes in testvector (i.e., classifcation results)
        % and the targetvector (the true +ve classification of the samples) and
        % displays the confusion matrix both the original 3-way classification
        % (of normal, disease, malignant) and the 2-way classification of
        % normal and disease combined, versus malignant.
        %
            if (size(testvector,1) == 1)
                testvector3= zeros(3,size(testvector,1));
                testvector3(1,(testvector==1)) = 1;
                testvector3(2,(testvector==2)) = 1;
                testvector3(3,(testvector==3)) = 1;
            else
                testvector3 = testvector;
            end

            if (size(targetvector,1) == 1)
                targetvector3= zeros(3,size(targetvector,1));
                targetvector3(1,(targetvector==1)) = 1;
                targetvector3(2,(targetvector==2)) = 1;
                targetvector3(3,(targetvector==3)) = 1;
            else
                targetvector3 = targetvector;
            end

            figure;
            cfmatrix3x3 = plotconfusion(targetvector3,testvector3);
            cfmatrix3x3.CurrentAxes.XTickLabel{1,1}='normal';
            cfmatrix3x3.CurrentAxes.XTickLabel{2,1}='disease';
            cfmatrix3x3.CurrentAxes.XTickLabel{3,1}='malignant';
            cfmatrix3x3.CurrentAxes.XTickLabelRotation=0;
            cfmatrix3x3.CurrentAxes.YTickLabel{1,1}='normal';
            cfmatrix3x3.CurrentAxes.YTickLabel{2,1}='disease';
            cfmatrix3x3.CurrentAxes.YTickLabel{3,1}='malignant';
            cfmatrix3x3.CurrentAxes.YTickLabelRotation=90;
            cfmatrix3x3.CurrentAxes.XLabel.String = 'Input Class';

            chartgrey = [217,217,217]./255;
            changegrey = [255,255,255]./255;

            set(findall(gcf,'-property','FontSize'),'FontSize',16);
            set(findobj(gcf,'facecolor',chartgrey),'facecolor',changegrey)
            set(findobj(gcf,'color',[0.94,0.94,0.94]),'color',changegrey);
            set(findobj(gcf,'facecolor',[240,240,240]/255),'facecolor',changegrey);
            set(findall(gcf,'LineWidth',2),'LineWidth',3); % these are the special boundary thick lines
            set(findall(gcf,'LineWidth',0.5),'LineWidth',2); % these are the default box lines
            set(findall(gcf,'LineWidth',1),'LineWidth',2); % these are the default box lines

            testvector2way = [testvector3(1,:)+testvector3(2,:);testvector3(3,:)];
            targetvector2way = [targetvector3(1,:)+targetvector3(2,:);targetvector3(3,:)];

            figure;
            cfmatrix2x2 = plotconfusion(targetvector2way,testvector2way);
            cfmatrix2x2.CurrentAxes.XTickLabel{1,1}='normal & disease';
            cfmatrix2x2.CurrentAxes.XTickLabel{2,1}='malignant';
            cfmatrix2x2.CurrentAxes.XTickLabelRotation=0;
            cfmatrix2x2.CurrentAxes.YTickLabel{1,1}='normal & disease';
            cfmatrix2x2.CurrentAxes.YTickLabel{2,1}='malignant';
            cfmatrix2x2.CurrentAxes.YTickLabelRotation=90;
            cfmatrix2x2.CurrentAxes.XLabel.String = 'Input Class';
            
            set(findall(gcf,'-property','FontSize'),'FontSize',16);
            set(findobj(gcf,'facecolor',chartgrey),'facecolor',changegrey)
            set(findobj(gcf,'color',[0.94,0.94,0.94]),'color',changegrey);
            set(findobj(gcf,'facecolor',[240,240,240]/255),'facecolor',changegrey);
            set(findall(gcf,'LineWidth',2),'LineWidth',3); % these are the special boundary thick lines
            set(findall(gcf,'LineWidth',0.5),'LineWidth',2); % these are the default box lines
            set(findall(gcf,'LineWidth',1),'LineWidth',2); % these are the default box lines

        end

        function [spectraout] = stdizedata (spectrain)
        %
        % the member function takes input spectra and standardize the dynamic
        % range of the data using the mean and standard deviation and additionally
        % suppresses the water and EDTA related signals
        %
            if (isstruct(spectrain))
                spectrawowater = spectrain.classx;
                spectrain = spectrain.classx;
            else
                spectrawowater = spectrain;
            end

            % index2ppm = linspace(10,0.5,30142);
            % ppm2index = @(x) int64(1+(10-x)/(10-0.5)*30141);

            watersig = 15866:17452; %ppm: 4.9996 - 4.4997
            edta1 = 20149:20419;    %ppm: 3.6496 - 3.5645
            edta2 = 21482:21640;    %ppm: 3.2295 - 3.1797
            edta3 = 23123:23182;    %ppm: 2.7123 - 2.6937
            edta4 = 23479:23641;    %ppm: 2.6001 - 2.5490

            spectrawowater(:,edta4) = [];
            spectrawowater(:,edta3) = [];
            spectrawowater(:,edta2) = [];
            spectrawowater(:,edta1) = [];
            spectrawowater(:,watersig) = [];

            meanspwowater = mean(spectrawowater,2);
            stdspwowater = std(spectrawowater,1,2);

            spectraout = (spectrain - meanspwowater)./stdspwowater;
            spectraout(:,edta1) = 0;
            spectraout(:,edta2) = 0;
            spectraout(:,edta3) = 0;
            spectraout(:,edta4) = 0;
            spectraout(:,watersig) = 0;

        end
        

    end

    methods       

        function t2Y = testspectra(obj, cvridx, testvector)
        %
        % this is an intermediary member function called by the higher level
        % resolvejointclasses member function that basically performs
        % calling lower level function that captures the feature vector
        % (based on the current classifier's spectral region) and the
        % the corresponding prediction function.
        %
            panzgprnet_cm = obj.netarray_cm(cvridx);
            panzgprnet_bm = obj.netarray_bm(cvridx);
            panzgprnet_cb = obj.netarray_cb(cvridx);

            if (~isempty(panzgprnet_cm.predictnn))
                testfv = spectranet.capturefeaturevec(testvector, panzgprnet_cm);
                t2zgpr_cm = panzgprnet_cm.predictnn(testfv);
                [Yt2zgpr_fcm,Yt2zgpr_cm] = max(t2zgpr_cm);
            else
                Yt2zgpr_cm = zeros(1,size(testvector,1));
            end

            if (~isempty(panzgprnet_bm.predictnn))
                testfv = spectranet.capturefeaturevec(testvector, panzgprnet_bm);
                t2zgpr_bm = panzgprnet_bm.predictnn(testfv);
                [Yt2zgpr_fbm,Yt2zgpr_bm] = max(t2zgpr_bm);
            else
                Yt2zgpr_bm = zeros(1,size(testvector,1));
            end

            testfv = spectranet.capturefeaturevec(testvector, panzgprnet_cb);
            t2zgpr_cb = panzgprnet_cb.predictnn(testfv);
            [Yt2zgpr_fcb,Yt2zgpr_cb] = max(t2zgpr_cb);
            t2Y = [Yt2zgpr_cm; Yt2zgpr_bm; Yt2zgpr_cb];            
        end


        function [test1flags, t1accu] = testoftrainspectra(obj, cvridx, testvec,testdim,splittestvec)
        %
        % this member function addresses the specific problem of testing classifiers
        % using the training samples. The function identifies in each crossvalidation 
        % the training samples and excludes them from consideration in the cross validation
        % results table. Thus the results table consists of only the classifcation results
        % from the test samples in each cross validation run. This helps to avoid any 
        % potential bias arising from samples that were included in the training in a
        % given cross validation instance.
        %
        % this intermediary level function is called by the trainconfusion function
        % included below.
        %
            panzgprnet_cm = obj.netarray_cm(cvridx);
            panzgprnet_bm = obj.netarray_bm(cvridx);
            panzgprnet_cb = obj.netarray_cb(cvridx);

            if (isempty(panzgprnet_cm.predictnn))
                Yt1zgpr_cm = [ones(1,testdim(1)),2*ones(1,testdim(2)),2*ones(1,testdim(3))];
            else
                testfv_cm = spectranet.capturefeaturevec(testvec, panzgprnet_cm);
                t1zgpr_cm = panzgprnet_cm.predictnn(testfv_cm);
                [Yt1zgpr_fcm,Yt1zgpr_cm] = max(t1zgpr_cm);
            end

            if (isempty(panzgprnet_bm.predictnn))
                Yt1zgpr_bm = [ones(1,testdim(1)),ones(1,testdim(2)),2*ones(1,testdim(3))];
            else
                testfv_bm = spectranet.capturefeaturevec(testvec, panzgprnet_bm);
                t1zgpr_bm = panzgprnet_bm.predictnn(testfv_bm);
                [Yt1zgpr_fbm,Yt1zgpr_bm] = max(t1zgpr_bm);
            end

            testfv_cb = spectranet.capturefeaturevec(testvec, panzgprnet_cb);
            t1zgpr_cb = panzgprnet_cb.predictnn(testfv_cb);
            [Yt1zgpr_fcb,Yt1zgpr_cb] = max(t1zgpr_cb);
            
           if (panzgprnet_cm.bestperf > obj.cmperfcutoff) 
                Yt1zgpr_cm = Yt1zgpr_cm.*0;
            end
            if (panzgprnet_bm.bestperf > obj.bmperfcutoff)
                Yt1zgpr_bm = Yt1zgpr_bm.*0;
            end
            if (panzgprnet_cb.bestperf > obj.cbperfcutoff)
                Yt1zgpr_cb = Yt1zgpr_cb.*0;
            end

            test1mask1 = [zeros(1,obj.ctrsplitsz),ones(1,obj.ctrpivotsz)]; 
            test1mask2 = [zeros(1,obj.bensplitsz),ones(1,obj.benpivotsz)];
            test1mask3 = [zeros(1,obj.malsplitsz),ones(1,obj.malpivotsz)];

            test1mask1(splittestvec{1}') = 1;
            test1mask2(splittestvec{2}') = 1;
            test1mask3(splittestvec{3}') = 1;

            test1masks = [test1mask1,test1mask2,test1mask3]; 

            nend = testdim(1);
            bstart = 1+nend;
            bend = sum(testdim(1:2));
            mstart = bend +1;
            mend  = sum(testdim(1:3));

            test1vecsz= [obj.ctrpivotsz+size(splittestvec{1},1), ...
                obj.benpivotsz+size(splittestvec{2},1), ...
                obj.malpivotsz+size(splittestvec{3},1)]; 

            test1_cm = test1masks.* Yt1zgpr_cm;
            t1cm_n = round(sum(test1_cm(1:nend)==1)/test1vecsz(1)*100,2);
            t1cm_m = round(sum(test1_cm(mstart:mend)==2)/test1vecsz(3)*100,2);
            t1accu(1,1) = t1cm_n;
            t1accu(1,2) = t1cm_m;

            test1_bm = test1masks.* Yt1zgpr_bm;
            t1bm_b = round(sum(test1_bm(bstart:bend)==1)/test1vecsz(2)*100,2);
            t1bm_m = round(sum(test1_bm(mstart:mend)==2)/test1vecsz(3)*100,2);
            t1accu(1,3) = t1bm_b;
            t1accu(1,4) = t1bm_m;

            test1_cb = test1masks.* Yt1zgpr_cb;
            t1cb_c = round(sum(test1_cb(1:nend)==1)/test1vecsz(1)*100,2);
            t1cb_b = round(sum(test1_cb(bstart:bend)==2)/test1vecsz(2)*100,2);
            t1accu(1,5) = t1cb_c;
            t1accu(1,6) = t1cb_b;

            test1flags = [test1_cm; test1_bm; test1_cb];
        end

        function [t1cmflags,t1bmflags,t1cbflags] = trainconfusion(obj,trainspectra, jointconfusion, classifiertype)
        %
        % the member function, using the input training spectra data, 
        % and based on the classifier of it is amember, constructs the 
        % confusion matrix of the classification results, calling
        % testoftrainspectra above, that takes care of making sure that
        % in the cross validation based runs only the results of test samples
        % are counted in, in evaluating the final classification result.
        %
            if (nargin < 3)
                jointconfusion = false;
                classifiertype = 1;
            end
                            
            trainspectra.malignant=[trainspectra.malignant(obj.malstock,:); ...
                trainspectra.malignant(obj.malpivot,:)];

            trainspectra.normal=[trainspectra.normal(obj.ctrstock,:); ...
                trainspectra.normal(obj.ctrpivot,:)];

            trainspectra.disease=[trainspectra.disease(obj.benstock,:); ...
                trainspectra.disease(obj.benpivot,:)];

            testtrainzgprvec =  [trainspectra.normal; trainspectra.disease; trainspectra.malignant];
            testtrainzgprvec = spectrasort.stdizedata(testtrainzgprvec);
        
            nsamples = obj.ctrdim + obj.bendim + obj.maldim;

            t1cmflags = zeros(obj.ncvrunactual,nsamples);
            t1bmflags = zeros(obj.ncvrunactual,nsamples);
            t1cbflags = zeros(obj.ncvrunactual,nsamples);

            gentblmsgname = "Generating classification table ...";
            switch (classifiertype)
                case 1
                classifiername = "normal vs. malignant classifier on training data";            
                case 2
                classifiername = "normal vs. disease classifier on training data";            
            end
            tstart = tic;
            wb = waitbar(0, classifiername, 'Name', gentblmsgname);
            
            for j=1:obj.ncvrunactual

                if (mod(j,10)==0)
                    waitbar(j/obj.ncvrunactual, wb);
                end
                
                [t1flags, t1accu] = obj.testoftrainspectra(j, testtrainzgprvec, ... 
                    trainspectra.classdim,obj.testvec(j,:));
                t1cmflags(j,:) = t1flags(1,:);
                t1bmflags(j,:) = t1flags(2,:);
                t1cbflags(j,:) = t1flags(3,:);

            end

            t1cmfmax = ones(1,nsamples);
            t1bmfmax = ones(1,nsamples);
            t1cbfmax = ones(1,nsamples);
            
            t1cmfmax(sum(t1cmflags==2) >= sum(t1cmflags==1)) = 2;          
            t1bmfmax(sum(t1bmflags==2) >= sum(t1bmflags==1)) = 2;
            t1cbfmax(sum(t1cbflags==2) >= sum(t1cbflags==1)) = 2;       

            delete(wb);
            if (jointconfusion)
                return;
            end

            t1_cfmx_3d = zeros(3,nsamples);
            t1_cfmx_3d(3,(t1cmfmax ==2) & (t1bmfmax==2)) = 1;
            t1_cfmx_3d(1,(~t1_cfmx_3d(3,:)) & (t1cbfmax==1)) = 1;
            t1_cfmx_3d(2,(~t1_cfmx_3d(3,:)) & (t1cbfmax==2)) = 1;

            t1_itarget = [repmat([1;0;0],[1,obj.ctrdim]),repmat([0;1;0],[1,obj.bendim]),...
                repmat([0;0;1],[1,obj.maldim])];

            spectrasort.classifierplotconfusion(t1_cfmx_3d,t1_itarget);

            findmisclassified = false;
            %if needed, this flag can be turned on to determine the misclassified samples
            
            if (findmisclassified)
                malorigseq =[obj.malstock,obj.malpivot];
                benorigseq =[obj.benstock,obj.benpivot];
                ctrorigseq =[obj.ctrstock,obj.ctrpivot];

                misclass_m = sort(union(intersect(find(t1_cfmx_3d(2,:)==1),find(t1_itarget(3,:)==1)),intersect(find(t1_cfmx_3d(1,:)==1),find(t1_itarget(3,:)==1))));
                if (~isempty(misclass_m))
                    fprintf ("\nMisclassified malignant samples:\n");
                    fprintf ("%d,",malorigseq(misclass_m-obj.ctrdim-obj.bendim) + obj.ctrdim+obj.bendim);
                end
                misclass_b = find(t1bmfmax(obj.ctrdim+1:obj.ctrdim+obj.bendim)==2)+obj.ctrdim;
                if (~isempty(misclass_b))
                    fprintf ("\nMisclassified benign samples (as malignant):\n");
                    fprintf ("%d,",benorigseq(misclass_b-obj.ctrdim) + obj.ctrdim);
                end
                misclass_c = find(t1cmfmax(1:obj.ctrdim)==2);
                if (~isempty(misclass_c))
                    fprintf ("\nMisclassified normal samples (as disease/malignant):\n");
                    fprintf ("%d,",ctrorigseq(misclass_c));
                end
                misclass_c_b = find(t1cbfmax(1:obj.ctrdim)==2);
                if (~isempty(misclass_c_b))
                    fprintf ("\nMisclassified normal samples as disease\n");
                    fprintf ("%d,",ctrorigseq(misclass_c_b));
                end
                misclass_b_c = find(t1cbfmax(obj.ctrdim+1:obj.ctrdim+obj.bendim)==1)+obj.ctrdim;
                if (~isempty(misclass_b_c))
                    fprintf ("\nMisclassified disease samples as normal\n");
                    fprintf ("%d,",benorigseq(misclass_b_c-obj.ctrdim) + obj.ctrdim);
                end
            end
            
            fprintf ("\nCompleted. Time taken = %4.2f sec.\n",toc(tstart));
           
        end

        function trainjointconfusion(obj,trainspectra, cbclassifier)
        %
        % this member function is similar to the trainconfusion above, with
        % the difference it includes the additional classifier that is optimized
        % for normal versus disease. Thus when invoked as a member function of the
        % optimized normal versus malignant classifier, it can combine optimization of
        % both classifiers (similar to resolvejointclasses) as applied to the
        % training spectral data 
        %
        
            nsamples = obj.ctrdim + obj.bendim + obj.maldim;
            
            jointconfusion = true;
            tstart = tic;

            [t1cmflags,t1bmflags,~] = obj.trainconfusion(trainspectra,jointconfusion,1);
            [~,~,t1cbflags] = cbclassifier.trainconfusion(trainspectra,jointconfusion,2);

            t1cmfmax = ones(1,nsamples);
            t1bmfmax = ones(1,nsamples);
            t1cbfmax = ones(1,nsamples);
            
            t1cmfmax(sum(t1cmflags==2) >= sum(t1cmflags==1)) = 2;          
            t1bmfmax(sum(t1bmflags==2) >= sum(t1bmflags==1)) = 2;
            t1cbfmax(sum(t1cbflags==2) >= sum(t1cbflags==1)) = 2;       


            t1_cfmx_3d = zeros(3,nsamples);
            t1_cfmx_3d(3,(t1cmfmax ==2) & (t1bmfmax==2)) = 1;
            t1_cfmx_3d(1,(~t1_cfmx_3d(3,:)) & (t1cbfmax==1)) = 1;
            t1_cfmx_3d(2,(~t1_cfmx_3d(3,:)) & (t1cbfmax==2)) = 1;

            t1_itarget = [repmat([1;0;0],[1,obj.ctrdim]),repmat([0;1;0],[1,obj.bendim]),...
                repmat([0;0;1],[1,obj.maldim])];
            t1y_target = repmat(t1_itarget,[1,obj.ncvrunactual]);

            spectrasort.classifierplotconfusion(t1_cfmx_3d,t1_itarget);

            malorigseq =[obj.malstock,obj.malpivot];
            benorigseq =[obj.benstock,obj.benpivot];
            ctrorigseq =[obj.ctrstock,obj.ctrpivot];

            findmisclassified = false;
            %if needed, this flag can be turned on to determine the misclassified samples

            if (findmisclassified)
                misclass_m = sort(union(intersect(find(t1_cfmx_3d(2,:)==1),find(t1_itarget(3,:)==1)),intersect(find(t1_cfmx_3d(1,:)==1),find(t1_itarget(3,:)==1))));
                if (~isempty(misclass_m))
                    fprintf ("\nMisclassified malignant samples:\n");
                    fprintf ("%d,",malorigseq(misclass_m-obj.ctrdim-obj.bendim) + obj.ctrdim+obj.bendim);
                end
                misclass_b = find(t1bmfmax(obj.ctrdim+1:obj.ctrdim+obj.bendim)==2)+obj.ctrdim;
                if (~isempty(misclass_b))
                    fprintf ("\nMisclassified benign samples (as malignant):\n");
                    fprintf ("%d,",benorigseq(misclass_b-obj.ctrdim) + obj.ctrdim);
                end
                misclass_c = find(t1cmfmax(1:obj.ctrdim)==2);
                if (~isempty(misclass_c))
                    fprintf ("\nMisclassified normal samples (as disease/malignant):\n");
                    fprintf ("%d,",ctrorigseq(misclass_c));
                end
                misclass_c_b = find(t1cbfmax(1:obj.ctrdim)==2);
                if (~isempty(misclass_c_b))
                    fprintf ("\nMisclassified normal samples as disease\n");
                    fprintf ("%d,",ctrorigseq(misclass_c_b));
                end
                misclass_b_c = find(t1cbfmax(obj.ctrdim+1:obj.ctrdim+obj.bendim)==1)+obj.ctrdim;
                if (~isempty(misclass_b_c))
                    fprintf ("\nMisclassified disease samples as normal\n");
                    fprintf ("%d,",benorigseq(misclass_b_c-obj.ctrdim) + obj.ctrdim);
                end
            end
            fprintf ("\nCompleted. Time taken = %4.2f sec.\n",toc(tstart));

        end

    end

end



