function fmia_CreateROIRregressor

global gh
load('..\FishAnalysisSummary.mat'); %load fish summary eg excluded session
drivename=['E:\'];
addpath(genpath('..\samuelfcn'));
datadirc = [drivename,'TransitFMCalcium\fmca_data\'];
for FishN = 1:length(fish)

    ROIdirc = [datadirc,'PreprocessedCA\fm',num2str(fish{FishN}.id),'\'];
    imdirc = [datadirc,'RegisteredIm\fm',num2str(fish{FishN}.id),'\'];
    gh.data.roi.region.left = []; gh.data.roi.region.right = [];
    gh.data.roi.region.coarse = []; gh.data.roi.region.fine = [];
    gh.data.roi.region.mask = []; gh.data.roi.region.planen =[];
    gh.data.roi.dfof_session = [];
    
    for ii=1:fish{FishN}.TotalImPlane
        % do this task plane-by-plane
        display([num2str(fish{FishN}.id),' :: ',num2str(ii)])
        ROIfile = [ROIdirc,'fm',num2str(fish{FishN}.id),'_plane',num2str(ii),'.mat'];
        XBlurfile = [drivename,'\FM_IntegratedAnalysis\script\Stage0-XBlur\fm',num2str(fish{FishN}.id),'_XBlur.mat'];
        imfile = [imdirc,'fm',num2str(fish{FishN}.id),'_plane',num2str(ii),'.tif'];
        
        % load the previously cleaned up region data of each ROI
        load(ROIfile);
        load(XBlurfile);
        gh.data.roi.region.left = [gh.data.roi.region.left, CAIm.region.left];
        gh.data.roi.region.right = [gh.data.roi.region.right, CAIm.region.right];
        gh.data.roi.region.coarse = [gh.data.roi.region.coarse, CAIm.region.coarse];
        gh.data.roi.region.fine = [gh.data.roi.region.fine, CAIm.region.fine];
        gh.data.roi.region.mask = [gh.data.roi.region.mask, CAIm.ROI];
        gh.data.roi.region.planen = [gh.data.roi.region.planen, repmat([ii],1,size(CAIm.ROI,2))];
        
        % load the corresponding image
        im = double(readtiffstack(imfile));
        PlaneROIdfof_session = []; 
        
        for jj=1:size(CAIm.ROI,2) 
            ROIdfof_session=[]; 
            ROIRaw = CAIm.ROI(:,jj);
            ROIFull = full(reshape(ROIRaw,[size(im,1),size(im,2)]));
            NormROIFull = ROIFull./sum(ROIFull(:));
            
            for kk=1:size(im,3)
                imNROI = im(:,:,kk).*NormROIFull;
                ROIPxVal(kk) = sum(imNROI(:)); %Total Px value of the ROI jj at time pt kk
            end
            
            for kk=1:size(im,3)
                ROIPxVal_XBlur(kk) = ROIPxVal(XBlur.FramenMtx(ii,kk));
            end
               
            nsession = size(im,3)./fish{FishN}.FramePerSession;
            
            %find the baseline pixel value of each session
            for kk=1:nsession
                if ~ismember(kk,fish{FishN}.ExcludedSession) && ~ismember(kk,XBlur.XTrial{ii})
                    sessionindx = mod(1:fish{FishN}.TotalImN,fish{FishN}.FramePerSession)<=(5*fish{FishN}.FinalImFreq) & mod(1:fish{FishN}.TotalImN,fish{FishN}.FramePerSession)>=2 & ...
                        (floor(([1:fish{FishN}.TotalImN]-1)./fish{FishN}.FramePerSession)+1)==kk;
                    sessionbaseline(kk) = min(ROIPxVal_XBlur(find(sessionindx)));
                else
                    sessionbaseline(kk) = nan;
                end
            end
            
            for kk=1:size(im,3)
                sessionn = floor((kk-1)./fish{FishN}.FramePerSession)+1;
                %dfof_session value of the ROI jj at timept kk
                ROIdfof_session(kk) = (ROIPxVal_XBlur(kk)-sessionbaseline(sessionn))./sessionbaseline(sessionn);

            end
            PlaneROIdfof_session = [PlaneROIdfof_session, ROIdfof_session'];
            clear ROIPxVal ROIPxVal_XBlur sessionbaseline
        end
        gh.data.roi.dfof_session = [gh.data.roi.dfof_session, PlaneROIdfof_session];
        clear PlaneROIdfof_session
    end
    %save ROI signal
    ROI_Regressor_XBlur = gh.data.roi;
    save([drivename,'FM_IntegratedAnalysis\regressors\fm',num2str(fish{FishN}.id),'_ROI_XBlur.mat'],'ROI_Regressor_XBlur')
end
end