function fmia_calculateStimMI

load('..\FishAnalysisSummary.mat'); %load fish summary eg excluded session


%run each fish 1-by-1
for FishN = 1:length(fish)
    
    global gh
    drivename=['E:\'];
    addpath('..\util'); %add path
    load([drivename,'FM_IntegratedAnalysis\regressors\fm',num2str(fish{FishN}.id),'_stim.mat']); %load stim regressor
    load([drivename,'FM_IntegratedAnalysis\regressors\fm',num2str(fish{FishN}.id),'_ROI_XBlur.mat']); %load ROI df/f
    load([drivename,'TransitFMlog\fm',num2str(fish{FishN}.id),'.mat']) %load log
    XBlurfile = [drivename,'\FM_IntegratedAnalysis\script\Stage0-XBlur\fm',num2str(fish{FishN}.id),'_XBlur.mat'];
    load(XBlurfile);
    % remove the first 5 second
    indx = find(mod(1:size(ROI_Regressor_XBlur.dfof_session,1),fish{FishN}.FramePerSession)<=5*fish{FishN}.FinalImFreq & ...
        mod(1:size(ROI_Regressor_XBlur.dfof_session,1),fish{FishN}.FramePerSession)>0);
    ROI_Regressor_XBlur.dfof_session(indx,:) = [];
    ROI_Regressor_XBlur.dfof_session = ROI_Regressor_XBlur.dfof_session';
    
    gh.data.Stim_Regressor=Stim_Regressor;    gh.data.ROI_Regressor=ROI_Regressor_XBlur;
    for ii=1:size(gh.data.ROI_Regressor.dfof_session,1)
        tic
        display([num2str(fish{FishN}.id),' :: ',num2str(ii)])
        
        % Find the MI for the 5 cases
        for jj=1:5
            stimvec = gh.data.Stim_Regressor.case(jj,:)';
            ROI_PxData = gh.data.ROI_Regressor.dfof_session(ii,:)';
            ExcludedSession = union(fish{FishN}.ExcludedSession,XBlur.XTrial{gh.data.ROI_Regressor.region.planen(ii)}');
            deleteindx = ismember(floor(([1:size(ROI_Regressor_XBlur.dfof_session,2)]-1)./(35*fish{FishN}.FinalImFreq))+1,ExcludedSession);
            % remove excluding trials from the vect.
            stimvec(deleteindx) = [];
            ROI_PxData(deleteindx) =[];
            ROI_PxData_perm = ROI_PxData(randperm(length(ROI_PxData)));
            mi_stim(ii,jj) = MutualInfo2(stimvec,ROI_PxData);
            stim_perm(ii,jj) = MutualInfo2(stimvec,ROI_PxData_perm);
        end
        toc
    end
    gh.data.cal.mi_stim = mi_stim;
    gh.data.cal.mi_stim_perm = stim_perm;
    MI_stim_PR = gh.data.cal;
    save([drivename,'FM_IntegratedAnalysis\calculatedMI\fm',num2str(fish{FishN}.id),'_MI_stim_PR.mat'],'MI_stim_PR')
    clear global gh
    clearvars -except fish
end
end