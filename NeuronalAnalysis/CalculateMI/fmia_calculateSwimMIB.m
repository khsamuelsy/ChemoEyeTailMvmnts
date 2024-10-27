function fmia_calculateMI


load('..\FishAnalysisSummary.mat'); %load fish summary eg excluded session
%run each fish 1-by-1
for FishN = 1:length(fish)
    global gh
    
    drivename=['E:\'];
    addpath('..\util'); %add path
    load([drivename,'FM_IntegratedAnalysis\regressors\fm',num2str(fish{FishN}.id),'_swimB.mat']); %load behav regressor
    load([drivename,'FM_IntegratedAnalysis\regressors\fm',num2str(fish{FishN}.id),'_ROI_XBlur.mat']); %load ROI df/f
    load([drivename,'TransitFMlog\fm',num2str(fish{FishN}.id),'.mat']) %load log  
    XBlurfile = [drivename,'\FM_IntegratedAnalysis\script\Stage0-XBlur\fm',num2str(fish{FishN}.id),'_XBlur.mat'];
    load(XBlurfile);
    % remove the first 5 second
    indx = find(mod(1:size(ROI_Regressor_XBlur.dfof_session,1),fish{FishN}.FramePerSession)<=5*fish{FishN}.FinalImFreq & ...
        mod(1:size(ROI_Regressor_XBlur.dfof_session,1),fish{FishN}.FramePerSession)>0);
    ROI_Regressor_XBlur.dfof_session(indx,:) = [];
    ROI_Regressor_XBlur.dfof_session = ROI_Regressor_XBlur.dfof_session';
    
    
    StimSession = find(fishlog.trialdetails.trial(:,1)~=8);

    gh.data.Swim_Regressor=Swim_Regressor;    gh.data.ROI_Regressor=ROI_Regressor_XBlur;

    for ii=1:size(gh.data.ROI_Regressor.dfof_session,1)
        tic
        display([num2str(fish{FishN}.id),' :: ',num2str(ii)])
        
        % Find the MI for the 5 cases
        for jj=1:3
            swimvec = gh.data.Swim_Regressor(jj,:)';
            ROI_PxData = gh.data.ROI_Regressor.dfof_session(ii,:)';
            ExcludedSession = union(fish{FishN}.ExcludedSession,XBlur.XTrial{gh.data.ROI_Regressor.region.planen(ii)}');
            
            %remove 1) excluded session and 2) stim period of stimsession 
            deleteindx = ismember(floor(([1:size(ROI_Regressor_XBlur.dfof_session,2)]-1)./(35*fish{FishN}.FinalImFreq))+1,ExcludedSession) | ...
                (ismember(floor(([1:size(ROI_Regressor_XBlur.dfof_session,2)]-1)./(35*fish{FishN}.FinalImFreq))+1,StimSession) & ...
                (mod(1:size(ROI_Regressor_XBlur.dfof_session,2),(35*fish{FishN}.FinalImFreq))>(15*fish{FishN}.FinalImFreq) | ...
                mod(1:size(ROI_Regressor_XBlur.dfof_session,2),(35*fish{FishN}.FinalImFreq))==0));
            
            swimvec(deleteindx) = [];
            ROI_PxData(deleteindx) =[];
            ROI_PxData_perm = ROI_PxData(randperm(length(ROI_PxData)));
            mi_swim(ii,jj) = MutualInfo2(swimvec,ROI_PxData);
            swim_perm(ii,jj) = MutualInfo2(swimvec,ROI_PxData_perm);
        end
        toc
    end
    gh.data.cal.mi_swim = mi_swim;
    gh.data.cal.mi_swim_perm = swim_perm;
    MI_swimB = gh.data.cal;
    save([drivename,'FM_IntegratedAnalysis\calculatedMI\fm',num2str(fish{FishN}.id),'_MI_swimB.mat'],'MI_swimB')
    clear global gh
    clearvars -except fish
end
end