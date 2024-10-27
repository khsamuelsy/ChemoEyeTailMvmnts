function fmia_LoadFishData

global gh

addpath('..\samuelfcn');
addpath('..\saveastiff_4.5');
drivename='E:\';

%load fish summary eg excluded session
load('..\FishAnalysisSummary.mat')
gh.fish = fish;
gh.param.includthres = 1.3;
for FishN = 1:length(fish)
    load([drivename,'FM_IntegratedAnalysis\calculatedMI\fm',num2str(fish{FishN}.id),'_MI_stim.mat']) % stimulus MI
    load([drivename,'FM_IntegratedAnalysis\calculatedMI\fm',num2str(fish{FishN}.id),'_MI_swimB.mat']) % motor MI
    
    load([drivename,'TransitFMLog\fm',num2str(fish{FishN}.id),'.mat']) % log
    load([drivename,'FM_IntegratedAnalysis\regressors\fm',num2str(fish{FishN}.id),'_ROI_XBlur.mat']) % Calcium data
    load([drivename,'FM_IntegratedAnalysis\func2anat\fm',num2str(fish{FishN}.id),'_func2anat.mat']) % Registration
    load([drivename,'TransitFMCalcium\FMAnatReg\fm',num2str(fish{FishN}.id),'_regdata.mat']) % Registration
    ImLs = double(loadtiff([drivename,'TransitFMCalcium\FMAnatReg\fm',num2str(fish{FishN}.id),'_anatbgm.tif']));
    
    
    gh.PulledData{FishN}.fishlog = fishlog;
    gh.PulledData{FishN}.MI_stim = MI_stim;
    gh.PulledData{FishN}.MI_swimB = MI_swimB;
    
    gh.PulledData{FishN}.ROI_Regressor_XBlur = ROI_Regressor_XBlur;
    gh.PulledData{FishN}.func2anat = func2anat;
    gh.PulledData{FishN}.regdata = regdata;
    gh.PulledData{FishN}.im.ImLs = ImLs;
    
    addpath([drivename,'TransitFMBehavior\fm_behavim\analysis'])
    load([drivename,'TransitFMBehavior\fm',num2str(fish{FishN}.id),'\Registered\',num2str(fish{FishN}.id),'_cleaned.mat']) %behavim
    gh.param.rootdir = [drivename,'TransitFMBehavior\'];
    
    gh.param.framepersession = 8000;
    gh.param.framepersecond = 200;
    gh.param.fishid = fish{FishN}.id;
    gh.param.fishseq = FishN;
    
    gh.param.fineregionID = [279,283,291,60,58,13,76,74,35,66,106, 108,219:225,0];

    gh.param.ExcludedSession =[]; %can set excluded session to zero here, for behavioral analysis but not needed here
    if gh.param.fishid==330
        gh.param.ExcludedSession = [6;9;20];
    elseif gh.param.fishid==365
        gh.param.ExcludedSession = [5;10;11];
    else
        gh.param.ExcludedSession = [];
    end

    gh.PulledData{FishN}.anglevect = behavim.cleaneddata.angle;
    gh.PulledData{FishN}.saccadevect = behavim.cleaneddata.angle_LE;
    gh.PulledData{FishN}.saccadevect_R = behavim.cleaneddata.angle_RE;
    gh.PulledData{FishN}.boutmtx = fm_behavim_boutid;
    gh.PulledData{FishN}.bout_details = fm_behavim_boutcharacterization;
    gh.PulledData{FishN}.saccademtx = fm_behavim_saccadeid;
    gh.PulledData{FishN}.simuMtx = fm_behavim_simuBoutSaccadeID;
    clear behavim fishlog MI_stim MI_swimB ROI_Regressor_XBlur func2anat regdata ImLs
    
    
end
end