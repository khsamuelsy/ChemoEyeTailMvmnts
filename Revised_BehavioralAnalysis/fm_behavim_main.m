function fm_behavim_main

global gh

drivename = 'E:\[Important] FishMix';
gh.param.rootdir = [drivename,'\TransitFMBehavior\'];

gh.param.framepersession = 8000;
gh.param.framepersecond = 200;

if gh.param.fishid==330
    gh.param.ExcludedSession = [6;9;20];
elseif gh.param.fishid==365
    gh.param.ExcludedSession = [5;10;11];
else
    gh.param.ExcludedSession = [];
end

load([gh.param.rootdir,'fm',num2str(gh.param.fishid),'\Registered\',num2str(gh.param.fishid),'_cleaned.mat'])
load([drivename,'\TransitFMlog\fm',num2str(gh.param.fishid),'.mat'])

gh.param.fishlog = fishlog;
gh.data.anglevect = behavim.cleaneddata.angle;
gh.data.saccadevect = behavim.cleaneddata.angle_LE;
gh.data.saccadevect_R = behavim.cleaneddata.angle_RE;

% Behavioral timestamp adjustment for revision-added fish (fish IDs >= 700).
% These three fish were acquired with a 100-ms StreamPix polling interval.
% To align behavioral traces to the midpoint of this interval, we apply a
% 50-ms adjustment (10 frames at 200 Hz), assuming a uniform timing offset
% within the polling interval. Tail and eye traces are shifted later by
% 10 frames within each trial; the initial sample is repeated to preserve
% the original trial length.

if gh.param.fishid>=700
    % parameters
    Ntrial = 8000;      % frames per trial
    shiftN = 10;        % 50 ms at 200 Hz

    % pull concatenated vectors
    ang_all  = behavim.cleaneddata.angle(:);
    sacL_all = behavim.cleaneddata.angle_LE(:);
    sacR_all = behavim.cleaneddata.angle_RE(:);

    % number of trials
    nTrials = numel(ang_all) / Ntrial;

    % reshape to [nTrials x Ntrial]
    angMat  = reshape(ang_all,  Ntrial, nTrials).';  % nTrials x 8000
    sacLMat = reshape(sacL_all, Ntrial, nTrials).';
    sacRMat = reshape(sacR_all, Ntrial, nTrials).';

    % pad front with first value per trial and shift within each trial
    angMat_shift  = [repmat(angMat(:,1),  1, shiftN), angMat(:,1:end-shiftN)];
    sacLMat_shift = [repmat(sacLMat(:,1), 1, shiftN), sacLMat(:,1:end-shiftN)];
    sacRMat_shift = [repmat(sacRMat(:,1), 1, shiftN), sacRMat(:,1:end-shiftN)];

    % reshape back to 1D vectors and assign into gh.data
    gh.data.anglevect      = reshape(angMat_shift.',  [], 1);
    gh.data.saccadevect    = reshape(sacLMat_shift.', [], 1);
    gh.data.saccadevect_R  = reshape(sacRMat_shift.', [], 1);
end


gh.data.boutmtx = fm_behavim_boutid;
gh.data.bout_details = fm_behavim_boutcharacterization;
gh.data.saccademtx = fm_behavim_saccadeid;
gh.data.simuMtx = fm_behavim_simuBoutSaccadeID;
