function fm_behavim_main

global gh

drivename = 'E';
gh.param.rootdir = [drivename,':\TransitFMBehavior\'];

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


% if  gh.param.fishid==365 || gh.param.fishid==366
%     load(['./fishlog/fm',num2str(330),'.mat'])
% else
%     load(['./fishlog/fm',num2str(gh.param.fishid),'.mat'])
% end


load([drivename,':\TransitFMlog\fm',num2str(gh.param.fishid),'.mat'])

gh.param.fishlog = fishlog;
gh.data.anglevect = behavim.cleaneddata.angle;
gh.data.saccadevect = behavim.cleaneddata.angle_LE;
gh.data.saccadevect_R = behavim.cleaneddata.angle_RE;
gh.data.boutmtx = fm_behavim_boutid;
gh.data.bout_details = fm_behavim_boutcharacterization;
gh.data.saccademtx = fm_behavim_saccadeid;
gh.data.simuMtx = fm_behavim_simuBoutSaccadeID;


