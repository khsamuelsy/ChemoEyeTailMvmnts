function fmia_CreateStimRegressor
% create stim regressor

global gh

drivename=['E:\'];
load('..\FishAnalysisSummary.mat'); %load fish summary eg excluded session

for FishN = 1:length(fish)
    load([drivename,'TransitFMlog\fm',num2str(fish{FishN}.id),'.mat']) %load log
    StimCase = unique(fishlog.trialdetails.trial,'rows');
    ncase = size(StimCase*2,1); % 5 for multiple odor expts
    nframepersession = fish{FishN}.FinalImFreq*35;
    % remove the first 15 sec  / 10 sec light evoked stim +5 sec baseline finding
    discardedframe = fish{FishN}.FinalImFreq*15;
    gh.data.stim.case = zeros(ncase,(sum(fishlog.trialdetails.timing)*fish{FishN}.FinalImFreq -discardedframe)*fish{FishN}.TotalImTrial);
    
    for ii=1:ncase
        %create a vector of total frame number
        totalframe = 1:(sum(fishlog.trialdetails.timing)-discardedframe./fish{FishN}.FinalImFreq)*fish{FishN}.FinalImFreq *fish{FishN}.TotalImTrial;
        %find trial that belongs to case ii
        trialn = find(fishlog.trialdetails.trial(:,1)== StimCase(ii,1) & fishlog.trialdetails.trial(:,2)== StimCase(ii,2));
        %set 1 for the time on interval
        on_indx = mod(totalframe,nframepersession)>(fishlog.trialdetails.timing(1).*fish{FishN}.FinalImFreq -discardedframe) ...
            & mod(totalframe,nframepersession)<=(sum(fishlog.trialdetails.timing(1:2)).*fish{FishN}.FinalImFreq -discardedframe) ...
            & ismember(floor((totalframe-1)./nframepersession)+1,trialn);
        %set 1 for the time off interval
        %     off_indx = mod(totalframe,nframepersession)>(sum(fishlog.trialdetails.timing(1:2)).*fish{FishN}.FinalImFreq -discardedframe) ...
        %         & mod(totalframe,nframepersession)<(sum(fishlog.trialdetails.timing(1:3)).*fish{FishN}.FinalImFreq -discardedframe) ...
        %         & ismember(floor((totalframe-1)./nframepersession)+1,trialn);
        gh.data.stim.case(ii,on_indx) = 1;
        %     gh.data.stim.case(ii+5,off_indx) = 1;
    end
    % save stim regressor
    Stim_Regressor = gh.data.stim;
    save([drivename,'FM_IntegratedAnalysis\regressors\fm',num2str(fish{FishN}.id),'_stim.mat'],'Stim_Regressor')
end
end