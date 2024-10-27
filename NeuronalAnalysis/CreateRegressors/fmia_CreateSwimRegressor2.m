function fmia_CreateBehavRegressor
% create behav regressor

global gh

drivename=['E:\'];
load('..\FishAnalysisSummary.mat'); %load fish summary eg excluded session
addpath([drivename,'TransitFMBehavior\fm_behavim\analysis'])

anglebias_overall = [2.7117, 1.9589, 1.5048, 4.6064, 2.2587];


STCounter2=zeros(5,2);
for FishN = 1:length(fish)

    load([drivename,'TransitFMlog\fm',num2str(fish{FishN}.id),'.mat']) %load log
    load([drivename,'TransitFMBehavior\fm',num2str(fish{FishN}.id),'\Registered\',num2str(fish{FishN}.id),'_cleaned.mat'])
    gh.param.rootdir = [drivename,'TransitFMBehavior\'];
    gh.param.framepersession = 8000;
    gh.param.framepersecond = 200;
    gh.param.fishid = num2str(fish{FishN}.id);
    gh.param.fishlog = fishlog;
    % fishlog.trialdetails.trial(23,1)
    % fish{FishN}.id
    if fish{FishN}.id==330
        gh.param.ExcludedSession = [6;9;20];
    elseif fish{FishN}.id==365
        gh.param.ExcludedSession = [5;10;11];
    else
        gh.param.ExcludedSession = [];
    end
    gh.data.anglevect = behavim.cleaneddata.angle;
    gh.data.saccadevect = behavim.cleaneddata.angle_LE;
    gh.data.saccadevect_R = behavim.cleaneddata.angle_RE;
    gh.data.boutmtx = fm_behavim_boutid;
    gh.data.bout_details = fm_behavim_boutcharacterization;
    gh.data.saccademtx = fm_behavim_saccadeid;
    gh.data.simuMtx = fm_behavim_simuBoutSaccadeID;

    discardedframe = fish{FishN}.FinalImFreq*10;
    gh.data.swim = zeros(3,(sum(fishlog.trialdetails.timing)*fish{FishN}.FinalImFreq -discardedframe)*fish{FishN}.TotalImTrial);
    angleifmove = zeros(size(gh.data.anglevect));
    angleifmove2 = zeros(size(gh.data.anglevect));
    angleifmove3 = zeros(size(gh.data.anglevect));

    % load('STCounter2.mat')
    % STprob = STCounter2(FishN,2)./(STCounter2(FishN,2)+STCounter2(FishN,1));

    for ii=1:size(gh.data.bout_details,1)
        if   ~ismember(gh.data.bout_details(ii,1),fish{FishN}.ExcludedSession) & gh.param.fishlog.trialdetails.trial(gh.data.bout_details(ii,1),1)==8 
            % ...
            %     & ~isnan(gh.data.bout_details(ii,5)) ...

            if isempty(find(gh.data.simuMtx(:,1)==gh.data.boutmtx(ii,1) & gh.data.simuMtx(:,9)==gh.data.boutmtx(ii,2)))
                STCounter2(FishN,1)= STCounter2(FishN,1)+1;
            else
                STCounter2(FishN,2)= STCounter2(FishN,2)+1;
            end
        end
        % if ~ismember(gh.data.bout_details(ii,1),fish{FishN}.ExcludedSession)
        %     if ~isnan(gh.data.bout_details(ii,5))
        %         angleifmove(gh.data.bout_details(ii,3):gh.data.bout_details(ii,5)) = 1;
        %     else
        %         angleifmove(gh.data.bout_details(ii,3):(gh.data.bout_details(ii,1))*8000) = 1;
        %     end
        %     if isempty(find(gh.data.simuMtx(:,1)==gh.data.boutmtx(ii,1) & gh.data.simuMtx(:,9)==gh.data.boutmtx(ii,2)))
        %         if ~isnan(gh.data.bout_details(ii,5))
        %             angleifmove2(gh.data.bout_details(ii,3):gh.data.bout_details(ii,5)) = 1;
        %         else
        %             angleifmove2(gh.data.bout_details(ii,3):(gh.data.bout_details(ii,1))*8000) = 1;
        %         end
        %         % STCounter2(FishN,1)= STCounter2(FishN,1)+1;
        %     else
        %         % STCounter2(FishN,2)= STCounter2(FishN,2)+1;
        %     end
        % 
        % 
        %     if rand>=STprob
        %         if ~isnan(gh.data.bout_details(ii,5))
        %             angleifmove3(gh.data.bout_details(ii,3):gh.data.bout_details(ii,5)) = 1;
        %         else
        %             angleifmove3(gh.data.bout_details(ii,3):(gh.data.bout_details(ii,1))*8000) = 1;
        %         end
        %     end
        % end
    end

    % reshaped_angleifmove = reshape(angleifmove,200./fish{FishN}.FinalImFreq,[]);
    % reduced_angleifmove = sum(reshaped_angleifmove,1);
    % gh.data.swim(1,:) = reduced_angleifmove(1:size(gh.data.swim,2));
    % reshaped_angleifmove2 = reshape(angleifmove2,200./fish{FishN}.FinalImFreq,[]);
    % reduced_angleifmove2 = sum(reshaped_angleifmove2,1);
    % gh.data.swim(2,:) = reduced_angleifmove2(1:size(gh.data.swim,2));
    % reshaped_angleifmove3 = reshape(angleifmove3,200./fish{FishN}.FinalImFreq,[]);
    % reduced_angleifmove3 = sum(reshaped_angleifmove3,1);
    % gh.data.swim(3,:) = reduced_angleifmove3(1:size(gh.data.swim,2));


    % 
    % all_indx = 1:size(gh.data.swim,2);
    % delete_indx = mod(all_indx,fish{FishN}.FramePerSession)>0 & mod(all_indx,fish{FishN}.FramePerSession)<=5*fish{FishN}.FinalImFreq;
    % gh.data.swim(:,delete_indx)=[];
    % gh.data.swim(1,:)= normalize(gh.data.swim(1,:),'range');
    % gh.data.swim(2,:)= normalize(gh.data.swim(2,:),'range');
    % gh.data.swim(3,:)= normalize(gh.data.swim(3,:),'range');
    % Swim_Regressor = gh.data.swim;
    save('STCounter3.mat','STCounter2')
    % save([drivename,'FM_IntegratedAnalysis\regressors\fm',num2str(fish{FishN}.id),'_swimB.mat'],'Swim_Regressor')

end
% STCounter2
end