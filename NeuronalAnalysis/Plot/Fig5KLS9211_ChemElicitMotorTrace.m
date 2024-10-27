function Fig5KLS910_ChemElicitMotorTrace

addpath('../util');close all;
global gh

xfish_tracemtx{1}=[]; xfish_tracemtx{2}=[]; xfish_tracemtx{3}=[];

for FishN = 1:length(gh.fish)

    regionlist = find(gh.PulledData{FishN}.ROI_Regressor_XBlur.region.fine==283);
    for ii=1:3
        max_perm{ii} = prctile(gh.PulledData{FishN}.MI_swimB.mi_swim_perm(:,ii),99);
        behavroi{ii} = find(gh.PulledData{FishN}.MI_swimB.mi_swim(:,ii)> max_perm{ii}*gh.param.includthres);
    end

    ROIlist = union(behavroi{1},behavroi{2});
    ROIlist = intersect(regionlist,ROIlist);

    ImFreq = gh.fish{FishN}.FinalImFreq;
    tracemtx{1}=[];    tracemtx{2}=[];    tracemtx{3}=[];

    for ii=1:size(gh.PulledData{FishN}.boutmtx,1)
        % check if it is a S-T event
        if isempty(find(gh.PulledData{FishN}.simuMtx(:,1)==gh.PulledData{FishN}.boutmtx(ii,1) & gh.PulledData{FishN}.simuMtx(:,9)==gh.PulledData{FishN}.boutmtx(ii,2)))
            starttime_ST = gh.PulledData{FishN}.boutmtx(ii,2);
            endtime_ST = gh.PulledData{FishN}.boutmtx(ii,4);
            rounded_ST = round(starttime_ST*2)./2;
            trialn = gh.PulledData{FishN}.boutmtx(ii,1);
            rown = find(gh.PulledData{FishN}.simuMtx(:,1)==gh.PulledData{FishN}.boutmtx(ii,1) & gh.PulledData{FishN}.simuMtx(:,9)==gh.PulledData{FishN}.boutmtx(ii,2));
            STdelay = gh.PulledData{FishN}.simuMtx(rown,9)-gh.PulledData{FishN}.simuMtx(rown,2);

            if ~ismember(gh.PulledData{FishN}.boutmtx(ii,1),gh.fish{FishN}.ExcludedSession) & gh.PulledData{FishN}.boutmtx(ii,1)<=gh.fish{FishN}.TotalImTrial
                if starttime_ST >20 & starttime_ST <=30 & endtime_ST <=30 & gh.PulledData{FishN}.fishlog.trialdetails.trial(trialn,1)~=8   
                    % & STdelay<0.10
                    for jj=1:length(ROIlist)
                        tracexrange = [round(rounded_ST-5)*ImFreq:round(rounded_ST+5)*ImFreq] + (trialn-1)*40*ImFreq;

                        rawtrace = gh.PulledData{FishN}.ROI_Regressor_XBlur.dfof_session(tracexrange,ROIlist(jj));
                        p = polyfit([1:4*ImFreq]',rawtrace(1:4*ImFreq),1);
                        y1 = polyval(p,1:length(rawtrace));
                        d_rawtrace = rawtrace-y1';
                        if (gh.PulledData{FishN}.fishlog.trialdetails.trial(trialn,1)==1) | (gh.PulledData{FishN}.fishlog.trialdetails.trial(trialn,1)==2)
                            tracemtx{2} = [tracemtx{2},d_rawtrace-d_rawtrace(1)];
                        elseif (gh.PulledData{FishN}.fishlog.trialdetails.trial(trialn,1)==4) | (gh.PulledData{FishN}.fishlog.trialdetails.trial(trialn,1)==5)
                            tracemtx{3} = [tracemtx{3},d_rawtrace-d_rawtrace(1)];
                        end
                        clear d_rawtrace rawtrace y1 p
                    end

                end
                if (((gh.PulledData{FishN}.fishlog.trialdetails.trial(trialn,1)==8) & (starttime_ST >10 & starttime_ST <35)) | ...
                        (starttime_ST >10 & starttime_ST <=20))  
                    % & STdelay<0.10
                    for jj=1:length(ROIlist)
                        tracexrange = [round(rounded_ST-5)*ImFreq:round(rounded_ST+5)*ImFreq] + (trialn-1)*40*ImFreq;
                        rawtrace = gh.PulledData{FishN}.ROI_Regressor_XBlur.dfof_session(tracexrange,ROIlist(jj));
                        p = polyfit([1:4*ImFreq]',rawtrace(1:4*ImFreq),1);
                        y1 = polyval(p,1:length(rawtrace));
                        d_rawtrace = rawtrace-y1';

                        tracemtx{1} = [tracemtx{1},d_rawtrace-d_rawtrace(1)];
                        clear d_rawtrace rawtrace y1 p
                    end
                end
            end
        end
        clear STdelay starttime_ST endtime_ST rounded_ST trialn rown
    end

    for ii=1:3
        add_fish_avgtrace=[];
        if FishN~=4
            if size(tracemtx{ii}',1)>1
                add_fish_avgtrace = nanmedian(tracemtx{ii}');
            else
                add_fish_avgtrace = (tracemtx{ii}');
            end
        else
            if size(tracemtx{ii}',1)>1
                averagetrace = nanmedian(tracemtx{ii}');
            else
                averagetrace = (tracemtx{ii}');
            end
            if ~isempty(tracemtx{ii})
                for jj=1:21
                    t=[1,1,2,2,3,3,4,4,5,5,6,7,7,8,8,9,9,10,10,11,11];
                    add_fish_avgtrace(jj) = averagetrace(t(jj));
                end
            end
        end
        xfish_tracemtx{ii} =[xfish_tracemtx{ii};add_fish_avgtrace];
        clear add_fish_avgtrace
        %         [FishN  length(ROIlist) size(xfish_tracemtx{ii})]
    end
end

%Supp 9
% ylim_max = 0.036;
% ylim_min = -0.012;

%Supp 10
ylim_max = 0.030;
ylim_min = -0.018;

%Supp 11
% ylim_max = 0.016;
% ylim_min = -0.008;

%Figure 5
% ylim_max = 0.016;
% ylim_min = -0.008;

%for the singal-trial analysis
% ylim_max = 0.008;
% ylim_min = -0.016;

%for pallium
% ylim_max = 0.012;
% ylim_min = -0.012;

plot([11 11],[ylim_min ylim_max],'color','k','linewidth',1,'linestyle',':');hold on
specifiedcolor{1}=[0.5 0.5 0.5];specifiedcolor{2}=[119 136 172]./255;specifiedcolor{3}=[238 129 114]./255;
stat_threshold = 2.5;
for ii=1:3
    t=xfish_tracemtx{ii};

    if size(t,1)>=3
        if ii==1
            plot(min(max(nanmean(xfish_tracemtx{ii}),ylim_min),ylim_max),'color',specifiedcolor{ii},'linewidth',4,'linestyle',':'); hold on
        else
            plot(min(max(nanmean(xfish_tracemtx{ii}),ylim_min),ylim_max),'color',specifiedcolor{ii},'linewidth',4,'linestyle','-'); hold on
            localmtx=xfish_tracemtx{ii};
            notnancount = sum(isnan(localmtx(:,1))==0);
            pgon = polyshape([1:21,fliplr([1:21])],[nanmean(localmtx)+nanstd(localmtx)./sqrt(notnancount),fliplr(nanmean(localmtx)-nanstd(localmtx)./sqrt(notnancount))]);
            plot(pgon,'FaceColor',specifiedcolor{ii},'FaceAlpha',0.3,'EdgeColor','none');
        end
    end

    % estimate population S.D. with Bessel's correction
    if ii==1
        varpool = reshape(t(:,1:8),[],1);
        estd_varpool = std(varpool).*sqrt(size(t,1))./sqrt(size(t,1)-1);
    end

    if ii==3
        [h3_9,p3_9]=ztest(t(:,9),stat_threshold*estd_varpool,estd_varpool,'Tail','right');
        [h3_10,p3_10]=ztest(t(:,10),stat_threshold*estd_varpool,estd_varpool,'Tail','right');
    elseif ii==2
        [h2_9,p2_9]=ztest(t(:,9),stat_threshold*estd_varpool,estd_varpool,'Tail','right');
        [h2_10,p2_10]=ztest(t(:,10),stat_threshold*estd_varpool,estd_varpool,'Tail','right');
    end
    clear t
end
clear estd_varpool varpool

%display stat
format shortE
[p2_9,p2_10]
[p3_9,p3_10]

set(gcf,'Position',[100 100 300 400]);xticklabels([]); yticklabels([]); box off; axis off
xlim([1 21]); ylim([ylim_min ylim_max])
