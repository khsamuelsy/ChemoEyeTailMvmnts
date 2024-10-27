function FIG5E_SponMotorTrace

% Generate map of interested ROI
% add path to get utilizable script in generating the brain map
anglebias_overall = [2.7117, 1.9589, 1.5048, 4.6064, 2.2587];
addpath('../util');close all

global gh

xfish_tracemtx =[];xfish_tracemtx2 =[];
cbr_perc = zeros(5,4); cbr = [275, 1 , 94, 114];

for FishN = 1:length(gh.fish)

    for ii=1:3
        max_perm{ii} = prctile(gh.PulledData{FishN}.MI_swimB.mi_swim_perm(:,ii),99);
        behavroi{ii} = find(gh.PulledData{FishN}.MI_swimB.mi_swim(:,ii)> max_perm{ii}*gh.param.includthres);
    end

    norm_mi = zeros(size(gh.PulledData{FishN}.MI_swimB.mi_swim));
    for ii=1:3
        norm_mi(:,ii) = gh.PulledData{FishN}.MI_swimB.mi_swim(:,ii)./prctile(gh.PulledData{FishN}.MI_swimB.mi_swim_perm(:,ii),99);
    end

    data.qty = (norm_mi(:,1)-norm_mi(:,2))-(norm_mi(:,1)-norm_mi(:,3));
    shortlisted = intersect(union(behavroi{1},behavroi{2}),find(data.qty>.175));
 % shortlisted = intersect(union(behavroi{1},behavroi{2}),find( data.qty<=-.35));
 % shortlisted = intersect(union(behavroi{1},behavroi{2}),find( data.qty>.175));
    for ii=1:4
        %find total shortlied ROI no
        totaln = length(shortlisted);
        %find roi number that belong to each coarse region
        regionn = length(find(gh.PulledData{FishN}.ROI_Regressor_XBlur.region.coarse(shortlisted)==cbr(ii)));
        cbr_perc(FishN,ii) = regionn./totaln;
    end


    ImFreq = gh.fish{FishN}.FinalImFreq;

    % extract traces for S-T events
    tracemtx = nan(ImFreq*10+1,10000);
    traceidx = 0;
    for ii=1:size(gh.PulledData{FishN}.simuMtx,1)
        starttime_ST = gh.PulledData{FishN}.simuMtx(ii,9);
        % round to the closest 0.5 seconds
        rounded_ST = round(starttime_ST*2)./2;
        trialn = gh.PulledData{FishN}.simuMtx(ii,1);
        %check if lies in the excluded session and fall into displayable
        %spontaneous range
        if ~ismember(gh.PulledData{FishN}.simuMtx(ii,1),gh.fish{FishN}.ExcludedSession)& gh.PulledData{FishN}.simuMtx(ii,1)<=gh.fish{FishN}.TotalImTrial
            if ((gh.PulledData{FishN}.fishlog.trialdetails.trial(gh.PulledData{FishN}.simuMtx(ii,1),1)==8) & (starttime_ST >10 & starttime_ST <35)) | ...
                    (starttime_ST >10 & starttime_ST <=20)
                for jj=1:length(shortlisted)
                    traceidx = traceidx+1;
                    tracexrange = [round(rounded_ST-5)*ImFreq:round(rounded_ST+5)*ImFreq] + (trialn-1)*40*ImFreq;
                    rawtrace = gh.PulledData{FishN}.ROI_Regressor_XBlur.dfof_session(tracexrange,shortlisted(jj));
                    p = polyfit([1:4*ImFreq]',rawtrace(1:4*ImFreq),1);
                    y1 = polyval(p,1:length(rawtrace));
                    d_rawtrace = rawtrace-y1';
                    tracemtx(:,traceidx) = d_rawtrace-d_rawtrace(1);
                end
            end
        end
    end
    add_fish_avgtrace=[];
    if FishN~=4
        add_fish_avgtrace = nanmedian(tracemtx');
    else
        averagetrace = nanmedian(tracemtx');
        for ii=1:21
            t=[1,1,2,2,3,3,4,4,5,5,6,7,7,8,8,9,9,10,10,11,11];
            add_fish_avgtrace(ii) = averagetrace(t(ii));
        end
    end
    xfish_tracemtx =[xfish_tracemtx;add_fish_avgtrace];
    clear tracemtx


    % extract traces for independent tail events
    tracemtx = nan(ImFreq*10+1,10000);
    traceidx = 0;
    for ii=1:size(gh.PulledData{FishN}.boutmtx,1)
        if isempty(find(gh.PulledData{FishN}.simuMtx(:,1)==gh.PulledData{FishN}.boutmtx(ii,1) & gh.PulledData{FishN}.simuMtx(:,9)==gh.PulledData{FishN}.boutmtx(ii,2)))
            starttime_ST = gh.PulledData{FishN}.boutmtx(ii,2);
            % round to the closest 0.5 seconds
            rounded_ST = round(starttime_ST*2)./2;
            trialn = gh.PulledData{FishN}.boutmtx(ii,1);
            if ~ismember(gh.PulledData{FishN}.boutmtx(ii,1),gh.fish{FishN}.ExcludedSession) & gh.PulledData{FishN}.boutmtx(ii,1)<=gh.fish{FishN}.TotalImTrial
                if ((gh.PulledData{FishN}.fishlog.trialdetails.trial(gh.PulledData{FishN}.boutmtx(ii,1),1)==8) & (starttime_ST >10 & starttime_ST <35)) | ...
                        (starttime_ST >10 & starttime_ST <=20)
                    for jj=1:length(shortlisted)
                        traceidx = traceidx+1;
                        tracexrange = [round(rounded_ST-5)*ImFreq:round(rounded_ST+5)*ImFreq] + (trialn-1)*40*ImFreq;
                        rawtrace = gh.PulledData{FishN}.ROI_Regressor_XBlur.dfof_session(tracexrange,shortlisted(jj));
                        p = polyfit([1:5*ImFreq]',rawtrace(1:5*ImFreq),1);
                        y1 = polyval(p,1:length(rawtrace));
                        d_rawtrace = rawtrace-y1';
                        tracemtx(:,traceidx) = d_rawtrace-d_rawtrace(1);
                    end
                end
            end
        end
    end
    %add roi-average trace to the fish population
    add_fish_avgtrace2=[];
    if FishN~=4
        add_fish_avgtrace2 = nanmedian(tracemtx');
    else
        averagetrace2 = nanmedian(tracemtx');
        for ii=1:21
            t=[1,1,2,2,3,3,4,4,5,5,6,7,7,8,8,9,9,10,10,11,11];
            add_fish_avgtrace2(ii) = averagetrace2(t(ii));
        end

    end
    xfish_tracemtx2 =[xfish_tracemtx2;add_fish_avgtrace2];
end



% plot figure
figure(1)
%figure 1 first plot the independent T events
plot([11 11],[-.002 .010],'color','k','linewidth',1,'linestyle',':');hold on
notnancount = sum(isnan(xfish_tracemtx2(:,1))==0)
pgon = polyshape([1:21,fliplr([1:21])],[nanmean(xfish_tracemtx2)+nanstd(xfish_tracemtx2)./sqrt(notnancount),fliplr(nanmean(xfish_tracemtx2)-nanstd(xfish_tracemtx2)./sqrt(notnancount))]);
plot(pgon,'FaceColor',[100 164 177]./255,'FaceAlpha',0.3,'EdgeColor','none');
plot(smooth(nanmean(xfish_tracemtx2),1),'color',[100 164 177]./255,'linewidth',4); hold on

%figure 1 first plot the independent S-T events
notnancount = sum(isnan(xfish_tracemtx(:,1))==0)
pgon = polyshape([1:21,fliplr([1:21])],[nanmean(xfish_tracemtx)+nanstd(xfish_tracemtx)./sqrt(notnancount),fliplr(nanmean(xfish_tracemtx)-nanstd(xfish_tracemtx)./sqrt(notnancount))]);
plot(pgon,'FaceColor',[191 152 149]./255,'FaceAlpha',0.3,'EdgeColor','none');
plot(smooth(nanmean(xfish_tracemtx),1),'color',[191 152 149]./255,'linewidth',4); hold on
ylim([-.002 .010]);yticks([-.02:.001:.010]);xlim([1 21]);xticks([1:21]);
xticklabels([]);yticklabels([]);box off
set(gcf,'Position',[100 700 200 400]);
axis off

rownisnan=find(isnan(cbr_perc(:,1)));
cbr_perc(rownisnan,:)=[];

figure(3)
ax = gca();

pieData = [median(cbr_perc)]./sum(median(cbr_perc));

h = pie(ax, pieData);
% Define 3 colors, one for each of the 3 wedges
newColors = [[[145 79 30]*.5 + [255 255 255]*.5];
    222 172 128;
    247 220 185;
    [181 193 142]]./255;

ax.Colormap = newColors;
for ii=2:2:length(ax.Children)
    ax.Children(ii).EdgeAlpha = 0;
    ax.Children(ii).FaceAlpha = 0.7;
end
delete(ax.Children([1:2:length(ax.Children)]));
set(gcf,'Position',[400 700 300 300])
end