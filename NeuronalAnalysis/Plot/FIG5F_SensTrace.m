function FIG5F_SensTrace

% Generate map of interested ROI
% add path to get utilizable script in generating the brain map
anglebias_overall = [2.7117, 1.9589, 1.5048, 4.6064, 2.2587];
addpath('../util')
close all

global gh

xfish_tracemtx =[];xfish_tracemtx2 =[];
cbr_perc = zeros(5,4);
cbr = [275, 1 , 94, 114];

for FishN = 1:length(gh.fish)
      
    for ii=1:5
        perm_stim(ii) = prctile(gh.PulledData{FishN}.MI_stim.mi_stim_perm(:,ii),99);
        sensoryroi{ii} = find(gh.PulledData{FishN}.MI_stim.mi_stim(:,ii)> perm_stim(ii)*gh.param.includthres);
    end

    data.qty1 = mean([gh.PulledData{FishN}.MI_stim.mi_stim(:,3)./perm_stim(3),gh.PulledData{FishN}.MI_stim.mi_stim(:,4)./perm_stim(4)]')';
    data.qty2 = mean([gh.PulledData{FishN}.MI_stim.mi_stim(:,1)./perm_stim(1),gh.PulledData{FishN}.MI_stim.mi_stim(:,2)./perm_stim(2)]')';
    data.qty3 = data.qty1-data.qty2;
    appegroup = intersect(sensoryroi{1},sensoryroi{2});
    avergroup = intersect(sensoryroi{3},sensoryroi{4});
    
    sensgroup = union(appegroup,avergroup);

    sensgroup = sensgroup(~ismember(sensgroup,sensoryroi{5}));
    shortlisted = intersect(sensgroup,find(data.qty3>2));

    for ii=1:4
        totaln = length(shortlisted);
        regionn = length(find(gh.PulledData{FishN}.ROI_Regressor_XBlur.region.coarse(shortlisted)==cbr(ii)));
        cbr_perc(FishN,ii) = regionn./totaln;
    end

    ImFreq = gh.fish{FishN}.FinalImFreq;
    tracemtx = nan(ImFreq*10+1,1000);
    traceidx = 0;
    for ii=1:gh.fish{FishN}.TotalImTrial
        
        if (gh.PulledData{FishN}.fishlog.trialdetails.trial(ii,1)==4 ||  gh.PulledData{FishN}.fishlog.trialdetails.trial(ii,1)==5) & ...
                ~ismember(ii,gh.fish{FishN}.ExcludedSession)
            starttime_ST = 20;
            for jj=1:length(shortlisted)
                traceidx = traceidx+1;
                tracexrange = [round(starttime_ST-5)*ImFreq:round(starttime_ST+5)*ImFreq] + (ii-1)*40*ImFreq;
                rawtrace = gh.PulledData{FishN}.ROI_Regressor_XBlur.dfof_session(tracexrange,shortlisted(jj));
                tracemtx(:,traceidx) = rawtrace-rawtrace(1);
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

    tracemtx = nan(ImFreq*10+1,1000);
    traceidx = 0;
    
    for ii=1:gh.fish{FishN}.TotalImTrial
        
        if (gh.PulledData{FishN}.fishlog.trialdetails.trial(ii,1)==1 ||  gh.PulledData{FishN}.fishlog.trialdetails.trial(ii,1)==2) & ...
                ~ismember(ii,gh.fish{FishN}.ExcludedSession)
            starttime_ST = 20;
            for jj=1:length(shortlisted)
                traceidx = traceidx+1;
                tracexrange = [round(starttime_ST-5)*ImFreq:round(starttime_ST+5)*ImFreq] + (ii-1)*40*ImFreq;
                rawtrace = gh.PulledData{FishN}.ROI_Regressor_XBlur.dfof_session(tracexrange,shortlisted(jj));
                tracemtx(:,traceidx) =  rawtrace-rawtrace(1);
            end
        end
    end
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
    clear data
end
figure
plot([11.5 11.5],[-.004 .016],'color','k','linewidth',1,'linestyle',':');hold on

notnancount = sum(isnan(xfish_tracemtx2(:,1))==0);
pgon = polyshape([1:21,fliplr([1:21])],[nanmean(xfish_tracemtx2)+nanstd(xfish_tracemtx2)./sqrt(notnancount),fliplr(nanmean(xfish_tracemtx2)-nanstd(xfish_tracemtx2)./sqrt(notnancount))]);
plot(pgon,'FaceColor',[119 136 172]./255,'FaceAlpha',0.3,'EdgeColor','none');
plot(nanmean(xfish_tracemtx2),'color',[119 136 172]./255,'linewidth',4); hold on

notnancount = sum(isnan(xfish_tracemtx(:,1))==0);
pgon = polyshape([1:21,fliplr([1:21])],[nanmean(xfish_tracemtx)+nanstd(xfish_tracemtx)./sqrt(notnancount),fliplr(nanmean(xfish_tracemtx)-nanstd(xfish_tracemtx)./sqrt(notnancount))]);
plot(pgon,'FaceColor',[238 129 114]./255,'FaceAlpha',0.3,'EdgeColor','none');
plot(nanmean(xfish_tracemtx),'color',[238 129 114]./255,'linewidth',4); hold on
ylim([-.004 .016]);
yticks([-.004:.001:.016])
xlim([1 21])
xticks([1:21]);
xticklabels([]);yticklabels([]);box off
set(gcf,'Position',[100 700 200 400]); 
axis off
cbr_perc



rownisnan = find(isnan(cbr_perc(:,1)))
cbr_perc(rownisnan,:)=[];


figure(3)
% plot the pie chart
ax = gca();
% normalize the data
pieData = [median(cbr_perc)]./sum(median(cbr_perc));
h = pie(ax, pieData);
% Define 3 colors, one for each of the 3 wedges
newColors = [[[145 79 30]*.5 + [255 255 255]*.5];...
    222 172 128;...
    247 220 185;...
    [181 193 142]]./255;  
ax.Colormap = newColors;

%beautify the pie chart
for ii=2:2:length(ax.Children)
    ax.Children(ii).EdgeAlpha = 0;
    ax.Children(ii).FaceAlpha = 0.7;
end
    delete(ax.Children([1:2:length(ax.Children)]));
    set(gcf,'Position',[400 700 300 300])
end