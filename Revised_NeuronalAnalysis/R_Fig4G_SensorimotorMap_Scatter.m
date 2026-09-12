function Fig4G_SensorimotorMap_Scatter

addpath('../util');close all

global gh

BrainRegion.qty = [];BrainRegion.fineid = [];
BrainRegion.coarseid = [];BrainRegion.percincludedroi =[];

fineregion_id = [279,283,291,60,58, 13,76,74, 35,66,106, 108,219:225,0];
threshold = 1.3;

senmotdatapair = [];
ROIquadrantCount = zeros(4,2);
capcount = zeros(1,2);
for FishN = 1:length(gh.fish)
    
    for ii=1:5
        perm_stim(ii) = prctile(gh.PulledData{FishN}.MI_stim.mi_stim_perm(:,ii),99);
        sensoryroi{ii} = find(gh.PulledData{FishN}.MI_stim.mi_stim(:,ii)> perm_stim(ii)*gh.param.includthres);
    end
    
    for ii=1:3
        max_perm{ii} = prctile(gh.PulledData{FishN}.MI_swimB.mi_swim_perm(:,ii),99);
        behavroi{ii} = find(gh.PulledData{FishN}.MI_swimB.mi_swim(:,ii)> max_perm{ii}*gh.param.includthres);
    end
    
    norm_mi = zeros(size(gh.PulledData{FishN}.MI_swimB.mi_swim));
    
    for ii=1:3
        norm_mi(:,ii) = gh.PulledData{FishN}.MI_swimB.mi_swim(:,ii)./prctile(gh.PulledData{FishN}.MI_swimB.mi_swim_perm(:,ii),99);
    end
    
    appe_list = intersect(sensoryroi{1},sensoryroi{2});
    aver_list = intersect(sensoryroi{3},sensoryroi{4});
    sens_list = union(appe_list,aver_list);
    sens_list = sens_list(~ismember(sens_list,sensoryroi{5}));
    
    moto_list = union(behavroi{1},behavroi{2});
    
    data1.qty = (norm_mi(:,1)-norm_mi(:,2))-(norm_mi(:,1)-norm_mi(:,3));
    
    data1.qty1 = mean([gh.PulledData{FishN}.MI_stim.mi_stim(:,3)./perm_stim(3),gh.PulledData{FishN}.MI_stim.mi_stim(:,4)./perm_stim(4)]')';
    data1.qty2 = mean([gh.PulledData{FishN}.MI_stim.mi_stim(:,1)./perm_stim(1),gh.PulledData{FishN}.MI_stim.mi_stim(:,2)./perm_stim(2)]')';
    data1.qty3 = data1.qty1-data1.qty2;
    senmotlist = union(sens_list,moto_list);

    pallist = intersect(find(gh.PulledData{FishN}.ROI_Regressor_XBlur.region.fine==283),senmotlist);
    nonpallist = intersect(find(gh.PulledData{FishN}.ROI_Regressor_XBlur.region.fine~=283),senmotlist);
    ROIquadrantCount = ROIquadrantCount + ...
        [length(find(data1.qty(pallist)>0 & data1.qty3(pallist)>0)), length(find(data1.qty(nonpallist)>0 & data1.qty3(nonpallist)>0));...
        length(find(data1.qty(pallist)<0 & data1.qty3(pallist)>0)), length(find(data1.qty(nonpallist)<0 & data1.qty3(nonpallist)>0));...
        length(find(data1.qty(pallist)<0 & data1.qty3(pallist)<0)), length(find(data1.qty(nonpallist)<0 & data1.qty3(nonpallist)<0));...
        length(find(data1.qty(pallist)>0 & data1.qty3(pallist)<0)), length(find(data1.qty(nonpallist)>0 & data1.qty3(nonpallist)<0))];
    
    entirelist = 1:length(data1.qty);
    senmotdatapair = [senmotdatapair; data1.qty(senmotlist), data1.qty3(senmotlist)];
    clear data1 appe_list aver_list sens_list perm_stim sensoryroi max_perm behavroi norm_mi
end

x = senmotdatapair(:,1);
y = senmotdatapair(:,2);

outsideMask = (x < -1) | (x > 1) | (y < -6) | (y > 6);

nOutside = sum(outsideMask)          % number of rows beyond limits
nRows    = size(senmotdatapair, 1)
perc_Outside = nOutside.*100./nRows


[min(senmotdatapair(:,1)) max(senmotdatapair(:,1)) median(senmotdatapair(:,1))]
figure
scatter(max(min(senmotdatapair(:,1),1),-1),max(min(senmotdatapair(:,2),6),-6),2,'MarkerFaceColor',[0.65 0.65 0.65],'MarkerEdgeColor','none');hold on
plot([-1 1],[0 0],'linestyle','--','color',[.25 .25 .25]); hold on
plot([0 0],[-6 6],'linestyle','--','color',[.25 .25 .25]); hold on
ylim([-6 6]); xlim([-1 1]); yticks([-6:6]);xticks([-1:0.2:1]);
xticklabels({});yticklabels({});
set(gcf,'Position',[100 100 400 400])




quadrants = {'Q1','Q2','Q3','Q4'};


% Plot
figure;
b = bar(ROIquadrantCount./sum(ROIquadrantCount), 'grouped');
ROIquadrantCount
% Optional colors
b(1).FaceColor = [167 114 75]./255;
b(2).FaceColor = [1 1 1];
b(1).FaceAlpha = 0.15;
b(2).FaceAlpha = 0.15;
% Axes formatting
box off
set(gca, 'XTickLabel', quadrants, 'FontSize', 12);
% xlabel('Quadrant');
% ylabel('ROI (%)');
xlabel('');
ylabel('');
yticks(0:0.1:0.6)
xticklabels([])
yticklabels([])
ylim([0 0.64])

set(gcf, 'Position', [500 100 800 200])







end