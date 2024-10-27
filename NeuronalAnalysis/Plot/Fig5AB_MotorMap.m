function Fig5AB_MotorMap

% Generate map of interested ROI
% add path to get utilizable script in generating the brain map
addpath('../util');close all;

global gh

for FishN = 1:length(gh.fish)

    % Rename quantities
    for ii=1:3
        motor_perm{ii} = prctile(gh.PulledData{FishN}.MI_swimB.mi_swim_perm(:,ii),99);
        behavroi{ii} = find(gh.PulledData{FishN}.MI_swimB.mi_swim(:,ii)> motor_perm{ii}*gh.param.includthres);
    end

    norm_mi = zeros(size(gh.PulledData{FishN}.MI_swimB.mi_swim));
    for ii=1:3
        norm_mi(:,ii) = gh.PulledData{FishN}.MI_swimB.mi_swim(:,ii)./prctile(gh.PulledData{FishN}.MI_swimB.mi_swim_perm(:,ii),99);
    end

    % Shortlisted ROI and their quantities
    data.list = union(behavroi{1},behavroi{2});

    data.qty = (norm_mi(:,1)-norm_mi(:,2))-(norm_mi(:,1)-norm_mi(:,3));
    % Uncomment if need to check range of qty
    % [max(data.qty(data.list)),min(data.qty(data.list))]
    % Uncomment if need to check ROI no.
    % length(data.list)

    % Specify map plane, mean or max and other options
    caloptions.type = 'y';    caloptions.meanflag = true;
    caloptions.FishN = FishN;    caloptions.start = true;
    caloptions.onefishonly = false;
 % Uncomment if to generate map
    % fmia_GenerateMapData(data,caloptions) 
   
    %Figure 5A

    for ii=1:length(gh.param.fineregionID)
        if gh.param.fineregionID(ii)~=0
            regionroilist = find(gh.PulledData{FishN}.ROI_Regressor_XBlur.region.fine==gh.param.fineregionID(ii));
            includregionroilist = intersect(regionroilist,data.list);
            regionqty(ii) = nanmean(data.qty(includregionroilist));
        else
            regionroilist = find(gh.PulledData{FishN}.ROI_Regressor_XBlur.region.fine==gh.param.fineregionID(ii) & ...
                gh.PulledData{FishN}.ROI_Regressor_XBlur.region.coarse==114 );
            includregionroilist = intersect(regionroilist,data.list);
            regionqty(ii) = nanmean(data.qty(includregionroilist));
        end
        BrainRegion.meanmi(FishN,ii) = regionqty(ii);
    end

    for ii=1:2
        if ii==1
            regionroilist = find(gh.PulledData{FishN}.ROI_Regressor_XBlur.region.fine==283);
        else
            regionroilist = find(gh.PulledData{FishN}.ROI_Regressor_XBlur.region.fine~=283);
        end
        includregionroilist = intersect(regionroilist,data.list);
        regionqty(ii) = nanmean(data.qty(includregionroilist));
        BrainRegion.meanmiC(FishN,ii) = regionqty(ii);
    end


    if FishN == length(gh.fish)

        global map count
        plotoptions.fign = 1;
        plotoptions.crange = [-0.2 0.2];
        greenredcolor;
        J = green_to_white_to_red;

        plotoptions.colormap = J;
        plotoptions.type=caloptions.type ;

        if caloptions.meanflag==false
            mapdata = map;
        else
            mapdata = map./count;
        end

        % max(max(mapdata))
        % min(min(mapdata))
        % fmia_BrainMap(mapdata,plotoptions)
        BrainRegion.plotyzero = true;
        % BrainRegion.ylim = [0 .2];
        BrainRegion.yticks = [-.6:.1:.6];
        BrainRegion.percylim = [0 1];
        fmia_PlotRegion(BrainRegion)
    end
end
% [h, p, ci,stats]=ttest(gh.test.meanmiC(:,1),gh.test.meanmiC(:,2),'Tail','right');

end

