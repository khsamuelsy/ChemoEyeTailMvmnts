function Fig5CD_ValencePreference

% Generate map of interested ROI
% add path to get utilizable script in generating the brain map
addpath('../util');close all

global gh

for FishN = 1:length(gh.fish)
    threshold = 1.3;

    for ii=1:5
        perm_stim(ii) = prctile(gh.PulledData{FishN}.MI_stim.mi_stim_perm(:,ii),99);
        sensoryroi{ii} = find(gh.PulledData{FishN}.MI_stim.mi_stim(:,ii)> perm_stim(ii)*threshold);
    end

    data.appe_list = intersect(sensoryroi{1},sensoryroi{2});
    data.aver_list = intersect(sensoryroi{3},sensoryroi{4});
    data.list = union(data.appe_list,data.aver_list);
    data.list = data.list(~ismember(data.list,sensoryroi{5}));

    data.qty1 = mean([gh.PulledData{FishN}.MI_stim.mi_stim(:,3)./perm_stim(3),gh.PulledData{FishN}.MI_stim.mi_stim(:,4)./perm_stim(4)]')';
    data.qty2 = mean([gh.PulledData{FishN}.MI_stim.mi_stim(:,1)./perm_stim(1),gh.PulledData{FishN}.MI_stim.mi_stim(:,2)./perm_stim(2)]')';
    data.qty = data.qty1-(data.qty2);
    length(data.list)

    caloptions.type = 'y';    caloptions.meanflag = true;
    caloptions.FishN = FishN;    caloptions.start = true;
    caloptions.onefishonly = false;

    % fmia_GenerateMapData(data,caloptions)

    for ii=1:length(gh.param.fineregionID)
        if gh.param.fineregionID(ii)~=0

            regionroilist = find(gh.PulledData{FishN}.ROI_Regressor_XBlur.region.fine==gh.param.fineregionID(ii));
            includregionroilist = intersect(regionroilist,data.list);
            regionqty(ii) = nanmean(data.qty(includregionroilist));
        else
            regionroilist = find(gh.PulledData{FishN}.ROI_Regressor_XBlur.region.fine==gh.param.fineregionID(ii) & gh.PulledData{FishN}.ROI_Regressor_XBlur.region.coarse==114 );
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
        plotoptions.crange = [-2 2];

        blueredcolor;
        J = blue_to_white_to_red;
        plotoptions.colormap = J;
        plotoptions.type=caloptions.type ;

        if caloptions.meanflag==false
            mapdata = map;
        else
            mapdata = map./count;
        end
        [max(max(mapdata)) min(min(mapdata))]
        % fmia_BrainMap(mapdata,plotoptions)
        BrainRegion.plotyzero = true;
        BrainRegion.yticks = [-2:0.5:2];
        BrainRegion.percylim = [0 0.4];
        fmia_PlotRegion(BrainRegion)
    end
end

end