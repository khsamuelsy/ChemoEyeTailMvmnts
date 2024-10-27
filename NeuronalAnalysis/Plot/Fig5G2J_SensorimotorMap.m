function Fig5G2J_SensorimotorMap

addpath('../util');close all

global gh

BrainRegion.qty = [];BrainRegion.fineid = [];
BrainRegion.coarseid = [];BrainRegion.percincludedroi =[];

fineregion_id = [279,283,291,60,58, 13,76,74, 35,66,106, 108,219:225,0];
threshold =1.3;

for FishN = 1:5
    
    for ii=1:5
        perm_stim(ii) = prctile(gh.PulledData{FishN}.MI_stim.mi_stim_perm(:,ii),99);
        sensoryroi{ii} = find(gh.PulledData{FishN}.MI_stim.mi_stim(:,ii)> perm_stim(ii)*threshold);
    end
    
    for ii=1:3
        max_perm{ii} = prctile(gh.PulledData{FishN}.MI_swimB.mi_swim_perm(:,ii),99);
        behavroi{ii} = find(gh.PulledData{FishN}.MI_swimB.mi_swim(:,ii)> max_perm{ii}*threshold);
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
    
    data.list = intersect(find(data1.qty<0 & data1.qty3<0),union(sens_list,moto_list));
    length(data.list)
    data.qty = sqrt(data1.qty.^2+data1.qty3.^2);
    % regionlist = find(gh.PulledData{FishN}.ROI_Regressor_XBlur.region.fine==279);
    
    % data.list = intersect(data.list,regionlist);
    caloptions.type = 'z';
    caloptions.meanflag = true;
    caloptions.FishN = FishN;
    caloptions.start = true;
    caloptions.onefishonly = false;
     
     
    % fmia_GenerateMapData(data,caloptions)
    
    for ii=1:length(fineregion_id)
        if fineregion_id(ii)~=0
            regionroilist = find(gh.PulledData{FishN}.ROI_Regressor_XBlur.region.fine==fineregion_id(ii));
            includregionroilist = intersect(regionroilist,data.list);
            regionqty(ii) = nanmean(data.qty(includregionroilist));
        else
            regionroilist = find(gh.PulledData{FishN}.ROI_Regressor_XBlur.region.fine==fineregion_id(ii) & gh.PulledData{FishN}.ROI_Regressor_XBlur.region.coarse==114 );
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
        plotoptions.crange = [0 6];
        J = magma(70);  J = J(7:70,:);
        plotoptions.colormap = flip(J,1);
        plotoptions.type=caloptions.type ;
        
        if caloptions.meanflag==false
            mapdata = map;
        else
            mapdata = map./count;
        end
        [max(max(mapdata)) min(min(mapdata))]
        % fmia_BrainMap(mapdata,plotoptions)
        % fmia_BrainMap_ForebrainZoomed(mapdata,plotoptions)
        BrainRegion.plotyzero = false;
        % BrainRegion.yticks = [0:.3:3.6];
        BrainRegion.yticks = [0:.2:2.4];
        BrainRegion.percylim = [0 0.75];
        fmia_PlotRegion(BrainRegion)
    end
%     ylim([0.5 1050.5])
end

end