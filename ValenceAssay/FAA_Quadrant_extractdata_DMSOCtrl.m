function datasummary=FAA_Quadrant_extractdata_DMSOCtrl
clc; clear all;
load('DMSOCtrl.mat');
nodor = numel(fieldnames(odor));
fnames = fieldnames(odor);

for ii=1:nodor
    odorname = fnames{ii};
    nassay = length(eval(['odor.',odorname]));
    ts =[];    nfish=0;

    for jj=1:nassay
        load(['E:\TransitFMValence\FAA_QuadrantV2V1\valence_results_20220304\WaterDMSOMerged\', ...
            eval(['odor.',odorname,'{',num2str(jj),'}']),'_filtered.mat']);
        Im_Px = zeros(400,400);
        data.FhCoor = FhData_filtered.FhCoor;        data.FhAng = FhData_filtered.FhAng;
        data.nfish = length(FhData_filtered.FhCoor);        nfish = data.nfish+nfish;
        for kk=1:data.nfish
            yy = FhData_filtered.FhCoor{1,kk}(:,1);            xx = FhData_filtered.FhCoor{1,kk}(:,2);
            for mm=1:36000
                if yy(mm)<=400 && yy(mm)>0 && ...
                        xx(mm)<=400 && xx(mm)>0
                    Im_Px(yy(mm),xx(mm)) = Im_Px(yy(mm),xx(mm))+1;

                end
            end
        end
        Im_Px(:,:)= Im_Px(:,:)/sum(sum(Im_Px(:,:)));
        ts(jj,1) = sum(sum(Im_Px(1:200,201:400)));
        ts(jj,2) = sum(sum(Im_Px(1:200,1:200)));
        ts(jj,3) = sum(sum(Im_Px(201:400,1:200)));
        ts(jj,4) = sum(sum(Im_Px(201:400,201:400)));
    end
    datasummary{ii}.odorname = odorname;    datasummary{ii}.nassay = nassay;
    datasummary{ii}.ts = ts;    datasummary{ii}.nfish = nfish;

end
end