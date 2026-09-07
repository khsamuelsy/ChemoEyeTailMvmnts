clear all; close all; global gh
LoadFishNColorSel;
count_S=0; count_T=0; count_ST=0;
count_sp = 0; count_av=0; count_ap=0;
STdelay=[];
addpath('./stattool');addpath('./disptool');
for fishsub = 1:length(totalfishsub)
    gh.param.fishid = totalfishsub(fishsub);
    fm_behavim_main
    for ii=1:sessionn(fishsub)
        if ~ismember(ii,gh.param.ExcludedSession) 

                totalbout_n = find(gh.data.boutmtx(:,1)==ii);
                totalsaccade_n = find(gh.data.saccademtx(:,1)==ii);
                totalsimu_n = find(gh.data.simuMtx(:,1)==ii);
                
                count_T = count_T+length(totalbout_n);
                count_S = count_S+length(totalsaccade_n);
                count_ST = count_ST+length(totalsimu_n);
                if gh.param.fishlog.trialdetails.trial(ii,1)==8
                    count_sp = count_sp+1;
                elseif gh.param.fishlog.trialdetails.trial(ii,1)==1 | gh.param.fishlog.trialdetails.trial(ii,1)==2 
                    count_ap = count_ap+1;
                elseif gh.param.fishlog.trialdetails.trial(ii,1)==4 | gh.param.fishlog.trialdetails.trial(ii,1)==5
                    count_av = count_av+1;
                end

        
        end
    end
end

