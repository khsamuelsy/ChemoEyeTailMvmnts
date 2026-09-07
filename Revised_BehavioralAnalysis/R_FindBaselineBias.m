%this script find biaseline angle bias and undulation int bias
clear all;close all;
global gh
LoadFishNColorSel;addpath('./stattool');addpath('./disptool');

Group1 = [];
Group2 = [];

for fishsub = 1:length(totalfishsub)
    bout_bias=[];      bout_und_int=[];
    gh.param.fishid = totalfishsub(fishsub);
    fm_behavim_main
    for ii=1:sessionn(fishsub)
        if ~ismember(ii,gh.param.ExcludedSession)
            if gh.param.fishlog.trialdetails.trial(ii,1)==8
                totalbout_n = find(gh.data.boutmtx(:,1)==ii & gh.data.boutmtx(:,4)<=40);
            else
                totalbout_n = find(gh.data.boutmtx(:,1)==ii & gh.data.boutmtx(:,4)<=20);
            end
            for boutn = 1:length(totalbout_n)
                bout_bias = [bout_bias;(gh.data.bout_details(totalbout_n(boutn),7)*-1) ];
                X = gh.data.anglevect(gh.data.bout_details(totalbout_n(boutn),3):gh.data.bout_details(totalbout_n(boutn),5));
                X0 = X-max(X);
                X0_inv = -X0;
                pks_X0 = findpeaks(X0);
                pks_X0_inv = findpeaks(X0_inv);
                und_int = gh.data.bout_details(totalbout_n(boutn),6)./((length(pks_X0)+length(pks_X0_inv))./2);
                bout_und_int = [bout_und_int; und_int];
            end
        end
    end
    m_bias = mean(bout_bias);
    m_und  = mean(bout_und_int);
    fprintf('fish %d angle bias is %.6f, und int is %.4f\n', ...
        totalfishsub(fishsub), m_bias, m_und);

    if gh.param.fishid<700
        Group1 =[Group1;m_bias];
    else
        Group2 =[Group2;m_bias];
    end
end