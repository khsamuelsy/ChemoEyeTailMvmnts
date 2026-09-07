clc;clear all;close all
global gh
addpath('./stattool');addpath('./disptool');
gh.param.fishid = 330;
fm_behavim_main
anglebias_overall = 1.9589; % data for fish 330
angle_bias = gh.data.bout_details(:,7).*-1-anglebias_overall;

selprctile =[1,25,50,75,99];
calibratedanglevect = gh.data.anglevect.*-1-anglebias_overall;
figure
for ii=length(selprctile)
    list = find(abs(angle_bias-prctile(angle_bias,selprctile(ii)))==min(abs(angle_bias-prctile(angle_bias,selprctile(ii)))));
    selidx(ii) = list(1);
 
    traceX = calibratedanglevect(gh.data.bout_details(selidx(ii),3):gh.data.bout_details(selidx(ii),5))+65*ii-calibratedanglevect(gh.data.bout_details(selidx(ii),3));
    X0 = traceX-max(traceX);
    X0_inv = -X0;
    [pks_X0 loc_pos] = findpeaks(X0);
    [pks_X0_inv loc_neg] = findpeaks(X0_inv);
    plot(traceX,'color','k','LineWidth',2.5); hold on
    plot([1 length(traceX)],[mean(traceX(loc_pos)) mean(traceX(loc_pos))],'LineStyle',':','color','k','linewidth',1)
    plot([1 length(traceX)],[mean(traceX(loc_neg)) mean(traceX(loc_neg))],'LineStyle',':','color','k','linewidth',1)
    plot([1 length(traceX)],[mean([traceX(1) traceX(end)]) mean([traceX(1) traceX(end)])],'LineStyle',':','color','k','linewidth',1)
    c=mean([traceX(1) traceX(end)])-(+65*ii-calibratedanglevect(gh.data.bout_details(selidx(ii),3)));
    
    a=mean(traceX(loc_pos))-(+65*ii-calibratedanglevect(gh.data.bout_details(selidx(ii),3)));
    
    b=mean(traceX(loc_neg))-(+65*ii-calibratedanglevect(gh.data.bout_details(selidx(ii),3)));
    dur = gh.data.bout_details(selidx(ii),6);
%     length(traceX)
%     angle_bias(list(1))
ncycle = (length(pks_X0)+length(pks_X0_inv))./2;

end

set(gcf,'Position',[300 300 300 100])
yticks([-65.*(length(selprctile)):65:-65]);
box off;xticklabels([]);yticklabels([])

h = gca;
h.XAxis.Visible = 'off';
h.YAxis.Visible = 'off';