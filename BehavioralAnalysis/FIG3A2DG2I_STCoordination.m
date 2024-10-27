clear all; close all; global gh
colorsel;
count_S=0; count_T=0; count_ST=0;
count_S_contra = 0; count_S_ipsi = 0; 
count_T_contra = 0; count_T_ipsi = 0;
count_ST_contra = 0; count_ST_ipsi = 0; 
STdelay=[];
addpath('./stattool');addpath('./disptool');
for fishsub = 1:length(totalfishsub)
    gh.param.fishid = totalfishsub(fishsub);
    fm_behavim_main
    for ii=1:sessionn(fishsub)
        if ~ismember(ii,gh.param.ExcludedSession) 
            if gh.param.fishlog.trialdetails.trial(ii,1)==8
                totalbout_n = find(gh.data.boutmtx(:,1)==ii & gh.data.boutmtx(:,4)<=40);
                totalsaccade_n = find(gh.data.saccademtx(:,1)==ii & gh.data.saccademtx(:,4)<=40);
                totalsimu_n = find(gh.data.simuMtx(:,1)==ii & gh.data.simuMtx(:,4)<=40);
            else
                totalbout_n = find(gh.data.boutmtx(:,1)==ii & gh.data.boutmtx(:,4)<=20);
                totalsaccade_n = find(gh.data.saccademtx(:,1)==ii & gh.data.saccademtx(:,4)<=20);
                totalsimu_n = find(gh.data.simuMtx(:,1)==ii & gh.data.simuMtx(:,4)<=20);
            end
            count_T = count_T+length(totalbout_n);
            count_S = count_S+length(totalsaccade_n);
            count_ST = count_ST+length(totalsimu_n);
            
            % find the proportion of ipsi-contra events
            for saccaden = 2:length(totalsaccade_n)
                if gh.data.saccademtx(totalsaccade_n(saccaden),7) == gh.data.saccademtx(totalsaccade_n(saccaden-1),7)
                    count_S_ipsi = count_S_ipsi+1;
                else
                    count_S_contra = count_S_contra+1;
                end
            end
            
            saved_lastboutdir=[];
            threshold_boutbias = 1.5;
            for boutn = 1:length(totalbout_n)
                boutdir = sign(gh.data.bout_details(totalbout_n(boutn),7)*-1 - anglebias_overall(fishsub));
                boutbias = abs(gh.data.bout_details(totalbout_n(boutn),7)*-1 - anglebias_overall(fishsub));
                if boutbias>=threshold_boutbias
                    if boutn>1 && ~isempty(saved_lastboutdir)
                        if boutdir == saved_lastboutdir
                            count_T_ipsi  = count_T_ipsi +1;
                        else
                            count_T_contra = count_T_contra +1;
                        end
                    end
                    saved_lastboutdir = boutdir;
                else
                    saved_lastboutdir = [];
                end
            end
            
            for simun = 1:length(totalsimu_n)
                boutdir = sign(gh.data.simuMtx(totalsimu_n(simun),14)*-1 - anglebias_overall(fishsub));
                boutbias = abs(gh.data.simuMtx(totalsimu_n(simun),14)*-1 - anglebias_overall(fishsub));
                saccadedir = gh.data.simuMtx(totalsimu_n(simun),7);
                STdelay= [STdelay;gh.data.simuMtx(totalsimu_n(simun),9)-gh.data.simuMtx(totalsimu_n(simun),2), ...
                    gh.data.simuMtx(totalsimu_n(simun),4)-gh.data.simuMtx(totalsimu_n(simun),9), ...
                    gh.data.simuMtx(totalsimu_n(simun),11)-gh.data.simuMtx(totalsimu_n(simun),4), ...
                    (gh.data.simuMtx(totalsimu_n(simun),9)-gh.data.simuMtx(totalsimu_n(simun),2)) - (gh.data.simuMtx(totalsimu_n(simun),11)-gh.data.simuMtx(totalsimu_n(simun),4)), ...
                    gh.data.simuMtx(totalsimu_n(simun),4)-gh.data.simuMtx(totalsimu_n(simun),2), ...
                    gh.data.simuMtx(totalsimu_n(simun),11)-gh.data.simuMtx(totalsimu_n(simun),9)];
                if boutbias>=threshold_boutbias
                    if saccadedir == boutdir
                        count_ST_ipsi = count_ST_ipsi+1;
                    else
                        count_ST_contra = count_ST_contra+1;
                    end
                end
            end
        end
    end
end
figure(1)
ax = gca();
pieData = [count_ST_contra count_ST_ipsi];
h = pie(ax, pieData);
% Define 3 colors, one for each of the 3 wedges
newColors = [...
    .65,       0.65, 0.65;   %hot pink
    0.9,       .9,       .9];  %dark orchid
ax.Colormap = newColors;
set(gcf,'Position',[100 1000 300 300])
figure(2)
ax = gca();
pieData = [count_S_contra count_S_ipsi];
h = pie(ax, pieData);
% Define 3 colors, one for each of the 3 wedges
newColors = [...
    .65,       0.65, 0.65;   %hot pink
    0.9,       .9,       .9];  %dark orchid
ax.Colormap = newColors;
set(gcf,'Position',[400 1000 300 300])
figure(3)
ax = gca();
pieData = [count_T_contra count_T_ipsi];
h = pie(ax, pieData);
% Define 3 colors, one for each of the 3 wedges
newColors = [...
    .65,       0.65, 0.65;   %hot pink
    0.9,       .9,       .9];  %dark orchid
ax.Colormap = newColors;
set(gcf,'Position',[700 1000 300 300])


for ii=1:4
    figure(3+ii)
    STdelay_cropped = max(min(STdelay(:,ii),.5),-.5);
    plot([0 0],[0 0.25],'color',[.65 .65 .65],'linestyle',':','linewidth',5); hold on;
    plot([median(STdelay(:,ii)) median(STdelay(:,ii))],[0 0.25],'color','k','linestyle',':','linewidth',5); hold on
    histogram(STdelay_cropped,[-0.5:0.05:0.5],'normalization','probability','displaystyle','stairs','LineWidth',8,'EdgeColor',STColor); hold on;
    xlim([-0.5 0.5]);ylim([0 0.25]); xticks([-0.5:0.1:0.5]);yticks([0:0.05:0.25]);
    yticklabels([]);xticklabels([]);box off;
    set(gcf,'Position',[1200 1300-ii*300 600 200])
    [min(STdelay(:,ii)) max(STdelay(:,ii)) median(STdelay(:,ii)) var(STdelay(:,ii))]
end

% disp('Signed Rank Stat Test results for a-d')
% signrank(STdelay(:,1));signrank(STdelay(:,2));
% signrank(STdelay(:,3));signrank(STdelay(:,4));

% disp('Bionomial Test results for h-j')
% [h, p] = binotest(num_successes, num_trials, p0);
