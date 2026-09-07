clear all; close all;
global gh
LoadFishNColorSel;rng(1)

addpath('./stattool');
addpath('./disptool');

count_ST=0;count_ST_contra = 0; count_ST_ipsi = 0;
STdelay = [];

stimcase = [1,2];
% stimcase = [8];

if mean(stimcase)==1.5
    plotcolor = appeColor;
elseif mean(stimcase)==4.5
    plotcolor = averColor;
else
    plotcolor = blankColor;
end

nFish = length(totalfishsub);

% bookkeeping
STdelay_fish = cell(nFish,2);
ST_ipsi_fish   = zeros(nFish,1);
ST_contra_fish = zeros(nFish,1);

DelayAll  = [];
FishIDAll = [];
TypeAll   = [];

Y_S_all = [];   FishID_S_all = [];
Y_T_all = [];   FishID_T_all = [];
Y_ST_all = [];  FishID_ST_all = [];

for fishsub = 1:nFish
    gh.param.fishid = totalfishsub(fishsub);
    fm_behavim_main
    ST_ipsi_this = 0; ST_contra_this = 0;

    for ii = 1:sessionn(fishsub)
        if ~ismember(ii,gh.param.ExcludedSession)

            if ismember(gh.param.fishlog.trialdetails.trial(ii,1),stimcase)
                totalsimu_n    = find(gh.data.simuMtx(:,1)==ii    & gh.data.simuMtx(:,2)>20    & gh.data.simuMtx(:,4)<=30);
            else
                totalsimu_n    = [];
            end

            count_ST = count_ST + length(totalsimu_n);
            threshold_boutbias = 1.5;
            %% ST events + delays
            for simun = 1:length(totalsimu_n)
                boutdir    = sign(gh.data.simuMtx(totalsimu_n(simun),14)*-1 - anglebias_overall(fishsub));
                boutbias   = abs(gh.data.simuMtx(totalsimu_n(simun),14)*-1 - anglebias_overall(fishsub));
                saccadedir = gh.data.simuMtx(totalsimu_n(simun),7);

                delayvec = [ gh.data.simuMtx(totalsimu_n(simun),9)-gh.data.simuMtx(totalsimu_n(simun),2), ...
                    gh.data.simuMtx(totalsimu_n(simun),4)-gh.data.simuMtx(totalsimu_n(simun),9)];

                STdelay = [STdelay; delayvec];

                for kk = 1:2
                    STdelay_fish{fishsub,kk} = [STdelay_fish{fishsub,kk}; delayvec(kk)];
                    DelayAll  = [DelayAll; delayvec(kk)];
                    FishIDAll = [FishIDAll; fishsub];
                    TypeAll   = [TypeAll; kk];
                end

                if boutbias >= threshold_boutbias
                    if saccadedir == boutdir
                        count_ST_ipsi = count_ST_ipsi + 1;
                        ST_ipsi_this  = ST_ipsi_this + 1;
                        Y_ST_all      = [Y_ST_all; 1];
                    else
                        count_ST_contra = count_ST_contra + 1;
                        ST_contra_this  = ST_contra_this + 1;
                        Y_ST_all        = [Y_ST_all; 0];
                    end
                    FishID_ST_all = [FishID_ST_all; fishsub];
                end
            end
        end
    end
    ST_ipsi_fish(fishsub)   = ST_ipsi_this;
    ST_contra_fish(fishsub) = ST_contra_this;
end

%% fish-level fractions
ST_frac_fish = ST_ipsi_fish ./ (ST_ipsi_fish + ST_contra_fish);

%% fish-level fraction plots 
minSep = 0.06;  step   = 0.04;  markerSize = 28;

%% ORIGINAL pie figures preserved
figure(1)
ax = gca();
pieData = [count_ST_contra count_ST_ipsi];
h = pie(ax, pieData);
newColors = [.65 0.65 0.65; 0.9 0.9 0.9];
ax.Colormap = newColors;
set(gcf,'Position',[100 1050 300 300])

% Fishlevel ST fraction
figure(2);
set(gcf,'Position',[100 650 220 300])
ax = axes('Position',[0.25 0.15 0.55 0.72]); hold(ax,'on')
valid = ~isnan(ST_frac_fish) & isfinite(ST_frac_fish);
y = ST_frac_fish(valid);
x = ones(size(y));

[ys, order] = sort(y);
xoffset = zeros(size(ys));
groupStart = 1;
while groupStart <= numel(ys)
    groupEnd = groupStart;
    while groupEnd < numel(ys) && abs(ys(groupEnd+1)-ys(groupEnd)) < minSep
        groupEnd = groupEnd + 1;
    end
    nGroup = groupEnd - groupStart + 1;
    if nGroup > 1
        offsets = ((1:nGroup) - mean(1:nGroup)) * step;
        xoffset(groupStart:groupEnd) = offsets;
    end
    groupStart = groupEnd + 1;
end
xplot = x;
xplot(order) = xplot(order) + xoffset;
scatter(ax, xplot, y, markerSize, 'k', 'filled', ...
    'MarkerFaceAlpha',0.8, 'MarkerEdgeAlpha',0.8);
plot(ax, [0.8 1.2],[0.5 0.5],':','Color',[.65 .65 .65],'LineWidth',2);
if ~isempty(y)
    plot(ax, [0.85 1.15],[median(y) median(y)], 'k-','LineWidth',3);
end
xlim(ax,[0.7 1.3]); ylim(ax,[0 1]);
xticks(ax,[1]); yticks(ax,[0 0.5 1]);
xticklabels(ax,[]); yticklabels(ax,[]);
box(ax,'off')

%% ORIGINAL histogram style preserved; fish medians added as overlay only
for ii = 1:2
    figure(2+ii)
    STdelay_cropped = max(min(STdelay(:,ii),.5),-.5);

    plot([0 0],[0 0.5],'color',[.65 .65 .65],'linestyle',':','linewidth',5); hold on;
    plot([median(STdelay(:,ii)) median(STdelay(:,ii))],[0 0.5],'color','k','linestyle',':','linewidth',5); hold on

    histogram(STdelay_cropped,[-0.5:0.1:0.5], ...
        'normalization','probability', ...
        'displaystyle','stairs', ...
        'LineWidth',8, ...
        'EdgeColor',plotcolor ); hold on;

    fishMedian = nan(nFish,1);
    for ff = 1:nFish
        if ~isempty(STdelay_fish{ff,ii})
            fishMedian(ff) = median(STdelay_fish{ff,ii});
        end
    end
    validFish = ~isnan(fishMedian);
    ydot = 0.455 + 0.02*rand(sum(validFish),1);
    scatter(fishMedian(validFish), ydot, 25, 'k', 'filled', ...
        'MarkerFaceAlpha',0.7, 'MarkerEdgeAlpha',0.7); hold on;

    xlim([-0.5 0.5]); ylim([0 0.5]);
    xticks([-0.5:0.1:0.5]); yticks([0:0.1:0.5]);
    yticklabels([]); xticklabels([]); box off;
    set(gcf,'Position',[500 1400-ii*300 600 200])

    format short
    % disp([min(STdelay(:,ii)) max(STdelay(:,ii)) median(STdelay(:,ii)) var(STdelay(:,ii))])
end

%% event-level binomial GLME for S/T/ST ipsi fractions
tbl_ST = table(Y_ST_all, categorical(FishID_ST_all), ...
    'VariableNames', {'Y','FishID'});
if height(tbl_ST) > 1
    glme_ST = fitglme(tbl_ST, 'Y ~ 1 + (1|FishID)', ...
        'Distribution','Binomial', 'Link','logit');
    beta0_ST = fixedEffects(glme_ST);
    p_ST = coefTest(glme_ST, [1]);
    fprintf('GLME ST ipsi fraction: logit=%.3f, p=%.4g\n', beta0_ST, p_ST);
end
%% event-level mixed model for delays
tbl_delay = table(DelayAll, categorical(FishIDAll), TypeAll, ...
    'VariableNames', {'Delay','FishID','DelayType'});

for kk = 1:2
    idx = (tbl_delay.DelayType == kk);
    tbl_k = tbl_delay(idx,:);
    if height(tbl_k) < 2
        fprintf('LME delay type %d: not enough events, skipping\n', kk);
        continue;
    end
    lme_k = fitlme(tbl_k, 'Delay ~ 1 + (1|FishID)');
    beta0 = fixedEffects(lme_k);
    [p_F, Fstat, df1, df2] = coefTest(lme_k, [1]);
    fprintf('LME delay type %d: mean=%.4f, F=%.3f, df=[%.1f %.1f], p=%.4g\n', ...
        kk, beta0, Fstat, df1, df2, p_F);
end

