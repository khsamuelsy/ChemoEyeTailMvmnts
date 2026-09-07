clear all; close all; global gh
LoadFishNColorSel;rng(1);

count_S=0; count_T=0; count_ST=0;
count_S_contra = 0; count_S_ipsi = 0; 
count_T_contra = 0; count_T_ipsi = 0;
count_ST_contra = 0; count_ST_ipsi = 0; 
STdelay=[];

addpath('./stattool');addpath('./disptool');

nFish = length(totalfishsub);STdelay_fish = cell(nFish,2);

S_ipsi_fish   = zeros(nFish,1);S_contra_fish = zeros(nFish,1);
T_ipsi_fish   = zeros(nFish,1);T_contra_fish = zeros(nFish,1);
ST_ipsi_fish  = zeros(nFish,1);ST_contra_fish= zeros(nFish,1);

DelayAll = [];
FishIDAll = [];
TypeAll = [];

Y_S_all = [];  FishID_S_all = [];
Y_T_all = [];  FishID_T_all = [];
Y_ST_all = []; FishID_ST_all = [];

gh.cum_saccaden=0;
for fishsub = 1:nFish
    gh.param.fishid = totalfishsub(fishsub);
    fm_behavim_main

    S_ipsi_this = 0; S_contra_this = 0;
    T_ipsi_this = 0; T_contra_this = 0;
    ST_ipsi_this = 0; ST_contra_this = 0;

    for ii=1:sessionn(fishsub)
        if ~ismember(ii,gh.param.ExcludedSession)

            if gh.param.fishlog.trialdetails.trial(ii,1)==8
                totalbout_n    = find(gh.data.boutmtx(:,1)==ii    & gh.data.boutmtx(:,4)<=40);
                totalsaccade_n = find(gh.data.saccademtx(:,1)==ii & gh.data.saccademtx(:,4)<=40);
                totalsimu_n    = find(gh.data.simuMtx(:,1)==ii    & gh.data.simuMtx(:,4)<=40);
            else
                totalbout_n    = find(gh.data.boutmtx(:,1)==ii    & gh.data.boutmtx(:,4)<=20);
                totalsaccade_n = find(gh.data.saccademtx(:,1)==ii & gh.data.saccademtx(:,4)<=20);
                totalsimu_n    = find(gh.data.simuMtx(:,1)==ii    & gh.data.simuMtx(:,4)<=20);
            end

            count_T  = count_T  + length(totalbout_n);
            count_S  = count_S  + length(totalsaccade_n);
            count_ST = count_ST + length(totalsimu_n);

            %% S events
            for saccaden = 2:length(totalsaccade_n)
                if gh.data.saccademtx(totalsaccade_n(saccaden),7) == ...
                   gh.data.saccademtx(totalsaccade_n(saccaden-1),7)
                    count_S_ipsi = count_S_ipsi + 1;
                    S_ipsi_this  = S_ipsi_this + 1;
                    Y_S_all      = [Y_S_all; 1];
                else
                    count_S_contra = count_S_contra + 1;
                    S_contra_this  = S_contra_this + 1;
                    Y_S_all        = [Y_S_all; 0];
                end
                FishID_S_all = [FishID_S_all; fishsub];
            end

            %% T events
            saved_lastboutdir = [];
            threshold_boutbias = 1.5;
            for boutn = 1:length(totalbout_n)
                boutdir  = sign(gh.data.bout_details(totalbout_n(boutn),7)*-1 - anglebias_overall(fishsub));
                boutbias = abs(gh.data.bout_details(totalbout_n(boutn),7)*-1 - anglebias_overall(fishsub));

                if boutbias >= threshold_boutbias
                    if boutn > 1 && ~isempty(saved_lastboutdir)
                        if boutdir == saved_lastboutdir
                            count_T_ipsi = count_T_ipsi + 1;
                            T_ipsi_this  = T_ipsi_this + 1;
                            Y_T_all      = [Y_T_all; 1];
                        else
                            count_T_contra = count_T_contra + 1;
                            T_contra_this  = T_contra_this + 1;
                            Y_T_all        = [Y_T_all; 0];
                        end
                        FishID_T_all = [FishID_T_all; fishsub];
                    end
                    saved_lastboutdir = boutdir;
                else
                    saved_lastboutdir = [];
                end
            end

            %% ST events + delays
            for simun = 1:length(totalsimu_n)
                boutdir = sign(gh.data.simuMtx(totalsimu_n(simun),14)*-1 - anglebias_overall(fishsub));
                boutbias = abs(gh.data.simuMtx(totalsimu_n(simun),14)*-1 - anglebias_overall(fishsub));
                saccadedir = gh.data.simuMtx(totalsimu_n(simun),7);

                delayvec = [gh.data.simuMtx(totalsimu_n(simun),9)  - gh.data.simuMtx(totalsimu_n(simun),2), ... %onset delay
                    gh.data.simuMtx(totalsimu_n(simun),4)  - gh.data.simuMtx(totalsimu_n(simun),9)]; ... %overlap duration

                STdelay = [STdelay; delayvec];

                for kk = 1:2
                    STdelay_fish{fishsub,kk} = [STdelay_fish{fishsub,kk}; delayvec(kk)];
                end

                for kk = 1:2
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

    S_ipsi_fish(fishsub)    = S_ipsi_this;
    S_contra_fish(fishsub)  = S_contra_this;
    T_ipsi_fish(fishsub)    = T_ipsi_this;
    T_contra_fish(fishsub)  = T_contra_this;
    ST_ipsi_fish(fishsub)   = ST_ipsi_this;
    ST_contra_fish(fishsub) = ST_contra_this;
end

%% fish-level fractions
S_frac_fish  = S_ipsi_fish  ./ (S_ipsi_fish  + S_contra_fish);
T_frac_fish  = T_ipsi_fish  ./ (T_ipsi_fish  + T_contra_fish);
ST_frac_fish = ST_ipsi_fish ./ (ST_ipsi_fish + ST_contra_fish);

%% helper settings for side plots
minSep = 0.06;
step   = 0.04;
markerSize = 28;

%% Figure 1: ST pie only
figure(1); clf
set(gcf,'Position',[100 1000 300 300])
ax = axes('Position',[0.08 0.08 0.84 0.84]);
pieData = [count_ST_contra count_ST_ipsi];
h = pie(ax, pieData);
ax.Colormap = [.65 .65 .65; 0.9 0.9 0.9];
axis(ax,'off')
delete(findobj(h,'Type','text'));

%% Figure 2: S pie only
figure(2); clf
set(gcf,'Position',[420 1000 300 300])
ax = axes('Position',[0.08 0.08 0.84 0.84]);
pieData = [count_S_contra count_S_ipsi];
h = pie(ax, pieData);
ax.Colormap = [.65 .65 .65; 0.9 0.9 0.9];
axis(ax,'off')
delete(findobj(h,'Type','text'));

%% Figure 3: T pie only
figure(3); clf
set(gcf,'Position',[740 1000 300 300])
ax = axes('Position',[0.08 0.08 0.84 0.84]);
pieData = [count_T_contra count_T_ipsi];
h = pie(ax, pieData);
ax.Colormap = [.65 .65 .65; 0.9 0.9 0.9];
axis(ax,'off')
delete(findobj(h,'Type','text'));

%% Figure 4: ST fish-level fraction only
figure(4); clf
set(gcf,'Position',[100 600 220 300])
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
plot(ax, [0.85 1.15],[median(y) median(y)], 'k-','LineWidth',3);
xlim(ax,[0.7 1.3]); ylim(ax,[0 1]);
xticks(ax,[1]); yticks(ax,[0 0.5 1]);
xticklabels(ax,[]); yticklabels(ax,[]);
box(ax,'off')

%% Figure 5: S fish-level fraction only
figure(5); clf
set(gcf,'Position',[360 600 220 300])
ax = axes('Position',[0.25 0.15 0.55 0.72]); hold(ax,'on')
valid = ~isnan(S_frac_fish) & isfinite(S_frac_fish);
y = S_frac_fish(valid);
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
plot(ax, [0.85 1.15],[median(y) median(y)], 'k-','LineWidth',3);
xlim(ax,[0.7 1.3]); ylim(ax,[0 1]);
xticks(ax,[1]); yticks(ax,[0 0.5 1]);
xticklabels(ax,[]); yticklabels(ax,[]);
box(ax,'off')

%% Figure 6: T fish-level fraction only
figure(6); clf
set(gcf,'Position',[620 600 220 300])
ax = axes('Position',[0.25 0.15 0.55 0.72]); hold(ax,'on')
valid = ~isnan(T_frac_fish) & isfinite(T_frac_fish);
y = T_frac_fish(valid);
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
plot(ax, [0.85 1.15],[median(y) median(y)], 'k-','LineWidth',3);
xlim(ax,[0.7 1.3]); ylim(ax,[0 1]);
xticks(ax,[1]); yticks(ax,[0 0.5 1]);
xticklabels(ax,[]); yticklabels(ax,[]);
box(ax,'off')

%% Figures 7-8: pooled histograms + fish-level medians
for ii=1:2
    figure(6+ii); 
    STdelay_cropped = max(min(STdelay(:,ii),.5),-.5);

    plot([0 0],[0 0.25], ...
        'color',[.65 .65 .65], ...
        'linestyle',':', ...
        'linewidth',5); hold on;

    plot([median(STdelay(:,ii)) median(STdelay(:,ii))],[0 0.25], ...
        'color','k', ...
        'linestyle',':', ...
        'linewidth',5); hold on

    histogram(STdelay_cropped, -0.5:0.05:0.5, ...
        'normalization','probability', ...
        'displaystyle','stairs', ...
        'LineWidth',8, ...
        'EdgeColor',STColor); hold on;

    fishMedian = nan(nFish,1);
    for fishsub = 1:nFish
        if ~isempty(STdelay_fish{fishsub,ii})
            fishMedian(fishsub) = median(STdelay_fish{fishsub,ii});
        end
    end

    validFish = ~isnan(fishMedian);
    ydot = 0.228 + 0.012*rand(sum(validFish),1);
    scatter(fishMedian(validFish), ydot, ...
        25, 'k', 'filled', ...
        'MarkerFaceAlpha',0.7, ...
        'MarkerEdgeAlpha',0.7); hold on;

    xlim([-0.5 0.5]); ylim([0 0.25]);
    xticks(-0.5:0.1:0.5); yticks(0:0.05:0.25);
    xticklabels([]); yticklabels([]);
    box off
    set(gcf,'Position',[1200 1300-ii*300 600 200])
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

tbl_S = table(Y_S_all, categorical(FishID_S_all), ...
    'VariableNames', {'Y','FishID'});
if height(tbl_S) > 1
    glme_S = fitglme(tbl_S, 'Y ~ 1 + (1|FishID)', ...
        'Distribution','Binomial', 'Link','logit');
    beta0_S = fixedEffects(glme_S);
    p_S = coefTest(glme_S, [1]);
    fprintf('GLME S ipsi fraction: logit=%.3f, p=%.4g\n', beta0_S, p_S);
end

tbl_T = table(Y_T_all, categorical(FishID_T_all), ...
    'VariableNames', {'Y','FishID'});
if height(tbl_T) > 1
    glme_T = fitglme(tbl_T, 'Y ~ 1 + (1|FishID)', ...
        'Distribution','Binomial', 'Link','logit');
    beta0_T = fixedEffects(glme_T);
    p_T = coefTest(glme_T, [1]);
    fprintf('GLME T ipsi fraction: logit=%.3f, p=%.4g\n', beta0_T, p_T);
end