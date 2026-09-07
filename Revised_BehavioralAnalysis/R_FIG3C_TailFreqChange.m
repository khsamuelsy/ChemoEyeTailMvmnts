clear all; close all;
global gh
LoadFishNColorSel;
addpath('./stattool'); addpath('./disptool');

nbin = 18;
alpha = 0.05;       % 95% CI
bandAlpha = 0.15;   % error band transparency

% containers for time-course plotting (trial-based)
for ii = 1:3
    datamtx{ii} = [];
end

% containers for mixed model / fish-based (long format)
all_y      = [];all_bin    = [];
all_group  = [];all_fishID = [];

for fishsub = 1:length(totalfishsub)
    gh.param.fishid = totalfishsub(fishsub);
    fm_behavim_main;

    for ii = 1:sessionn(fishsub)
        tailn = [];

        if ~ismember(ii, gh.param.ExcludedSession)

            % ---------- binning ----------
            for jj = 1:nbin
                t1 = 10 + (jj-1)*30/nbin;
                t2 = 10 +  jj   *30/nbin;

                tailn(jj) = length(find(gh.data.bout_details(:,1)==ii & ...
                                        gh.data.bout_details(:,2) > t1 & ...
                                        gh.data.bout_details(:,2) <= t2));
            end

            % ---------- z-prestim & baseline subtraction ----------
            varval = std(tailn(1:(nbin/3)));
            if varval == 0
                continue
            end

            n1_tailn = tailn ./ varval;
            n2_tailn = n1_tailn - mean(n1_tailn(1:(nbin/3)));

            % ---------- group code ----------
            trcode = gh.param.fishlog.trialdetails.trial(ii,1);
            if trcode == 8
                g = 1;
            elseif trcode == 1 || trcode == 2
                g = 2;
            elseif trcode == 4 || trcode == 5
                g = 3;
            else
                continue
            end

            % ---------- trial-based matrix ----------
            datamtx{g} = [datamtx{g}; n2_tailn];

            % ---------- long-format for fish-based plot & mixed model ----------
            for jj = 1:nbin
                all_y      = [all_y;      n2_tailn(jj)];
                all_bin    = [all_bin;    jj];
                all_group  = [all_group;  g];
                all_fishID = [all_fishID; gh.param.fishid];
            end

            clear n1_tailn n2_tailn
        end

        clear tailn 
    end
end

%% ---------- line styles ----------
linecolor{1} = blankColor;
linecolor{2} = appeColor;
linecolor{3} = averColor;

linestyle{1} = '--';
linestyle{2} = ':';
linestyle{3} = '-';


%% ---------- Figure 1: fish-based time-course (mean ± 95% CI across fish) ----------
Tall = table(all_y, all_bin, ...
             categorical(all_group), categorical(all_fishID), ...
             'VariableNames', {'y','bin','group','fishID'});

% collapse within fish × group × bin (fish-level means)
[G1, fish_u, group_u, bin_u] = findgroups(Tall.fishID, Tall.group, Tall.bin);
y_fishbin = splitapply(@mean, Tall.y, G1);

figure(1); clf; hold on;
set(gcf,'Position',[920 800 100 250]);

for g = 1:3
    idx_g = (double(group_u) == g);
    fish_ids_g = unique(fish_u(idx_g));

    % rows = fish, cols = bins
    tempmtx = nan(numel(fish_ids_g), nbin);

    for f = 1:numel(fish_ids_g)
        idx_fg = idx_g & (fish_u == fish_ids_g(f));
        bins_f = double(bin_u(idx_fg));
        y_f    = y_fishbin(idx_fg);

        [bins_f, ord] = sort(bins_f);
        y_f = y_f(ord);

        tempmtx(f, bins_f) = y_f;
    end

    % mean across fish
    m_fish = mean(tempmtx,1,'omitnan');

    % n fish per bin
    nfish = sum(~isnan(tempmtx),1);

    % SEM across fish
    se_fish = std(tempmtx,0,1,'omitnan') ./ sqrt(nfish);

    % 95% CI using t-multiplier matched to nfish at each bin
    tcrit_bin = nan(1, nbin);
    for jj = 1:nbin
        if nfish(jj) > 1
            tcrit_bin(jj) = tinv(1 - alpha/2, nfish(jj)-1);
        end
    end

    ci_upper = m_fish + tcrit_bin .* se_fish;
    ci_lower = m_fish - tcrit_bin .* se_fish;

    valid = ~isnan(ci_upper) & ~isnan(ci_lower) & ~isnan(m_fish);
    x_valid = find(valid);

    % 95% CI shading
    if ~isempty(x_valid)
        pgon = polyshape([x_valid, fliplr(x_valid)], ...
                         [ci_upper(valid), fliplr(ci_lower(valid))]);
        plot(pgon, ...
             'FaceColor', linecolor{g}, ...
             'FaceAlpha', bandAlpha, ...
             'EdgeColor', 'none');
    end

    % smoothed mean across fish
    m_fish_smooth = smooth(m_fish,3);
    plot(1:nbin, m_fish_smooth, ...
         'Color',     linecolor{g}, ...
         'LineWidth', 1.5, ...
         'LineStyle', linestyle{g}); hold on

    % raw mean across fish
    plot(1:nbin, m_fish, ...
         'Color',     linecolor{g}, ...
         'LineWidth', 1, ...
         'LineStyle', linestyle{g}); hold on
end

xlim([0.5 nbin+0.5]);
ylim([-1.6 1.6]);
yticks(-1.6:0.4:1.6);
plot([0.5 nbin+0.5],[0 0],'k','LineStyle','--');

xticks(1:nbin);
xticklabels({});
yticklabels({});
box off;

yl = ylim;
pgon = polyshape([nbin/3+0.5 nbin*2/3+0.5 nbin*2/3+0.5 nbin/3+0.5], ...
                 [yl(1) yl(1) yl(2) yl(2)]);
plot(pgon,'FaceColor','none','EdgeColor','k','LineStyle','--');

h = gca;
h.XAxis.Visible = 'off';

%% ---------- mixed model: stim + 5 s post (bins 7–15, fish random) ----------
idx_stim_post = ismember(Tall.bin, 7:15);
Tstimpost = Tall(idx_stim_post, :);

lme_tail_stimpost = fitlme(Tstimpost, 'y ~ 1 + group + (1|fishID)');

disp('=== Mixed model (tail metric, stim + 5 s post, bins 7–15, fish random) ===')
anova(lme_tail_stimpost)

[beta_sp, SE_sp, stats_sp] = fixedEffects(lme_tail_stimpost,'DFMethod','Residual');
disp(stats_sp)

% 95% CI for fixed effects using model DF
tcrit_fix = tinv(1 - alpha/2, stats_sp.DF);
CI_low_sp  = stats_sp.Estimate - tcrit_fix .* stats_sp.SE;
CI_high_sp = stats_sp.Estimate + tcrit_fix .* stats_sp.SE;

disp(table(stats_sp.Name, stats_sp.Estimate, stats_sp.SE, CI_low_sp, CI_high_sp, stats_sp.pValue, ...
    'VariableNames', {'Name','Beta','SE','CI_low','CI_high','pValue'}))

% group 2 vs group 3 over stim + post window
H = [0 1 -1];
[p_c_sp, F_c_sp, DF1_c_sp, DF2_c_sp] = coefTest(lme_tail_stimpost, H);

disp('Group 2 vs Group 3 (tail metric, stim + 5 s post, bins 7–15):')
disp(['p = ' num2str(p_c_sp)])
disp(['F = ' num2str(F_c_sp)])
disp(['DF1 = ' num2str(DF1_c_sp) ', DF2 = ' num2str(DF2_c_sp)])
