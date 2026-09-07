clear all; close all; 
global gh
LoadFishNColorSel;
addpath('./stattool'); addpath('./disptool');

nbin      = 8;      % 40 s total in 5 s bins
alpha_CI  = 0.05;   % 95% CI
bandAlpha = 0.15;   % shading alpha

% containers for time-course plotting (trial-based)
for ii = 1:3
    datamtx{ii} = [];
end

% containers for mixed model / fish-based (long format)
all_y      = [];
all_bin    = [];
all_group  = [];
all_fishID = [];

for fishsub = 1:length(totalfishsub)
    gh.param.fishid = totalfishsub(fishsub);
    fm_behavim_main;

    for ii = 1:sessionn(fishsub)
        saccaden = [];

        if ~ismember(ii, gh.param.ExcludedSession)

            % --- bin counts (saccades), 0 to 40 s in 5 s bins ---
            for jj = 1:nbin
                t1 = (jj-1)*5;
                t2 = jj*5;
                saccaden(jj) = length(find(gh.data.saccademtx(:,1)==ii & ...
                                           gh.data.saccademtx(:,2) > t1 & ...
                                           gh.data.saccademtx(:,2) <= t2));
            end

            % --- minus baseline: first 10 s = bins 1–2 ---
            n_saccaden = saccaden - mean(saccaden(1:2));

            % --- group code ---
            trcode = gh.param.fishlog.trialdetails.trial(ii,1);
            if     trcode == 8
                g = 1;              % blank
            elseif trcode == 1 || trcode == 2
                g = 2;              % appetitive
            elseif trcode == 4 || trcode == 5
                g = 3;              % aversive
            else
                continue
            end

            % --- matrix for trial-based plotting ---
            datamtx{g} = [datamtx{g}; n_saccaden./5];

            % --- long format for mixed model / fish-based plot ---
            trialID = ii;
            for jj = 1:nbin
                all_y      = [all_y;      n_saccaden(jj)/5];
                all_bin    = [all_bin;    jj];
                all_group  = [all_group;  g];
                all_fishID = [all_fishID; gh.param.fishid];
            end

            clear n_saccaden
        end

        clear saccaden
    end
end

%% ---------- colors & styles ----------
linecolor{1} = blankColor;
linecolor{2} = appeColor;
linecolor{3} = averColor;

linestyle{1} = '--';
linestyle{2} = ':';
linestyle{3} = '-';

%% ---------- fish-based time-course plot (mean ± 95% CI across fish) ----------
Tall = table(all_y, categorical(all_bin), ...
             categorical(all_group), categorical(all_fishID), ...
             'VariableNames', {'y','bin','group','fishID'});

% collapse events/trials within each fish for each group × bin
[G1, fish_u, group_u, bin_u] = findgroups(Tall.fishID, Tall.group, Tall.bin);
y_fishbin = splitapply(@mean, Tall.y, G1);

figure(2); clf; hold on;
set(gcf,'Position',[920 800 130 250]);

for g = 1:3
    idx_g      = (double(group_u) == g);
    fish_ids_g = unique(fish_u(idx_g));

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
    m = mean(tempmtx, 1, 'omitnan');

    % n fish per bin
    nfish = sum(~isnan(tempmtx), 1);

    % SEM across fish
    se = std(tempmtx, 0, 1, 'omitnan') ./ sqrt(nfish);

    % 95% CI t-multiplier per bin (df = nfish-1)
    tcrit_fish = nan(1, nbin);
    for jj = 1:nbin
        if nfish(jj) > 1
            tcrit_fish(jj) = tinv(1 - alpha_CI/2, nfish(jj)-1);
        end
    end

    ci_upper = m + tcrit_fish .* se;
    ci_lower = m - tcrit_fish .* se;

    valid   = ~isnan(ci_upper) & ~isnan(ci_lower) & ~isnan(m);
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

    % mean curve
    plot(1:nbin, m, ...
         'Color',     linecolor{g}, ...
         'LineWidth', 1.5, ...
         'LineStyle', linestyle{g});
end

ylim([-0.14 0.14]);
yticks(-0.14:0.035:0.14);
xlim([0.5 nbin+0.5]);
xticks(1:nbin);
plot([0.5 nbin+0.5],[0 0],'k','linestyle','--');

xticklabels({});
yticklabels({});
box off;

yl = ylim;
pgon = polyshape([4.5 6.5 6.5 4.5], ...
                 [yl(1) yl(1) yl(2) yl(2)]);
plot(pgon,'FaceColor','none','EdgeColor','k','linestyle','--');

h = gca;
h.XAxis.Visible = 'off';

%% ---------- mixed model on stimulus bins 5–7 ----------
idx_stim = ismember(all_bin, [5 6 7]);

T = table(all_y(idx_stim), ...
          categorical(all_group(idx_stim)), ...
          categorical(all_fishID(idx_stim)), ...
          categorical(all_bin(idx_stim)), ...
          'VariableNames', {'y','group','fishID','bin'});

lme = fitlme(T, 'y ~ 1 + group + (1|fishID)');
anova(lme)

[beta, SE, stats] = fixedEffects(lme,'DFMethod','Residual');
stats

H = [0 1 -1];
[p_c, F_c, DF1_c, DF2_c] = coefTest(lme, H);

beta_2 = stats.Estimate(2);
beta_3 = stats.Estimate(3);
beta_c = beta_2 - beta_3;

p_c, F_c, DF1_c, DF2_c, beta_c