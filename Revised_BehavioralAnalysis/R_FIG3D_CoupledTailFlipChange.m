clear all; close all; 
global gh
LoadFishNColorSel;
addpath('./stattool'); addpath('./disptool');

nbin = 6;
alpha = 0.05;       % 95% CI
bandAlpha = 0.15;   % transparency of CI band

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
        tailn = []; simu = []; saccaden = []; ration = [];

        if ~ismember(ii, gh.param.ExcludedSession)

            % ---------- binning ----------
            for jj = 1:nbin
                t1 = 10 + (jj-1)*30/nbin;
                t2 = 10 +  jj   *30/nbin;

                tailn(jj) = length(find(gh.data.bout_details(:,1)==ii & ...
                                        gh.data.bout_details(:,2) > t1 & ...
                                        gh.data.bout_details(:,2) <= t2));

                simu(jj)  = length(find(gh.data.simuMtx(:,1)==ii & ...
                                        gh.data.simuMtx(:,9) > t1 & ...
                                        gh.data.simuMtx(:,9) <= t2));

                saccaden(jj) = length(find(gh.data.saccademtx(:,1)==ii & ...
                                           gh.data.saccademtx(:,2) > t1 & ...
                                           gh.data.saccademtx(:,2) <= t2));

                ration(jj) = simu(jj) - (tailn(jj) - simu(jj));
            end

            % ---------- baseline subtraction (pre-stim bins 1–2) ----------
            n_ration = ration - mean(ration(1:(nbin/3)));

            % ---------- group code ----------
            trcode = gh.param.fishlog.trialdetails.trial(ii,1);
            if     trcode == 8
                g = 1;                 % blank
            elseif trcode == 1 || trcode == 2
                g = 2;                 % appetitive
            elseif trcode == 4 || trcode == 5
                g = 3;                 % aversive
            else
                continue
            end

            % ---------- trial-based matrix ----------
            datamtx{g} = [datamtx{g}; n_ration/5];

            % ---------- long-format for fish-based plot & mixed model ----------
            for jj = 1:nbin
                all_y      = [all_y;      n_ration(jj)/5];
                all_bin    = [all_bin;    jj];
                all_group  = [all_group;  g];
                all_fishID = [all_fishID; gh.param.fishid];
            end

            clear n_ration
        end

        clear tailn simu saccaden ration
    end
end

%% ---------- style ----------
linecolor{1}  = blankColor;
linecolor{2}  = appeColor;
linecolor{3}  = averColor;

linestyle{1}  = '--';
linestyle{2}  = ':';
linestyle{3}  = '-';

%% ---------- Figure 1: fish-based time-course (mean ± 95% CI across larvae) ----------
Tall = table(all_y, categorical(all_bin), ...
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
    m = mean(tempmtx, 1, 'omitnan');

    % n fish per bin
    nfish = sum(~isnan(tempmtx), 1);

    % SEM across fish
    se = std(tempmtx, 0, 1, 'omitnan') ./ sqrt(nfish);

    % 95% CI using t-multiplier matched to nfish at each bin
    tcrit_fish = nan(1, nbin);
    for jj = 1:nbin
        if nfish(jj) > 1
            tcrit_fish(jj) = tinv(1 - alpha/2, nfish(jj)-1);
        end
    end

    ci_upper = m + tcrit_fish .* se;
    ci_lower = m - tcrit_fish .* se;

    valid = ~isnan(ci_upper) & ~isnan(ci_lower) & ~isnan(m);
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
         'LineStyle', linestyle{g}); hold on
end

ylim([-0.25 0.25]); yticks(-0.25:0.0625:0.25);
xlim([0.5 nbin+0.5]); xticks(1:nbin);
plot([0.5 nbin+0.5],[0 0],'k','LineStyle','--');

xticklabels({}); 
yticklabels({});
box off;

yl = ylim;
pgon = polyshape([nbin/3+0.5 nbin*2/3+0.5 nbin*2/3+0.5 nbin/3+0.5], ...
                 [yl(1) yl(1) yl(2) yl(2)]);
plot(pgon,'FaceColor','none','EdgeColor','k','LineStyle','--');

h = gca; 
h.XAxis.Visible = 'off';

%% ---------- Mixed model on bins 3–5 (stim + 5 s post, fish random) ----------
idx_stim = ismember(all_bin, [3 4 5]);

Tmm = table(all_y(idx_stim), ...
            categorical(all_group(idx_stim)), ...
            categorical(all_fishID(idx_stim)), ...
            categorical(all_bin(idx_stim)), ...
            'VariableNames', {'y','group','fishID','bin'});

lme_ratio = fitlme(Tmm, 'y ~ 1 + group + (1|fishID)');

disp('=== Mixed model (S–T vs independent ratio, stim + 5 s post, bins 3–5) ===');
anova(lme_ratio);

[beta, SE, stats] = fixedEffects(lme_ratio,'DFMethod','Residual');
stats

% appetitive vs aversive contrast
H = [0 1 -1];
[p_c, F_c, DF1_c, DF2_c] = coefTest(lme_ratio, H);

fprintf('Time-course mixed model, av vs app: F=%.3f, p=%.4f\n', ...
        F_c, p_c);