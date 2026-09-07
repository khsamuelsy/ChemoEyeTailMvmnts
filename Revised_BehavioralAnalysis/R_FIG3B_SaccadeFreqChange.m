clear all; close all; global gh
LoadFishNColorSel;
addpath('./stattool');addpath('./disptool');

nbin = 6;

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
        saccaden = [];
        if ~ismember(ii, gh.param.ExcludedSession)

            % --- bin counts (saccades) ---
            for jj = 1:nbin
                t1 = 10 + (jj-1)*30/nbin;
                t2 = 10 +  jj   *30/nbin;
                saccaden(jj) = length(find(gh.data.saccademtx(:,1)==ii & ...
                                           gh.data.saccademtx(:,2) > t1 & ...
                                           gh.data.saccademtx(:,2) <= t2));
            end

            % --- minus prestim: bins 1–2 ---
            n_saccaden = saccaden - mean(saccaden(1:2));

            % --- group code ---
            trcode = gh.param.fishlog.trialdetails.trial(ii,1);
            if trcode == 8
                group(ii) = 1;              % blank
            elseif trcode == 1 || trcode == 2
                group(ii) = 2;              % appetitive
            elseif trcode == 4 || trcode == 5
                group(ii) = 3;              % aversive
            else
                continue
            end

            % --- matrix for trial-based plotting ---
            datamtx{group(ii)} = [datamtx{group(ii)}; n_saccaden./5];

            % --- long format for mixed model / fish-based plot ---
            trialID = ii;
            for jj = 1:nbin
                all_y      = [all_y;      n_saccaden(jj)/5];
                all_bin    = [all_bin;    jj];
                all_group  = [all_group;  group(ii)];
                all_fishID = [all_fishID; gh.param.fishid];
            end

            clear n_saccaden
        end
        clear saccaden
    end
end

%% ---------- plotting style ----------
linecolor{1} = blankColor;
linecolor{2} = appeColor;
linecolor{3} = averColor;
linestyle{1} = '--';
linestyle{2} = ':';
linestyle{3} = '-';

%% ---------- fish-based time-course plot ----------
Tall = table(all_y, categorical(all_bin), ...
             categorical(all_group), categorical(all_fishID), ...
             'VariableNames', {'y','bin','group','fishID'});

% collapse events/trials within each fish for each group x bin
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

    % SEM across fish
    nfish = sum(~isnan(tempmtx), 1);
    se = std(tempmtx, 0, 1, 'omitnan') ./ sqrt(nfish);

    % 95% CI across fish: use t critical value reflecting per-bin n
    tcrit = nan(1, nbin);
    for b = 1:nbin
        if nfish(b) > 1
            tcrit(b) = tinv(0.975, nfish(b)-1);
        else
            tcrit(b) = NaN;
        end
    end
    ci = tcrit .* se;

    upper = m + ci;
    lower = m - ci;

    % CI shading
    pgon = polyshape([1:nbin, fliplr(1:nbin)], ...
                     [upper, fliplr(lower)]);
    plot(pgon, ...
         'FaceColor', linecolor{g}, ...
         'FaceAlpha', 0.15, ...
         'EdgeColor', 'none');

    % mean curve
    plot(1:nbin, m, ...
         'Color', linecolor{g}, ...
         'LineWidth', 1.5, ...
         'LineStyle', linestyle{g});
end

% axis style
ylim([-0.14 0.14]);
yticks(-0.14:0.035:0.14);
xlim([0.5 nbin+0.5]);
xticks(1:nbin);
plot([0.5 nbin+0.5],[0 0],'k','linestyle','--');
xticklabels({});
yticklabels({});
box off;

yl = ylim;
pgon = polyshape([nbin/3+0.5 nbin*2/3+0.5 nbin*2/3+0.5 nbin/3+0.5], ...
                 [yl(1) yl(1) yl(2) yl(2)]);
plot(pgon,'FaceColor','none','EdgeColor','k','linestyle','--');

h = gca;
h.XAxis.Visible = 'off';

%% ---------- mixed model on bins 3–5 (trial + fish in one model) ----------
idx_stim = ismember(all_bin, [3 4 5]);
T = table(all_y(idx_stim), ...
          categorical(all_group(idx_stim)), ...
          categorical(all_fishID(idx_stim)), ...
          categorical(all_bin(idx_stim)), ...
          'VariableNames', {'y','group','fishID','bin'});

% simplest: same structure as before, just more bins
lme = fitlme(T, 'y ~ 1 + group + (1|fishID)');
anova(lme)                                      % overall group effect
[beta, SE, stats] = fixedEffects(lme,'DFMethod','Residual');
stats                                            % group_2, group_3 vs blank

H = [0 1 -1];
[p_c, F_c, DF1_c, DF2_c] = coefTest(lme, H);

fprintf('Time-course mixed model (bins 3–5), av vs app: F=%.3f, p=%.4f\n', ...
        F_c, p_c);