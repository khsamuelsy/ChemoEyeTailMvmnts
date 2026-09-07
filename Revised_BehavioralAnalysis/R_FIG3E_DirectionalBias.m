clear all; close all;
global gh
LoadFishNColorSel;
addpath('./stattool'); addpath('./disptool');

%% ---------- collect bout-wise data ----------
all_y_bout       = [];   % bout-level absolute bias-corrected angle
all_group3_bout  = [];   % 1=blank, 2=app pooled, 3=av pooled
all_group5_bout  = [];   % 1=blank, 2=app1, 3=app2, 4=av1, 5=av2
all_fishID_bout  = [];


for fishsub = 1:length(totalfishsub)
    gh.param.fishid = totalfishsub(fishsub);
    fm_behavim_main
    for ii = 1:sessionn(fishsub)
        if ~ismember(ii, gh.param.ExcludedSession)

            % bouts during 20–30 s window
            stim_boutn = find(gh.data.boutmtx(:,1) == ii & ...
                              gh.data.boutmtx(:,2) > 20 & ...
                              gh.data.boutmtx(:,4) <= 30);
            prestim_boutn = find(gh.data.boutmtx(:,1) == ii & ...
                                 gh.data.boutmtx(:,2) > 10 & ...
                                 gh.data.boutmtx(:,4) <= 20);

            if isempty(stim_boutn)
                continue
            end

            trcode = gh.param.fishlog.trialdetails.trial(ii,1);

            % 5-group labelling (chemicals)
            if     trcode == 8
                grp5 = 1;           % blank
            elseif trcode == 1
                grp5 = 2;           % appetitive 1
            elseif trcode == 2
                grp5 = 3;           % appetitive 2
            elseif trcode == 4
                grp5 = 4;           % aversive 1
            elseif trcode == 5
                grp5 = 5;           % aversive 2
            else
                continue;
            end

            % 3-group pooling (blank / app / av)
            if     trcode == 8
                grp3 = 1;           % blank
            elseif trcode == 1 || trcode == 2
                grp3 = 2;           % appetitive
            elseif trcode == 4 || trcode == 5
                grp3 = 3;           % aversive
            else
                continue;
            end

            % collect bout-wise values
            sumval = 0;
            if length(prestim_boutn) > 0
                for jj = 1:length(prestim_boutn)
                    sumval = sumval + abs( gh.data.bout_details(prestim_boutn(jj),7)*(-1) ...
                        - anglebias_overall(fishsub));
                end
                avgval = sumval ./ length(prestim_boutn);
            else
                avgval = 0;
            end

            for jj = 1:length(stim_boutn)
                val = abs( gh.data.bout_details(stim_boutn(jj),7)*(-1) ...
                           - anglebias_overall(fishsub) ) - avgval;

                all_y_bout       = [all_y_bout;       val];
                all_group3_bout  = [all_group3_bout;  grp3];
                all_group5_bout  = [all_group5_bout;  grp5];
                all_fishID_bout  = [all_fishID_bout;  gh.param.fishid];
            end
        end
    end
end

%% ---------- assemble bout-level table ----------
T_bout = table(all_y_bout, categorical(all_group3_bout), ...
               categorical(all_group5_bout), categorical(all_fishID_bout), ...
               'VariableNames', {'y','group3','group5','fishID'});

%% ---------- bout-median per fish × 3-group (RAW, uncapped) ----------
[G3_f, fish_u3, group3_u] = findgroups(T_bout.fishID, T_bout.group3);
y_fish3 = splitapply(@median, T_bout.y, G3_f);  % RAW fish medians
T_fish3 = table(y_fish3, categorical(group3_u), categorical(fish_u3), ...
                'VariableNames', {'y','group','fishID'});

fish_ids   = categories(T_fish3.fishID);
nFish      = numel(fish_ids);
y_fish_mat = nan(nFish, 3);  % col1=blank, col2=app, col3=av

for f = 1:nFish
    thisFish = fish_ids{f};
    idx_f    = T_fish3.fishID == thisFish;
    g_f      = double(T_fish3.group(idx_f));
    y_f      = T_fish3.y(idx_f);   % RAW fish medians
    for k = 1:numel(g_f)
        y_fish_mat(f, g_f(k)) = y_f(k);
    end
end

%% ---------- bout-level distributions per chemical (RAW for stats) ----------
datavec_bout5 = cell(1,5);
for k = 1:5
    idx_k = (double(T_bout.group5) == k);
    datavec_bout5{k} = T_bout.y(idx_k);
end

%% ---------- figure: bout-level violin + bout medians/IQR + per-fish dots/lines ----------
figure(1); clf; hold on;
set(gcf,'Position',[800 800 300 400])

capLow  = -4;
capHigh =  4;

xpos5 = nan(1,5);
xpos5(1) = 1;
xpos5(2) = 2.2;
xpos5(3) = 2.8;
xpos5(4) = 3.7;
xpos5(5) = 4.3;

xpos3 = [1 2.5 4];

dotcolor{1} = blankColor;
dotcolor{2} = appeColor; dotcolor{3} = appeColor;
dotcolor{4} = averColor; dotcolor{5} = averColor;
pool_col{1} = blankColor;
pool_col{2} = appeColor;
pool_col{3} = averColor;

% 1) 5-group violin-style spread on bout data + bout medians/IQR
for ii = 1:5
    xpos = xpos5(ii);
    this_data_raw = datavec_bout5{ii};
    if isempty(this_data_raw)
        continue
    end

    this_data_plot = max(min(this_data_raw, capHigh), capLow);

    violin(xpos, this_data_plot, ...
           'facecolor', dotcolor{ii}, ...
           'scaling', 1.5, ...
           'facealpha', 0.075, ...
           'style', 2); hold on;

    limiter_vec      = capLow:1:capHigh;
    limiter_vec_indx = zeros(size(this_data_plot));
    prop             = zeros(size(limiter_vec));
    for jj = 1:length(this_data_plot)
        [~, limiter_vec_indx(jj)] = min(abs(this_data_plot(jj) - limiter_vec));
    end
    for jj = 1:length(limiter_vec)
        prop(jj) = sum(limiter_vec_indx == jj) ./ length(this_data_plot);
    end
    weight = prop(limiter_vec_indx);
    randx  = xpos + rand(1, length(this_data_plot)) .* weight - weight./2;
    scatter(randx, this_data_plot, 20, ...
            'MarkerFaceColor', dotcolor{ii}, ...
            'MarkerEdgeColor', 'none', ...
            'MarkerFaceAlpha', 0.1); hold on

    med_g = nanmedian(this_data_raw);
    q25   = quantile(this_data_raw, 0.25);
    q75   = quantile(this_data_raw, 0.75);

    med_g_plot = max(min(med_g, capHigh), capLow);
    q25_plot   = max(min(q25, capHigh), capLow);
    q75_plot   = max(min(q75, capHigh), capLow);

    plot([xpos-0.15 xpos+0.15],[med_g_plot med_g_plot], ...
         'color', dotcolor{ii}, 'linewidth',1.5); hold on
    plot([xpos-0.06 xpos+0.06],[q25_plot q25_plot], ...
         'color', dotcolor{ii}, 'linewidth',1.5); hold on
    plot([xpos-0.06 xpos+0.06],[q75_plot q75_plot], ...
         'color', dotcolor{ii}, 'linewidth',1.5); hold on
end

% 2) per-fish linking lines using RAW fish medians, capped only for display
for f = 1:nFish
    y_raw = y_fish_mat(f,:);
    valid = ~isnan(y_raw);
    if sum(valid) >= 2
        y_plot = max(min(y_raw, capHigh), capLow);
        plot(xpos3(valid), y_plot(valid), '-', ...
             'Color', [0.75 0.75 0.75], ...
             'LineWidth', 0.8); hold on
    end
end

% 3) pooled fish median dots, capped only for display
for g3 = 1:3
    this_y_raw  = y_fish_mat(:,g3);
    this_y_raw  = this_y_raw(~isnan(this_y_raw));
    this_y_plot = max(min(this_y_raw, capHigh), capLow);

    scatter(xpos3(g3)*ones(size(this_y_plot)), this_y_plot, 30, ...
            'MarkerFaceColor', pool_col{g3}, ...
            'MarkerEdgeColor', [0.3 0.3 0.3], ...
            'MarkerFaceAlpha', 0.45, ...
            'MarkerEdgeAlpha', 0.45, ...
            'LineWidth', 0.5); hold on
end

% --- axis appearance ---
yticklabels([]);
h = gca;
h.XAxis.Visible = 'off';
xticklabels({});
xticks([1 2.5 4]);
xlim([0 5]);
ylim([capLow capHigh]);
yticks(capLow:1:capHigh);

% horizontal zero reference line
yline(0, ':', 'Color', [0.25 0.25 0.25], 'LineWidth', 0.8); hold on

box off;

%% ---------- bout-level mixed model (fish random intercept) ----------
Tb = T_bout;
Tb = Tb(~isnan(Tb.y), :);

lme_bout = fitlme(Tb, 'y ~ 1 + group3 + (1|fishID)');
anova(lme_bout)
[beta_b, SE_b, stats_b] = fixedEffects(lme_bout,'DFMethod','Residual');
stats_b

H = [0 -1 1];
[p_32_b, F_32_b, df1_b, df2_b] = coefTest(lme_bout, H);
fprintf('Bout-level mixed model, av vs app: F=%.3f, p=%.4f\n', F_32_b, p_32_b);