clear all; close all; 
global gh
LoadFishNColorSel;
addpath('./stattool'); 
addpath('./disptool');

for ii = 1:5
    datavec{ii} = [];
end

% containers for mixed model + per-fish layer
all_y       = [];   % bout-wise metric: amplitude / peak count
all_group5  = [];   % 1..5
all_group3  = [];   % 1=blank, 2=app pooled, 3=av pooled
all_fishID  = [];

for fishsub = 1:length(totalfishsub)
    gh.param.fishid = totalfishsub(fishsub);
    fm_behavim_main
    for ii = 1:sessionn(fishsub)
        if ~ismember(ii, gh.param.ExcludedSession)

            stim_boutn = find(gh.data.boutmtx(:,1) == ii & ...
                              gh.data.boutmtx(:,2) > 20 & ...
                              gh.data.boutmtx(:,4) <= 30);
            prestim_boutn = find(gh.data.boutmtx(:,1) == ii & ...
                                 gh.data.boutmtx(:,2) > 10 & ...
                                 gh.data.boutmtx(:,4) <= 20);

            trcode = gh.param.fishlog.trialdetails.trial(ii,1);

            % 5-group labelling
            if     trcode == 8
                grp5 = 1;
            elseif trcode == 1
                grp5 = 2;
            elseif trcode == 2
                grp5 = 3;
            elseif trcode == 4
                grp5 = 4;
            elseif trcode == 5
                grp5 = 5;
            else
                continue;
            end

            % 3-group pooled labels
            if     trcode == 8
                grp3 = 1;               % blank
            elseif trcode == 1 || trcode == 2
                grp3 = 2;               % appetitive
            elseif trcode == 4 || trcode == 5
                grp3 = 3;               % aversive
            else
                continue;
            end

            sumval    = 0;
            validbout = 0;
            for jj = 1:length(prestim_boutn)


                X      = gh.data.anglevect(gh.data.bout_details(prestim_boutn(jj),3): ...
                                           gh.data.bout_details(prestim_boutn(jj),5));
                X0     = X - max(X);
                X0_inv = -X0;
                pks_X0     = findpeaks(X0);
                pks_X0_inv = findpeaks(X0_inv);
                pkcount    = (length(pks_X0) + length(pks_X0_inv)) / 2;
                if pkcount >= 4
                    validbout = validbout + 1;
                    val = gh.data.bout_details(prestim_boutn(jj),6) / pkcount;
                    sumval = sumval + val;
                end
            end

            if validbout > 0
                avgval = sumval ./ validbout - boutundint_overall(fishsub);
            else
                avgval = 0;
            end

            for jj = 1:length(stim_boutn)
                X      = gh.data.anglevect(gh.data.bout_details(stim_boutn(jj),3): ...
                                           gh.data.bout_details(stim_boutn(jj),5));
                X0     = X - max(X);
                X0_inv = -X0;
                pks_X0     = findpeaks(X0);
                pks_X0_inv = findpeaks(X0_inv);
                pkcount    = (length(pks_X0) + length(pks_X0_inv)) / 2;
                if pkcount >= 4 
                    val = (gh.data.bout_details(stim_boutn(jj),6) / pkcount) - avgval - boutundint_overall(fishsub);

                    % original per-bout plotting data (5 chems)
                    datavec{grp5} = [datavec{grp5}; val];

                    % data for mixed model and per-fish summaries
                    all_y       = [all_y;       val];
                    all_group5  = [all_group5;  grp5];
                    all_group3  = [all_group3;  grp3];
                    all_fishID  = [all_fishID;  gh.param.fishid];
                end
            end
        end
    end
end

%% ---------- bout-level table ----------
T_bout = table(all_y, categorical(all_group3), ...
               categorical(all_group5), categorical(all_fishID), ...
               'VariableNames', {'y','group3','group5','fishID'});
T_bout = T_bout(~isnan(T_bout.y), :);

%% ---------- bout-median per fish × 3-group ----------
[G3_f, fish_u3, group3_u] = findgroups(T_bout.fishID, T_bout.group3);
y_fish3 = splitapply(@median, T_bout.y, G3_f);

T_fish3 = table(y_fish3, categorical(group3_u), categorical(fish_u3), ...
                'VariableNames', {'y','group','fishID'});

fish_ids   = categories(T_fish3.fishID);
nFish      = numel(fish_ids);
y_fish_mat = nan(nFish, 3);  % col1=blank, col2=app, col3=av

for f = 1:nFish
    thisFish = fish_ids{f};
    idx_f    = T_fish3.fishID == thisFish;
    g_f      = double(T_fish3.group(idx_f));  % 1,2,3
    y_f      = T_fish3.y(idx_f);              % medians per fish×group
    for k = 1:numel(g_f)
        y_fish_mat(f, g_f(k)) = y_f(k);
    end
end

%% ---------- bout-level distributions per chemical (RAW source) ----------
datavec_bout5 = cell(1,5);
for k = 1:5
    datavec_bout5{k} = datavec{k};
end

%% ---------- visual cap settings ----------
capLow  = -0.004;
capHigh =  0.004;

%% ---------- figure: bout-level violin + medians/IQR + per-fish dots/lines ----------
figure(1); clf; hold on;
set(gcf,'Position',[800 800 300 400])

% x positions for 5 chemicals
xpos5 = nan(1,5);
xpos5(1) = 1;
xpos5(2) = 2.2;
xpos5(3) = 2.8;
xpos5(4) = 3.7;
xpos5(5) = 4.3;

% x positions for pooled 3 groups (blank/app/av)
xpos3 = [1 2.5 4];

% colors
dotcolor{1} = blankColor;
dotcolor{2} = appeColor; dotcolor{3} = appeColor;
dotcolor{4} = averColor; dotcolor{5} = averColor;

pool_col{1} = blankColor;
pool_col{2} = appeColor;
pool_col{3} = averColor;

% 1) 5-group violin-style bout spread + bout medians/IQR
for ii = 1:5
    xpos = xpos5(ii);
    this_data_raw = datavec_bout5{ii};
    if isempty(this_data_raw)
        continue
    end

    % visual capping only
    this_data = max(min(this_data_raw, capHigh), capLow);

    % jitter width based on local density using capped data
    dmin = min(this_data);
    dmax = max(this_data);
    if dmin == dmax
        randx = xpos + (rand(size(this_data))-0.5)*0.05;
    else
        limiter_vec = linspace(dmin, dmax, 60);
        limiter_vec_indx = zeros(size(this_data));
        prop = zeros(size(limiter_vec));
        for jj = 1:length(this_data)
            [~, limiter_vec_indx(jj)] = min(abs(this_data(jj) - limiter_vec));
        end
        for jj = 1:length(limiter_vec)
            prop(jj) = sum(limiter_vec_indx == jj) ./ length(this_data);
        end
        weight = prop(limiter_vec_indx) * 1.5;
        randx  = xpos + rand(1, length(this_data)) .* weight - weight./2;
    end

    scatter(randx, this_data, 20, ...
            'MarkerFaceColor', dotcolor{ii}, ...
            'MarkerEdgeColor', 'none', ...
            'MarkerFaceAlpha', .12); hold on

mid_raw = nanmedian(this_data_raw);
u_raw   = quantile(this_data_raw, 0.75);
l_raw   = quantile(this_data_raw, 0.25);

% Cap only the displayed locations of the summary bars
mid_plot = max(min(mid_raw, capHigh), capLow);
u_plot   = max(min(u_raw, capHigh), capLow);
l_plot   = max(min(l_raw, capHigh), capLow);

plot([xpos-0.15 xpos+0.15], [mid_plot mid_plot], ...
     'color', dotcolor{ii}, 'linewidth', 1.5); hold on
plot([xpos-0.06 xpos+0.06], [u_plot u_plot], ...
     'color', dotcolor{ii}, 'linewidth', 1.5); hold on
plot([xpos-0.06 xpos+0.06], [l_plot l_plot], ...
     'color', dotcolor{ii}, 'linewidth', 1.5); hold on

    violin(xpos, this_data, ...
           'facecolor', dotcolor{ii}, ...
           'scaling', 0.0025, ...
           'facealpha', 0.075, ...
           'style', 2); hold on
end

% 2) per-fish linking lines using fish MEDIANS
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

% 3) overlay pooled fish MEDIAN dots
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

%% ---------- axis cosmetics with visual cap ----------
ylim([capLow capHigh]);
yticks(linspace(capLow, capHigh, 9));
xticks([1 2.5 4]);
xlim([0 5]);

% dotted zero line
plot([0 5], [0 0], ':', 'Color', [0.25 0.25 0.25], 'LineWidth', 0.8);

box off;
h = gca;
h.XAxis.Visible = 'off';
yticklabels([]);
xticklabels({});

%% ---------- bout-level mixed models (fish random, RAW values) ----------
T3 = T_bout(:,{'y','group3','fishID'});
lme3 = fitlme(T3, 'y ~ 1 + group3 + (1|fishID)');
disp('=== 3-group pooled model ===');
anova(lme3)
[beta3, SE3, stats3] = fixedEffects(lme3,'DFMethod','Residual');
stats3

% av vs app (group3_3 - group3_2)
H = [0 -1 1];
[p_32, F_32, df1_32, df2_32] = coefTest(lme3, H);
fprintf('Bout-level mixed model, av vs app (this metric): F=%.3f, p=%.4f\n', F_32, p_32);