clear; close all; global gh

LoadFishNColorSel;
addpath('./stattool');
addpath('./disptool');

% =========================================================
% SETTINGS
% =========================================================
T_ONLY      = 1; ST_COMBINED = 2;
APPETITIVE  = 1; AVERSIVE    = 2;
nSubgroups  = 2; nStimGroups = 2;

metricName = 'boutbias';

% Plot appearance
mkrsz           = 20;       dotAlpha        = 0.12;
violinAlpha     = 0.075;    medianLineWidth = 1.5;
iqrLineWidth    = 1.2;      

violinScaling   = 0.50;     displayMin      = 0;
displayMax      = 10;       intVal          = 2;
weightVal       = 0.4;      yticksVal       = 0:1:10;

medianHalfWidth = 0.08;     iqrHalfWidth    = 0.03;

specifiedColor = cell(1, nStimGroups);
specifiedColor{APPETITIVE} = appeColor;
specifiedColor{AVERSIVE}   = averColor;

fishLineColor = [0.68 0.68 0.68];
fishDotFace   = [0.78 0.78 0.78];
fishDotEdge   = [0.30 0.30 0.30];
fishLineWidth = 0.45;
fishDotSize   = 34;
fishDotAlpha  = 0.55;

% Inverted subgroup positions: S-T left, T right.
subgroupXOffset = nan(1, nSubgroups);
subgroupXOffset(T_ONLY)      =  0.125;
subgroupXOffset(ST_COMBINED) = -0.125;

% =========================================================
% PREALLOCATE
% =========================================================
% Each cell stores a column vector of bout-bias values.
DataGroup = cell(nSubgroups, nStimGroups);

% One row per bout for the mixed-effects model.
all_y        = []; all_stim     = [];
all_subgroup = []; all_fishID   = [];

% =========================================================
% COLLECT BOUT-LEVEL DATA
% =========================================================
for fishsub = 1:numel(totalfishsub)

    gh.param.fishid = totalfishsub(fishsub);
    fm_behavim_main

    for sess = 1:sessionn(fishsub)

        if ismember(sess, gh.param.ExcludedSession)
            continue
        end

        trialtype = gh.param.fishlog.trialdetails.trial(sess, 1);

        if ismember(trialtype, [1 2])
            stimGroup = APPETITIVE;
        elseif ismember(trialtype, [4 5])
            stimGroup = AVERSIVE;
        else
            continue
        end

        stim_boutn = find( ...
            gh.data.boutmtx(:,1) == sess & ...
            gh.data.boutmtx(:,2) > 20 & ...
            gh.data.boutmtx(:,4) <= 30);

        for boutn = 1:numel(stim_boutn)
            thisBout = stim_boutn(boutn);
            boutbias = abs( ...
                -gh.data.bout_details(thisBout, 7) - ...
                anglebias_overall(fishsub));
            isChemosaccadeCoupled = any( ...
                gh.data.simuMtx(:,9) == gh.data.boutmtx(thisBout,2) & ...
                gh.data.simuMtx(:,1) == gh.data.boutmtx(thisBout,1));
            if isChemosaccadeCoupled
                subgroup = ST_COMBINED;
            else
                subgroup = T_ONLY;
            end
            % For group-wise plotting.
            DataGroup{subgroup, stimGroup} = ...
                [DataGroup{subgroup, stimGroup}; boutbias];
            % For the bout-level mixed-effects model.
            all_y(end+1,1)        = boutbias;
            all_stim(end+1,1)     = stimGroup;
            all_subgroup(end+1,1) = subgroup;
            all_fishID(end+1,1)   = gh.param.fishid;
        end
    end
end

% =========================================================
% BUILD ANALYSIS TABLE
% =========================================================
tbl = table( ...
    all_y, ...
    categorical(all_stim, [APPETITIVE AVERSIVE], ...
    {'Appetitive','Aversive'}), ...
    categorical(all_subgroup, [T_ONLY ST_COMBINED], ...
    {'T','ST'}), ...
    categorical(all_fishID), ...
    'VariableNames', {'y','StimGroup','Subgroup','FishID'});

tbl = tbl(~isnan(tbl.y), :);

% Define reference categories explicitly:
% Intercept = Appetitive / T-only.
tbl.StimGroup = reordercats(tbl.StimGroup, {'Appetitive','Aversive'});
tbl.Subgroup  = reordercats(tbl.Subgroup, {'T','ST'});

% =========================================================
% PLOTTING
% =========================================================
rng(1)

figure(1); clf; hold on
set(gcf, 'Position', [1 650 420 620])

fish_ids = categories(tbl.FishID);
nFish    = numel(fish_ids);

% Each entry is the median bout-bias for one fish in one
% stimulus-group × subgroup combination.
y_fish_mat = nan(nFish, nStimGroups, nSubgroups);

for f = 1:nFish
    thisFish = fish_ids{f};
    for stimGroup = 1:nStimGroups

        if stimGroup == APPETITIVE
            stimLabel = 'Appetitive';
        else
            stimLabel = 'Aversive';
        end

        for subgroup = 1:nSubgroups

            if subgroup == T_ONLY
                subgroupLabel = 'T';
            else
                subgroupLabel = 'ST';
            end
            idx = tbl.FishID == thisFish & ...
                tbl.StimGroup == stimLabel & ...
                tbl.Subgroup == subgroupLabel;
            if any(idx)
                y_fish_mat(f, stimGroup, subgroup) = median(tbl.y(idx));
            end
        end
    end
end

% Bout-level violin distributions, dots, medians, and IQRs.
for stimGroup = 1:nStimGroups
    for subgroup = 1:nSubgroups

        datavec = DataGroup{subgroup, stimGroup};
        datavec = datavec(~isnan(datavec));

        if isempty(datavec)
            continue
        end

        xpos = stimGroup + subgroupXOffset(subgroup);

        % Clipping is for display only.
        % The mixed model below uses the original, unclipped values.
        tempdata = max(min(datavec, displayMax), displayMin);

        limiter_vec      = displayMin:intVal:displayMax;
        limiter_vec_indx = zeros(numel(tempdata), 1);

        for mm = 1:numel(tempdata)
            [~, limiter_vec_indx(mm)] = ...
                min(abs(tempdata(mm) - limiter_vec));
        end

        prop = accumarray( ...
            limiter_vec_indx, ...
            1, ...
            [numel(limiter_vec), 1]) ./ numel(tempdata);

        weight = prop(limiter_vec_indx) * weightVal;
        randx  = xpos + rand(size(tempdata)) .* weight - weight ./ 2;

        violin(xpos, tempdata, ...
            'facecolor', specifiedColor{stimGroup}, ...
            'scaling', violinScaling, ...
            'facealpha', violinAlpha, ...
            'style', 2);

        scatter(randx, tempdata, mkrsz, ...
            'MarkerFaceColor', specifiedColor{stimGroup}, ...
            'MarkerEdgeColor', 'none', ...
            'MarkerFaceAlpha', dotAlpha);

        mid_disp = median(tempdata);
        u_disp   = quantile(tempdata, 0.75);
        l_disp   = quantile(tempdata, 0.25);

        plot([xpos-medianHalfWidth, xpos+medianHalfWidth], ...
            [mid_disp, mid_disp], ...
            'Color', specifiedColor{stimGroup}, ...
            'LineWidth', medianLineWidth);

        plot([xpos-iqrHalfWidth, xpos+iqrHalfWidth], ...
            [u_disp, u_disp], ...
            'Color', specifiedColor{stimGroup}, ...
            'LineWidth', iqrLineWidth);

        plot([xpos-iqrHalfWidth, xpos+iqrHalfWidth], ...
            [l_disp, l_disp], ...
            'Color', specifiedColor{stimGroup}, ...
            'LineWidth', iqrLineWidth);

    end
end

% Join fish-level medians between T and S-T in each stimulus category.
for stimGroup = 1:nStimGroups

    xpos_T  = stimGroup + subgroupXOffset(T_ONLY);
    xpos_ST = stimGroup + subgroupXOffset(ST_COMBINED);

    for f = 1:nFish
        yT  = y_fish_mat(f, stimGroup, T_ONLY);
        yST = y_fish_mat(f, stimGroup, ST_COMBINED);
        if ~isnan(yT) && ~isnan(yST)
            ypair_plot = max(min([yT yST], displayMax), displayMin);
            plot([xpos_T xpos_ST], ypair_plot, '-', ...
                'Color', [fishLineColor 0.50], ...
                'LineWidth', fishLineWidth);
        end
    end
end

% Fish-level median dots.
for stimGroup = 1:nStimGroups
    for subgroup = 1:nSubgroups
        xpos = stimGroup + subgroupXOffset(subgroup);
        this_y = squeeze(y_fish_mat(:, stimGroup, subgroup));
        this_y = this_y(~isnan(this_y));
        if isempty(this_y)
            continue
        end

        this_y_plot = max(min(this_y, displayMax), displayMin);

        scatter(xpos * ones(size(this_y_plot)), this_y_plot, fishDotSize, ...
            'MarkerFaceColor', fishDotFace, ...
            'MarkerEdgeColor', fishDotEdge, ...
            'MarkerFaceAlpha', fishDotAlpha, ...
            'MarkerEdgeAlpha', fishDotAlpha, ...
            'LineWidth', 0.5);
    end
end

xticks([0.875 1.125 1.875 2.125])
xticklabels([])
yticks(yticksVal)
yticklabels([])
xlim([0.5 2.5])
ylim([displayMin displayMax])
pbaspect([1 1.8 1])
box off

% =========================================================
% MIXED-EFFECTS MODEL
% =========================================================
disp(' ')
disp('========== Mixed-effects model: bout bias ==========')

% Fixed effects:
%   StimGroup, Subgroup, and StimGroup × Subgroup interaction
%
% Random effect:
%   One random intercept per fish

lme = fitlme(tbl, ...
    'y ~ StimGroup*Subgroup + (1|FishID)');

disp(anova(lme))

[beta, SE, stats] = fixedEffects(lme, 'DFMethod', 'Residual');

disp(table(beta, SE, stats.tStat, stats.pValue, ...
    'VariableNames', {'Estimate','SE','tStat','pValue'}, ...
    'RowNames', lme.CoefficientNames))

% =========================================================
% PLANNED MIXED-MODEL CONTRASTS
% =========================================================
names = lme.CoefficientNames;

% Appetitive: S-T minus T.
H_app = zeros(1, numel(names));
H_app(strcmp(names, 'Subgroup_ST')) = 1;

[p_app, F_app, df1_app, df2_app] = coefTest(lme, H_app);

% Aversive: S-T minus T.
H_av = zeros(1, numel(names));
H_av(strcmp(names, 'Subgroup_ST')) = 1;
H_av(strcmp(names, 'StimGroup_Aversive:Subgroup_ST')) = 1;

[p_av, F_av, df1_av, df2_av] = coefTest(lme, H_av);

% Does the S-T-versus-T effect differ by stimulus category?
H_int = zeros(1, numel(names));
H_int(strcmp(names, 'StimGroup_Aversive:Subgroup_ST')) = 1;

[p_int, F_int, df1_int, df2_int] = coefTest(lme, H_int);

fprintf('\nS-T versus T contrasts from the mixed-effects model:\n');

fprintf('  Appetitive: F(%g,%g) = %.3f, p = %.6g\n', ...
    df1_app, df2_app, F_app, p_app);

fprintf('  Aversive:   F(%g,%g) = %.3f, p = %.6g\n', ...
    df1_av, df2_av, F_av, p_av);

fprintf('  Interaction: F(%g,%g) = %.3f, p = %.6g\n', ...
    df1_int, df2_int, F_int, p_int);