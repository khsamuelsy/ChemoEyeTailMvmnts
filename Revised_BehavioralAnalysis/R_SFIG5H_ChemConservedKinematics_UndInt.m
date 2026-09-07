clear; close all;global gh
LoadFishNColorSel;addpath('./stattool');addpath('./disptool');
% =========================================================
% SETTINGS
% =========================================================
T_ONLY      = 1;    ST_COMBINED = 2;
APPETITIVE  = 1;    AVERSIVE    = 2;
nSubgroups  = 2;    nStimGroups = 2;
metricName = 'Average undulation interval';

% =========================================================
% DISPLAY SETTINGS
% =========================================================
mkrsz           = 20;   dotAlpha        = 0.12; violinAlpha     = 0.075;
medianLineWidth = 1.5;  iqrLineWidth    = 1.2;  
violinScaling   = 0.0008;

% Visualization limits only.
% No raw observations are capped before calculating statistics
% or fitting the mixed-effects model.
displayMin      = 0.032;
displayMax      = 0.048;

intVal          = 0.002;weightVal       = 0.5;
yticksVal       = 0.032:0.002:0.048;

medianHalfWidth = 0.08;iqrHalfWidth    = 0.03;

specifiedColor = cell(1, nStimGroups);
specifiedColor{APPETITIVE} = appeColor;
specifiedColor{AVERSIVE}   = averColor;

fishLineColor = [0.68 0.68 0.68];
fishDotFace   = [0.78 0.78 0.78];
fishDotEdge   = [0.30 0.30 0.30];
fishLineWidth = 0.45;
fishDotSize   = 34;
fishDotAlpha  = 0.55;

% Subgroup locations: S-T on left, T-only on right.
subgroupXOffset = nan(1, nSubgroups);
subgroupXOffset(T_ONLY)      =  0.125;
subgroupXOffset(ST_COMBINED) = -0.125;

% =========================================================
% PREALLOCATE
% =========================================================
% Each cell holds raw, uncapped interval values:
% DataGroup{subgroup, stimGroup}.
DataGroup = cell(nSubgroups, nStimGroups);

% One raw row per included bout for LME fitting.
all_y        = [];  all_stim     = [];
all_subgroup = [];  all_fishID   = [];

% =========================================================
% COLLECT BOUT-LEVEL UNDULATION INTERVALS
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

            % -------------------------------------------------
            % Calculate raw average undulation interval
            % -------------------------------------------------
            frameStart = gh.data.bout_details(thisBout, 3);
            frameEnd   = gh.data.bout_details(thisBout, 5);

            X = gh.data.anglevect(frameStart:frameEnd);

            % Estimate the number of undulation cycles from the
            % detected extrema in the trace and in its inverse.
            X0     = X - max(X);
            X0_inv = -X0;

            pks_X0     = findpeaks(X0);
            pks_X0_inv = findpeaks(X0_inv);

            pkcount = (numel(pks_X0) + numel(pks_X0_inv)) / 2;

            % Retain only bouts containing at least four cycles.
            if ~isfinite(pkcount) || pkcount < 4
                continue
            end

            dur = gh.data.bout_details(thisBout, 6);
            avgUndulationInterval = dur / pkcount;

            % -------------------------------------------------
            % Determine subgroup: T-only or S-T
            % -------------------------------------------------
            isChemosaccadeCoupled = any( ...
                gh.data.simuMtx(:,9) == gh.data.boutmtx(thisBout, 2) & ...
                gh.data.simuMtx(:,1) == gh.data.boutmtx(thisBout, 1));

            if isChemosaccadeCoupled
                subgroup = ST_COMBINED;
            else
                subgroup = T_ONLY;
            end

            % -------------------------------------------------
            % Store raw, uncapped values
            % -------------------------------------------------
            DataGroup{subgroup, stimGroup} = ...
                [DataGroup{subgroup, stimGroup}; avgUndulationInterval];

            all_y(end+1,1)        = avgUndulationInterval;
            all_stim(end+1,1)     = stimGroup;
            all_subgroup(end+1,1) = subgroup;
            all_fishID(end+1,1)   = gh.param.fishid;
        end
    end
end

% =========================================================
% BUILD TABLE FOR THE MIXED-EFFECTS MODEL
% =========================================================
tbl = table( ...
    all_y, ...
    categorical(all_stim, [APPETITIVE AVERSIVE], ...
        {'Appetitive', 'Aversive'}), ...
    categorical(all_subgroup, [T_ONLY ST_COMBINED], ...
        {'T', 'ST'}), ...
    categorical(all_fishID), ...
    'VariableNames', {'y', 'StimGroup', 'Subgroup', 'FishID'});

% Retain valid, raw responses.
tbl = tbl(isfinite(tbl.y) & tbl.y > 0, :);

% Set the reference condition:
% Intercept = Appetitive T-only.
tbl.StimGroup = reordercats( ...
    tbl.StimGroup, {'Appetitive', 'Aversive'});

tbl.Subgroup = reordercats( ...
    tbl.Subgroup, {'T', 'ST'});

% =========================================================
% RAW PER-FISH MEDIANS FOR PLOTTING
% =========================================================
fishIDs = categories(tbl.FishID);
nFish   = numel(fishIDs);

% Dimensions: FishID × StimGroup × Subgroup.
% Each value is a median calculated from raw, uncapped data.
yFishMedian = nan(nFish, nStimGroups, nSubgroups);

for f = 1:nFish

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

            idx = tbl.FishID == fishIDs{f} & ...
                  tbl.StimGroup == stimLabel & ...
                  tbl.Subgroup == subgroupLabel;

            if any(idx)

                % Raw values are used for the per-fish median.
                yFishMedian(f, stimGroup, subgroup) = ...
                    median(tbl.y(idx), 'omitnan');
            end
        end
    end
end

% =========================================================
% PLOTTING
% =========================================================
rng(1)

figure(1); clf; hold on
set(gcf, 'Position', [1 650 420 620])

% ---------------------------------------------------------
% Individual bouts, violin, raw group summaries
% ---------------------------------------------------------
for stimGroup = 1:nStimGroups
    for subgroup = 1:nSubgroups

        % Raw observations for this condition.
        datavec = DataGroup{subgroup, stimGroup};
        datavec = datavec(isfinite(datavec) & datavec > 0);

        if isempty(datavec)
            continue
        end

        xpos = stimGroup + subgroupXOffset(subgroup);

        % A capped copy is used only for the rendered violin and dots.
        yDisplay = max(min(datavec, displayMax), displayMin);

        % Density-dependent horizontal jitter using displayed y positions.
        limiterVec = displayMin:intVal:displayMax;

        [~, binIdx] = min(abs(yDisplay - limiterVec), [], 2);

        binCount = accumarray( ...
            binIdx, ...
            1, ...
            [numel(limiterVec), 1]);

        binProp = binCount / numel(yDisplay);
        localWidth = binProp(binIdx) * weightVal;

        xJitter = xpos + rand(size(yDisplay)) .* localWidth - ...
                  localWidth ./ 2;

        violin(xpos, yDisplay, ...
            'facecolor', specifiedColor{stimGroup}, ...
            'scaling', violinScaling, ...
            'facealpha', violinAlpha, ...
            'style', 2);

        scatter(xJitter, yDisplay, mkrsz, ...
            'MarkerFaceColor', specifiedColor{stimGroup}, ...
            'MarkerEdgeColor', 'none', ...
            'MarkerFaceAlpha', dotAlpha);

        % Calculate group summaries from raw, uncapped observations.
        % qRaw = [Q1, median, Q3].
        qRaw = quantile(datavec, [0.25 0.50 0.75]);

        % Cap only their visible vertical locations.
        qDisplay = max(min(qRaw, displayMax), displayMin);

        % Median across all raw bouts in this condition.
        plot([xpos-medianHalfWidth, xpos+medianHalfWidth], ...
            [qDisplay(2), qDisplay(2)], ...
            'Color', specifiedColor{stimGroup}, ...
            'LineWidth', medianLineWidth);

        % Q1 across all raw bouts in this condition.
        plot([xpos-iqrHalfWidth, xpos+iqrHalfWidth], ...
            [qDisplay(1), qDisplay(1)], ...
            'Color', specifiedColor{stimGroup}, ...
            'LineWidth', iqrLineWidth);

        % Q3 across all raw bouts in this condition.
        plot([xpos-iqrHalfWidth, xpos+iqrHalfWidth], ...
            [qDisplay(3), qDisplay(3)], ...
            'Color', specifiedColor{stimGroup}, ...
            'LineWidth', iqrLineWidth);

    end
end

% ---------------------------------------------------------
% Paired fish lines: raw per-fish medians
% ---------------------------------------------------------
for stimGroup = 1:nStimGroups

    xpos_T  = stimGroup + subgroupXOffset(T_ONLY);
    xpos_ST = stimGroup + subgroupXOffset(ST_COMBINED);

    for f = 1:nFish

        % Raw per-fish medians.
        yT_raw  = yFishMedian(f, stimGroup, T_ONLY);
        yST_raw = yFishMedian(f, stimGroup, ST_COMBINED);
        if ~isnan(yT_raw) && ~isnan(yST_raw)
            % Cap only the plotted endpoints.
            yPairDisplay = max( ...
                min([yT_raw, yST_raw], displayMax), ...
                displayMin);

            plot([xpos_T, xpos_ST], yPairDisplay, '-', ...
                'Color', [fishLineColor, 0.50], ...
                'LineWidth', fishLineWidth);
        end
    end
end

% ---------------------------------------------------------
% Fish-level median dots: raw medians, display-capped only
% ---------------------------------------------------------
for stimGroup = 1:nStimGroups
    for subgroup = 1:nSubgroups

        xpos = stimGroup + subgroupXOffset(subgroup);

        % Raw per-fish median values.
        yFishRaw = squeeze(yFishMedian(:, stimGroup, subgroup));
        yFishRaw = yFishRaw(~isnan(yFishRaw));

        if isempty(yFishRaw)
            continue
        end

        % Only their visual y-positions are capped.
        yFishDisplay = max(min(yFishRaw, displayMax), displayMin);

        scatter(xpos * ones(size(yFishDisplay)), ...
            yFishDisplay, fishDotSize, ...
            'MarkerFaceColor', fishDotFace, ...
            'MarkerEdgeColor', fishDotEdge, ...
            'MarkerFaceAlpha', fishDotAlpha, ...
            'MarkerEdgeAlpha', fishDotAlpha, ...
            'LineWidth', 0.5);
    end
end

xticks([0.875 1.125 1.875 2.125])
xticklabels([])
yticks(yticksVal);yticklabels([])
xlim([0.5 2.5]);ylim([displayMin displayMax])
pbaspect([1 1.8 1])
box off

% =========================================================
% MIXED-EFFECTS MODEL
% =========================================================
disp(' ')
fprintf('========== Mixed-effects model: %s ==========\n', metricName)

% Fixed effects:
% - StimGroup
% - Subgroup
% - StimGroup × Subgroup
%
% Random effects:
% - Fish-specific random intercept
% The model uses raw, uncapped interval values from tbl.y.
lme = fitlme(tbl, ...
    'y ~ StimGroup*Subgroup + (1|FishID)');

disp(anova(lme))

[beta, SE, stats] = fixedEffects(lme, 'DFMethod', 'Residual');

disp(table(beta, SE, stats.tStat, stats.pValue, ...
    'VariableNames', {'Estimate', 'SE', 'tStat', 'pValue'}, ...
    'RowNames', lme.CoefficientNames))

% =========================================================
% PLANNED MIXED-MODEL CONTRASTS
% =========================================================
names = lme.CoefficientNames;

% S-T minus T-only within appetitive trials.
H_app = zeros(1, numel(names));
H_app(strcmp(names, 'Subgroup_ST')) = 1;

[p_app, F_app, df1_app, df2_app] = coefTest(lme, H_app);

% S-T minus T-only within aversive trials.
H_av = zeros(1, numel(names));
H_av(strcmp(names, 'Subgroup_ST')) = 1;
H_av(strcmp(names, ...
    'StimGroup_Aversive:Subgroup_ST')) = 1;

[p_av, F_av, df1_av, df2_av] = coefTest(lme, H_av);

fprintf('\nS-T versus T-only contrasts from mixed-effects model:\n');
fprintf('  Appetitive: F(%g,%g) = %.3f, p = %.6g\n', ...
    df1_app, df2_app, F_app, p_app);
fprintf('  Aversive:   F(%g,%g) = %.3f, p = %.6g\n', ...
    df1_av, df2_av, F_av, p_av);
