%% R_SFIG5IJ_ChemConservedSTCoordination_Shuffled.m
%
% Within-trial saccade-time circular shuffle and nested-model surrogate test.
%
% Real and shuffled data are analyzed using the same event-level models:
%
% Direction:
%   fitglme(Y ~ 1 + (1|FishID), Distribution = Binomial, Link = logit)
%
% Delay:
%   fitlme(Delay ~ 1 + (1|FishID))
%
% The empirical shuffle test compares the observed real fixed intercept
% with the distribution of fixed intercepts calculated from shuffled,
% reidentified S-T events.
%
% Saccade matrix:
%   col 1 = session
%   col 2 = onset time
%   col 3 = onset frame
%   col 4 = offset time
%   col 5 = offset frame
%   col 7 = direction
%
% Reidentified shuffled simuMtx:
%   cols 1:7 = shuffled saccademtx row
%   cols 8:15 = closest matching bout_details row
%   col 9 = tail onset
%   col 11 = tail offset
%   col 14 = tail-bias measure / bout_details(:,7)
%
% -------------------------------------------------------------------------
clear;
close all;
clc;

global gh
LoadFishNColorSel;
rng(1);

addpath('./stattool');
addpath('./disptool');

%% Settings
stimcase = [4,5];

stimulusStart = 20;
stimulusEnd = 30;

STwindow = 0.5;
threshold_boutbias = 1.5;

% Saccade/eye video sampling rate.
% The temporal grid is 1/8 s = 0.125 s per saccade-camera frame.
saccadeFrameRate = 8;

% Use 100 for debugging. Use 1000–5000 for final analysis.
nShuffle = 400;
plotcolor = averColor;
conditionName = 'Aversive';
nFish = length(totalfishsub);

%% ------------------------------------------------------------------------
%% 1. BUILD REAL EVENT-LEVEL TABLES
%% ------------------------------------------------------------------------

realY = [];realYFishID = [];
realDelay =[];realDelayFishID = [];


for fishsub = 1:nFish

    gh.param.fishid = totalfishsub(fishsub);
    fm_behavim_main;

    for ii = 1:sessionn(fishsub)

        if ismember(ii,gh.param.ExcludedSession)
            continue
        end

        if ~ismember(gh.param.fishlog.trialdetails.trial(ii,1),stimcase)
            continue
        end

        totalRealSimu = find( ...
            gh.data.simuMtx(:,1) == ii & ...
            gh.data.simuMtx(:,2) > stimulusStart & ...
            gh.data.simuMtx(:,4) <= stimulusEnd);

        for simun = 1:numel(totalRealSimu)

            rr = totalRealSimu(simun);

            boutdir = sign( ...
                gh.data.simuMtx(rr,14) * -1 - ...
                anglebias_overall(fishsub));

            boutbias = abs( ...
                gh.data.simuMtx(rr,14) * -1 - ...
                anglebias_overall(fishsub));

            saccadedir = gh.data.simuMtx(rr,7);

            % Same delay convention as original:
            % delay  = T_on - S_on
            delayvec = [gh.data.simuMtx(rr,9) - gh.data.simuMtx(rr,2)];

            realDelay(end+1,1) = delayvec; 
            realDelayFishID(end+1,1) = totalfishsub(fishsub);

            % Direction analysis follows the original threshold:
            % only tail flips with sufficiently large turning bias.
            if boutbias >= threshold_boutbias && ...
                    ismember(boutdir,[-1,1]) && ...
                    ismember(saccadedir,[-1,1])
                isIpsi = double(saccadedir == boutdir);
                realY(end+1,1) = isIpsi;
                realYFishID(end+1,1) = totalfishsub(fishsub); 
            end
            clear delayvec
        end
    end
end

%% ------------------------------------------------------------------------
%% 2. FIT ORIGINAL EVENT-LEVEL MODELS TO REAL DATA
%% ------------------------------------------------------------------------

% Direction: ipsilateral probability compared with 0.5.
tblRealDirection = table( ...
    realY, ...
    categorical(realYFishID), ...
    'VariableNames',{'Y','FishID'});

glmeRealDirection = fitglme(tblRealDirection, ...
    'Y ~ 1 + (1|FishID)', ...
    'Distribution','Binomial', ...
    'Link','logit');

betaRealDirection = fixedEffects(glmeRealDirection);

% Delay 1: T_on - S_on.
tblRealDelay = table( ...
    realDelay, ...
    categorical(realDelayFishID), ...
    'VariableNames',{'Delay','FishID'});

lmeRealDelay = fitlme(tblRealDelay, ...
    'Delay ~ 1 + (1|FishID)');
betaRealDelay = fixedEffects(lmeRealDelay);

fprintf('\n====================================================\n');
fprintf('Real nested-model results: %s\n',conditionName);
fprintf('====================================================\n');

fprintf('Direction GLME intercept [logit ipsilateral probability]: %.6f\n', ...
    betaRealDirection(1));

fprintf('Direction GLME implied ipsilateral probability: %.4f\n', ...
    1 ./ (1 + exp(-betaRealDirection(1))));

fprintf('Delay LME intercept, T_on - S_on: %.6f s\n', ...
    betaRealDelay(1));

%% ------------------------------------------------------------------------
%% 3. SHUFFLE SACCADE TIME AND FIT IDENTICAL MODELS
%% ------------------------------------------------------------------------

shuffleDirectionBeta = nan(nShuffle,1);
shuffleDelayBeta = nan(nShuffle,1);

for shufflen = 1:nShuffle

    shuffleY = [];
    shuffleYFishID = [];

    shuffleDelay = [];
    shuffleDelayFishID = [];

    for fishsub = 1:nFish

        gh.param.fishid = totalfishsub(fishsub);
        fm_behavim_main;

        saccademtxOriginal = gh.data.saccademtx;
        boutDetails = gh.data.bout_details;

        % Reset from the real data on every iteration.
        saccademtxShuffled = saccademtxOriginal;

        %% 3A. Shift saccades within every trial using 8-fps steps

        for ii = 1:sessionn(fishsub)

            if ismember(ii,gh.param.ExcludedSession)
                continue
            end

            if ~ismember( ...
                    gh.param.fishlog.trialdetails.trial(ii,1),stimcase)
                continue
            end

            saccadeRows = find( ...
                saccademtxOriginal(:,1) == ii & ...
                saccademtxOriginal(:,2) > stimulusStart & ...
                saccademtxOriginal(:,2) <= stimulusEnd);

            if isempty(saccadeRows)
                continue
            end

            originalOnset = saccademtxOriginal(saccadeRows,2);
            originalOffset = saccademtxOriginal(saccadeRows,4);

            % Preserve each saccade's original duration.
            saccadeDuration = originalOffset - originalOnset;

            % A 20–30 s epoch contains 80 frames at 8 fps.
            nStimulusFrames = round( ...
                (stimulusEnd - stimulusStart) * saccadeFrameRate);

            % One integer-frame circular shift per trial.
            % Thus shiftAmount is one of:
            % 0, 0.125, 0.250, ..., 9.875 seconds.
            shiftFrames = randi([0, nStimulusFrames - 1]);
            shiftAmount = shiftFrames / saccadeFrameRate;

            % Circularly shift all selected saccade onsets in this trial.
            shuffledOnset = mod( ...
                originalOnset - stimulusStart + shiftAmount, ...
                stimulusEnd - stimulusStart) + stimulusStart;

            shuffledOffset = shuffledOnset + saccadeDuration;

            % Store shuffled timestamps and matching 8-fps frame indices.
            saccademtxShuffled(saccadeRows,2) = shuffledOnset;
            saccademtxShuffled(saccadeRows,3) = nan;
            saccademtxShuffled(saccadeRows,4) = shuffledOffset;
            saccademtxShuffled(saccadeRows,5) = nan;
        end

        %% 3B. Reidentify shuffled S-T events
        % Exact event-identification logic:
        % same session, closest tail flip within +/- 0.5 s.

        simuMtxShuffled = [];
        simuN = 0;

        for saccadeN = 1:size(saccademtxShuffled,1)

            sessionN = saccademtxShuffled(saccadeN,1);
            saccadeOnset = saccademtxShuffled(saccadeN,2);
            saccadeOffset = saccademtxShuffled(saccadeN,4);

            boutRows = find( ...
                boutDetails(:,1) == sessionN & ...
                abs(boutDetails(:,2) - saccadeOnset) <= STwindow);

            if isempty(boutRows)
                continue
            end

            if numel(boutRows) > 1
                [~,closestIdx] = min(abs( ...
                    boutDetails(boutRows,2) - saccadeOnset));
                selectedBoutRow = boutRows(closestIdx);

            else
                selectedBoutRow = boutRows;
            end

            simuN = simuN + 1;

            simuMtxShuffled(simuN,1:7) = ...
                saccademtxShuffled(saccadeN,:);

            simuMtxShuffled(simuN,8:15) = ...
                boutDetails(selectedBoutRow,:);
        end

        %% 3C. Construct shuffled event-level direction and delay tables

        for ii = 1:sessionn(fishsub)

            if ismember(ii,gh.param.ExcludedSession)
                continue
            end

            if ~ismember( ...
                    gh.param.fishlog.trialdetails.trial(ii,1),stimcase)
                continue
            end

            shuffledSimuRows = find( ...
                simuMtxShuffled(:,1) == ii & ...
                simuMtxShuffled(:,2) > stimulusStart & ...
                simuMtxShuffled(:,4) <= stimulusEnd);

            for simun = 1:numel(shuffledSimuRows)

                rr = shuffledSimuRows(simun);

                boutdir = sign( ...
                    simuMtxShuffled(rr,14) * -1 - ...
                    anglebias_overall(fishsub));

                boutbias = abs( ...
                    simuMtxShuffled(rr,14) * -1 - ...
                    anglebias_overall(fishsub));

                saccadedir = simuMtxShuffled(rr,7);

                delayvec = [simuMtxShuffled(rr,9) - simuMtxShuffled(rr,2)];

                shuffleDelay(end+1,1) = delayvec; 

                shuffleDelayFishID(end+1,1) = ...
                    totalfishsub(fishsub); 

                if boutbias >= threshold_boutbias && ...
                        ismember(boutdir,[-1,1]) && ...
                        ismember(saccadedir,[-1,1])
                    isIpsi = double(saccadedir == boutdir);
                    shuffleY(end+1,1) = isIpsi; 
                    shuffleYFishID(end+1,1) = ...
                        totalfishsub(fishsub); 
                end
            end
        end
    end

    %% 3D. Fit exactly the same nested models to current shuffled events

    try

        tblShuffleDirection = table( ...
            shuffleY, ...
            categorical(shuffleYFishID), ...
            'VariableNames',{'Y','FishID'});

        glmeShuffleDirection = fitglme(tblShuffleDirection, ...
            'Y ~ 1 + (1|FishID)', ...
            'Distribution','Binomial', ...
            'Link','logit');

        betaShuffleDirection = fixedEffects(glmeShuffleDirection);
        shuffleDirectionBeta(shufflen) = betaShuffleDirection(1);

    catch ME

        warning('Direction GLME failed in shuffle %d: %s', ...
            shufflen,ME.message);
    end

    try

        tblShuffleDelay = table( ...
            shuffleDelay, ...
            categorical(shuffleDelayFishID), ...
            'VariableNames',{'Delay','FishID'});

        lmeShuffleDelay = fitlme(tblShuffleDelay, ...
            'Delay ~ 1 + (1|FishID)');

        betaShuffleDelay = fixedEffects(lmeShuffleDelay);

        shuffleDelayBeta(shufflen) = betaShuffleDelay(1);


    catch ME

        warning('Delay LME failed in shuffle %d: %s', ...
            shufflen,ME.message);
    end

    if mod(shufflen,100) == 0 || shufflen == nShuffle

        fprintf('Completed %d / %d shuffled nested models\n', ...
            shufflen,nShuffle);
    end
end

%% ------------------------------------------------------------------------
%% 4. EMPIRICAL PERMUTATION TESTS
%% ------------------------------------------------------------------------

validDirection = shuffleDirectionBeta(isfinite(shuffleDirectionBeta));
validDelay = shuffleDelayBeta(isfinite(shuffleDelayBeta));

% Direction hypothesis: real events are more ipsilateral.
% Positive GLME intercept = P(ipsilateral) > 0.5.
pDirection = (1 + sum(validDirection >= betaRealDirection(1))) / ...
    (1 + numel(validDirection));

% Delay 1 hypothesis: real tail flips occur later than real saccade onset.
% Positive intercept = T_on - S_on > 0.
pDelay = (1 + sum(validDelay >= betaRealDelay(1))) / ...
    (1 + numel(validDelay));


%% ------------------------------------------------------------------------
%% 5. PRINT RESULTS
%% ------------------------------------------------------------------------

fprintf('\n====================================================\n');
fprintf('Nested-model within-trial saccade-time shuffle test\n');
fprintf('Condition: %s\n',conditionName);
fprintf('====================================================\n');


fprintf('\nDelay 1: T_on - S_on\n');
fprintf('  Real LME intercept: %.6f s\n',betaRealDelay(1));

fprintf('  Shuffled median intercept: %.6f s\n',median(validDelay));

fprintf('  Shuffled 95%% interval: [%.6f, %.6f] s\n', ...
    prctile(validDelay,2.5), ...
    prctile(validDelay,97.5));

fprintf('  One-sided empirical P: %.6g\n',pDelay);
fprintf('\nDirection: ipsilateral S-T probability\n');
fprintf('  Real GLME intercept: %.6f\n',betaRealDirection(1));

fprintf('  Real implied ipsilateral probability: %.4f\n', ...
    1 ./ (1 + exp(-betaRealDirection(1))));

fprintf('  Shuffled median intercept: %.6f\n',median(validDirection));

fprintf('  Shuffled 95%% interval: [%.6f, %.6f]\n', ...
    prctile(validDirection,2.5), ...
    prctile(validDirection,97.5));
fprintf('  One-sided empirical P: %.6g\n',pDirection);




%% ------------------------------------------------------------------------
%% 6. SUPPLEMENTARY FIGURE: OBSERVED COORDINATION VS CIRCULAR-SHUFFLE NULL
%% ------------------------------------------------------------------------

figure(1);
clf;

set(gcf, ...
    'Color', 'w', ...
    'Position', [100 100 1040 400], ...
    'Renderer', 'painters');

shuffleColor = [0.70 0.70 0.70];
referenceColor = [0 0 0];
observedColor = plotcolor;

%% Panel a | Saccade-leading onset delay

figure(1);
hold on;

hDelay = histogram(validDelay, 35, ...
    'Normalization', 'probability', ...
    'FaceColor', shuffleColor, ...
    'EdgeColor', 'none', ...
    'FaceAlpha', 1); 

xline(median(validDelay), ':', ...
    'Color', referenceColor, ...
    'LineWidth', 1.5);

xline(0, ':', ...
    'Color', [.65 .65 .65], ...
    'LineWidth', 1.5);

xline(betaRealDelay(1), '-', ...
    'Color', observedColor, ...
    'LineWidth', 3.2);

xlim([-0.25 0.25]);
xticks(-0.25:0.05:0.25);

box off;

set(gca, ...
    'FontSize', 11, ...
    'LineWidth', 1.1, ...
    'TickDir', 'in', ...
    'Layer', 'top');

set(gcf,'Position',[100 100 400 400]);
ylim([0 0.10])
yticks([0:0.02:0.10])
yticklabels([]);
xticklabels([]);

%% Panel b | Ipsilateral eye-tail coordination

figure(2);
hold on;

shuffledIpsiPercent = 100 ./ (1 + exp(-validDirection));
observedIpsiPercent = 100 ./ (1 + exp(-betaRealDirection(1)));

hDirection = histogram(shuffledIpsiPercent, 35, ...
    'Normalization', 'probability', ...
    'FaceColor', shuffleColor, ...
    'EdgeColor', 'none', ...
    'FaceAlpha', 1); 

xline(50, ':', ...
    'Color', [.65 .65 .65], ...
    'LineWidth', 1.5);

xline(median(shuffledIpsiPercent), ':', ...
    'Color', referenceColor, ...
    'LineWidth', 1.5);

xline(observedIpsiPercent, '-', ...
    'Color', observedColor, ...
    'LineWidth', 3.2);

xlim([20 80]);
xticks(20:10:80);
ylim([0 0.10])
yticks([0:0.02:0.10])
yticklabels([]);
xticklabels([]);

box off;

set(gca, ...
    'FontSize', 11, ...
    'LineWidth', 1.1, ...
    'TickDir', 'in', ...
    'Layer', 'top');

set(gcf,'Position',[600 100 400 400]);

[min(validDelay) max(validDelay)]
[min(shuffledIpsiPercent) max(shuffledIpsiPercent)]

