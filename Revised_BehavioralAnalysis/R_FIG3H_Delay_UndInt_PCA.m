%% R_FIG3H_Delay_UndInt_PCA.m
clear all; close all;
global gh
LoadFishNColorSel;
addpath('./stattool');
addpath('./disptool');

all_und    = [];   % undulation interval metric (s)
all_x      = [];   % temporal predictors x1..x5
all_group3 = [];   % 1=blank, 2=app, 3=av
all_fishID = [];

%% ---------- build event-level data (20–30 s S–T events) ----------
for fishsub = 1:length(totalfishsub)
    gh.param.fishid = totalfishsub(fishsub);
    fm_behavim_main

    for ii = 1:sessionn(fishsub)
        if ~ismember(ii, gh.param.ExcludedSession)

            stimsimu_n = find(gh.data.simuMtx(:,1) == ii & ...
                gh.data.simuMtx(:,2) > 20 & ...
                gh.data.simuMtx(:,4) <= 30);
            prestim_boutn = find(gh.data.boutmtx(:,1) == ii & ...
                gh.data.boutmtx(:,2) > 10 & ...
                gh.data.boutmtx(:,4) <= 20);

            trcode = gh.param.fishlog.trialdetails.trial(ii,1);
            if     trcode == 8
                group3 = 1;   % blank
            elseif trcode == 1 || trcode == 2
                group3 = 2;   % appetitive
            elseif trcode == 4 || trcode == 5
                group3 = 3;   % aversive
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

            for jj = 1:length(stimsimu_n)
                irow = stimsimu_n(jj);

                % temporal predictors x1..x5
                x1 = gh.data.simuMtx(irow,9)  - gh.data.simuMtx(irow,2);
                x2 = gh.data.simuMtx(irow,11) - gh.data.simuMtx(irow,4);
                x3 = gh.data.simuMtx(irow,4)  - gh.data.simuMtx(irow,9);
                x4 = -gh.data.simuMtx(irow,2) + gh.data.simuMtx(irow,4);
                x5 = -gh.data.simuMtx(irow,9) + gh.data.simuMtx(irow,11);

                % undulation interval metric
                X      = gh.data.anglevect(gh.data.simuMtx(irow,10):gh.data.simuMtx(irow,12));
                X0     = X - max(X);
                X0_inv = -X0;
                pks_X0     = findpeaks(X0);
                pks_X0_inv = findpeaks(X0_inv);

                n_halfcycles = (length(pks_X0) + length(pks_X0_inv))./2;

                if n_halfcycles >=4
                    undul_interval = gh.data.simuMtx(irow,13) ./ n_halfcycles - avgval - boutundint_overall(fishsub);   % seconds
                else
                    undul_interval = NaN;
                end

                newentry = [x1 x2 x3 x4 x5 undul_interval];

                % long-format arrays
                all_x      = [all_x;      newentry([1,2,3,4,5])];
                all_und    = [all_und;    newentry(6)];
                all_group3 = [all_group3; group3];
                all_fishID = [all_fishID; gh.param.fishid];
            end
        end
    end
end

%% ---------- optional filter on x1 ----------
all_keep   = all_x(:,1) >= -25/400;
all_x      = all_x(all_keep,:);
all_und    = all_und(all_keep);
all_group3 = all_group3(all_keep);
all_fishID = all_fishID(all_keep);

%% ---------- remove invalid rows ----------
valid_keep = ~isnan(all_und) & all(isfinite(all_x),2);
% all_x      = all_x(valid_keep,:);
% all_und    = all_und(valid_keep);
% all_group3 = all_group3(valid_keep);
% all_fishID = all_fishID(valid_keep);

%% ---------- PCA on temporal predictors x1..x5 ----------
[coeff, score, latent, tsq, explained] = pca(all_x);
PC1 = score(:,1);
PC2 = score(:,2);

fprintf('Explained variance by PCs (%%): PC1=%.1f, PC2=%.1f, cum(1:2)=%.1f\n', ...
    explained(1), explained(2), sum(explained(1:2)));

%% ---------- Mixed model: undulation vs PC1/PC2 ----------
Tpc = table(all_und, PC1, PC2, ...
    categorical(all_group3), ...
    categorical(all_fishID), ...
    'VariableNames', {'und','PC1','PC2','group','fishID'});

lme_und_PC12 = fitlme(Tpc, 'und ~ 1 + PC1*group + PC2*group + (1|fishID)');

%% ---------- read p-values: ANOVA + fixed effects ----------
disp('=== UNDULATION model ANOVA ===');
tbl_und = anova(lme_und_PC12);
% disp(tbl_und);

%% ---------- slope tests and contrasts ----------
names_und = lme_und_PC12.CoefficientNames;

i_PC1_und   = find(strcmp(names_und,'PC1'));
i_PC1g2_und = find(strcmp(names_und,'PC1:group_2'));
i_PC1g3_und = find(strcmp(names_und,'PC1:group_3'));

% per-group PC1 slopes
H = zeros(1,length(names_und));
H(i_PC1_und) = 1;
[p_g1_und, F_g1_und] = coefTest(lme_und_PC12,H);

H = zeros(1,length(names_und));
H(i_PC1_und) = 1;
H(i_PC1g2_und) = 1;
[p_g2_und, F_g2_und] = coefTest(lme_und_PC12,H);

H = zeros(1,length(names_und));
H(i_PC1_und) = 1;
H(i_PC1g3_und) = 1;
[p_g3_und, F_g3_und] = coefTest(lme_und_PC12,H);

fprintf('UND:   PC1 slope p-values: blank=%.4g, app=%.4g, av=%.4g\n', ...
    p_g1_und, p_g2_und, p_g3_und);


% app vs av PC1 slope
H_app_vs_av_und_PC1 = zeros(1, numel(names_und));
H_app_vs_av_und_PC1(i_PC1g2_und) =  1;
H_app_vs_av_und_PC1(i_PC1g3_und) = -1;
[p_app_vs_av_und_PC1, F_app_vs_av_und_PC1] = coefTest(lme_und_PC12, H_app_vs_av_und_PC1);

fprintf('UND:   app vs av PC1 slope p=%.4g, F=%.4g\n', ...
    p_app_vs_av_und_PC1, F_app_vs_av_und_PC1);

%% ---------- colors and setup ----------
blankColor = [0.5 0.5 0.5];
appeColor  = [119 136 172]./255;
averColor  = [238 129 114]./255;

g1 = all_group3 == 1;
g2 = all_group3 == 2;
g3 = all_group3 == 3;

typicalFish = all_fishID(1);
alphaCI = 0.05;

%% ---------- PC1 range with padding ----------
pc1_min = min(PC1);
pc1_max = max(PC1);
pad = 0.02 * (pc1_max - pc1_min);
pc1_range = linspace(pc1_min - pad, pc1_max + pad, 50)';

Tpred_PC1 = table( ...
    nan(size(pc1_range)), ...
    pc1_range, ...
    zeros(size(pc1_range)), ...
    categorical(ones(size(pc1_range))), ...
    categorical(repmat(typicalFish,size(pc1_range))), ...
    'VariableNames', {'und','PC1','PC2','group','fishID'});

%% ---------- figure: undulation interval vs PC1 ----------
figure('Name','Undulation Interval vs PC1','Position',[320 550 200 400]);
hold on;

% visible zero reference line
yline(0, ':', 'Color', [0.6 0.6 0.6], 'LineWidth', 1);

% scatter, same style as angle plot
scatter(PC1(g1), all_und(g1)*1000, 20, blankColor, 'filled', 'MarkerFaceAlpha',0.4);
scatter(PC1(g2), all_und(g2)*1000, 20, appeColor,  'filled', 'MarkerFaceAlpha',0.4);
scatter(PC1(g3), all_und(g3)*1000, 20, averColor,  'filled', 'MarkerFaceAlpha',0.4);

for grp = 1:3
    Tpred = Tpred_PC1;
    Tpred.group = categorical(grp*ones(size(pc1_range)));
    [yhat_und, CI_und] = predict(lme_und_PC12, Tpred, 'Conditional', false, 'Alpha', alphaCI);

    switch grp
        case 1
            col = blankColor;
        case 2
            col = appeColor;
        case 3
            col = averColor;
    end

    yhat_und_ms = yhat_und * 1000;
    CI_und_ms   = CI_und * 1000;

    fill([pc1_range; flipud(pc1_range)], ...
        [CI_und_ms(:,1); flipud(CI_und_ms(:,2))], ...
        col, 'FaceAlpha',0.05, 'EdgeColor','none');
    if grp==3
        plot(pc1_range, yhat_und_ms, '-', 'Color', col, 'LineWidth', 2);
    elseif grp==2
        plot(pc1_range, yhat_und_ms, ':', 'Color', col, 'LineWidth', 2);
    else
        plot(pc1_range, yhat_und_ms, '--', 'Color', col, 'LineWidth', 2);
    end
end

% ylim([0 80]);  % includes y=0 so the reference line is visible
xlim([pc1_min - pad, pc1_max + pad]);
ylabel('');yticklabels({});
xlabel('');xticklabels({});
box off;title('');
ylim([-20 20]);
set(gca,'XTickMode','auto','YTickMode','auto');