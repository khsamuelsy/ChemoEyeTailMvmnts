%% R_FIG3G_Delay_Bias_PCA.m
clear all; close all;
global gh
LoadFishNColorSel;
addpath('./stattool');
addpath('./disptool');

all_ang    = [];   % angle bias metric
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

            trcode = gh.param.fishlog.trialdetails.trial(ii,1);
            if     trcode == 8
                group3 = 1;   % blank
            elseif trcode == 1 || trcode == 2
                group3 = 2;   % appetitive
            elseif trcode == 4 || trcode == 5
                group3 = 3;   % aversive
            end

            % prestim angle-bias baseline (10–20 s bouts)
            prestim_boutn = find(gh.data.boutmtx(:,1) == ii & ...
                gh.data.boutmtx(:,2) > 10 & ...
                gh.data.boutmtx(:,4) <= 20);

            sumval = 0;
            if ~isempty(prestim_boutn)
                for jj = 1:length(prestim_boutn)
                    sumval = sumval + abs(gh.data.bout_details(prestim_boutn(jj),7)*(-1) ...
                        - anglebias_overall(fishsub));
                end
                avgval = sumval ./ length(prestim_boutn);
            else
                avgval = 0;
            end

            for jj = 1:length(stimsimu_n)

                % for jj = 1:min(length(stimsimu_n),1)
                irow = stimsimu_n(jj);

                % temporal predictors x1..x5
                x1 = gh.data.simuMtx(irow,9)  - gh.data.simuMtx(irow,2);
                x2 = gh.data.simuMtx(irow,11) - gh.data.simuMtx(irow,4);
                x3 = gh.data.simuMtx(irow,4)  - gh.data.simuMtx(irow,9);
                x4 = -gh.data.simuMtx(irow,2) + gh.data.simuMtx(irow,4);
                x5 = -gh.data.simuMtx(irow,9) + gh.data.simuMtx(irow,11);

                % angle bias metric
                ang_metric = abs(gh.data.simuMtx(irow,14)*-1 - anglebias_overall(fishsub)) - avgval;

                newentry = [x1 x2 x3 x4 x5 ang_metric];

                % long-format arrays
                all_x      = [all_x;      newentry([1,2,3,4,5])];
                all_ang    = [all_ang;    newentry(6)];
                all_group3 = [all_group3; group3];
                all_fishID = [all_fishID; gh.param.fishid];
            end
        end
    end
end

%% ---------- optional filter on x1 ----------
all_keep   = all_x(:,1) >= -25/400;
all_x      = all_x(all_keep,:);
all_ang    = all_ang(all_keep);
all_group3 = all_group3(all_keep);
all_fishID = all_fishID(all_keep);

%% ---------- PCA on temporal predictors x1..x5 ----------
[coeff, score, latent, tsq, explained] = pca(all_x);
PC1 = score(:,1);
PC2 = score(:,2);

fprintf('Explained variance by PCs (%%): PC1=%.1f, PC2=%.1f, cum(1:2)=%.1f\n', ...
        explained(1), explained(2), sum(explained(1:2)));

%% ---------- Mixed model: angle vs PC1/PC2 ----------
Tpc = table(all_ang, PC1, PC2, ...
            categorical(all_group3), ...
            categorical(all_fishID), ...
            'VariableNames', {'ang','PC1','PC2','group','fishID'});

lme_ang_PC12 = fitlme(Tpc, 'ang ~ 1 + PC1*group + PC2*group + (1|fishID)');

%% ---------- read p-values: ANOVA + fixed effects ----------
disp('=== ANGLE model ANOVA ===');
tbl_ang = anova(lme_ang_PC12);
% disp(tbl_ang);

%% ---------- per-group PC1 slope tests ----------
names_ang = lme_ang_PC12.CoefficientNames;
i_PC1_ang   = find(strcmp(names_ang,'PC1'));
i_PC1g2_ang = find(strcmp(names_ang,'PC1:group_2'));
i_PC1g3_ang = find(strcmp(names_ang,'PC1:group_3'));

H = zeros(1,length(names_ang));
H(i_PC1_ang) = 1;
[p_g1_ang, F_g1_ang] = coefTest(lme_ang_PC12,H);

H = zeros(1,length(names_ang));
H(i_PC1_ang) = 1;
H(i_PC1g2_ang) = 1;
[p_g2_ang, F_g2_ang] = coefTest(lme_ang_PC12,H);

H = zeros(1,length(names_ang));
H(i_PC1_ang) = 1;
H(i_PC1g3_ang) = 1;
[p_g3_ang, F_g3_ang] = coefTest(lme_ang_PC12,H);

fprintf('ANGLE: PC1 slope p-values: blank=%.4g, app=%.4g, av=%.4g\n', ...
        p_g1_ang, p_g2_ang, p_g3_ang);

%% ---------- app vs av PC1 slope only ----------
H_app_vs_av_ang_PC1 = zeros(1, numel(names_ang));
H_app_vs_av_ang_PC1(i_PC1g2_ang) =  1;
H_app_vs_av_ang_PC1(i_PC1g3_ang) = -1;
[p_app_vs_av_ang_PC1, F_app_vs_av_ang_PC1] = coefTest(lme_ang_PC12, H_app_vs_av_ang_PC1);

fprintf('ANGLE: app vs av PC1 slope p=%.4g, F=%.4g\n', ...
        p_app_vs_av_ang_PC1, F_app_vs_av_ang_PC1);

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
    'VariableNames', {'ang','PC1','PC2','group','fishID'});

%% ---------- figure: angle bias vs PC1 ----------
figure('Name','Angle Bias vs PC1','Position',[320 550 200 400]);
hold on;

yline(0, ':', 'Color', [0 0 0], 'LineWidth', 1);

scatter(PC1(g1), all_ang(g1), 20, blankColor, 'filled', 'MarkerFaceAlpha',0.4);
scatter(PC1(g2), all_ang(g2), 20, appeColor,  'filled', 'MarkerFaceAlpha',0.4);
scatter(PC1(g3), all_ang(g3), 20, averColor,  'filled', 'MarkerFaceAlpha',0.4);

for grp = 1:3
    Tpred = Tpred_PC1;
    Tpred.group = categorical(grp*ones(size(pc1_range)));
    [yhat_ang, CI_ang] = predict(lme_ang_PC12, Tpred, 'Conditional', false, 'Alpha', alphaCI);
    
    switch grp
        case 1
            col = blankColor;
            lstyle = '--';   % blank group dotted
        case 2
            col = appeColor;
            lstyle = ':';
        case 3
            col = averColor;
            lstyle = '-';
    end
    
    fill([pc1_range; flipud(pc1_range)], ...
         [CI_ang(:,1); flipud(CI_ang(:,2))], ...
         col, 'FaceAlpha',0.05, 'EdgeColor','none');
    plot(pc1_range, yhat_ang, lstyle, 'Color', col, 'LineWidth', 2);
end

ylim([-20 20]);
xlim([pc1_min - pad, pc1_max + pad]);
ylabel('');yticklabels({});
xlabel('');xticklabels({});
box off;title('');
set(gca,'XTickMode','auto','YTickMode','auto');