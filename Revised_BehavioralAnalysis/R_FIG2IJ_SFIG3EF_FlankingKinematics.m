% this script plots event delays [figure 2ij, sfig 3fg] 
% adjust line 32 for
% different panels
%adjust lag position for p value check
clear all; close all;
global gh
LoadFishNColorSel;
timeframe = -10:11;
addpath('./stattool'); addpath('./disptool');

% per-bin containers + fish IDs
for ii = 1:length(timeframe)-1
    bout_bias{ii}    = [];
    ifsimu_bout{ii}  = [];
    bout_und_int{ii} = [];
    fishid_bout{ii}  = [];
end

for fishsub = 1:length(totalfishsub)
    gh.param.fishid = totalfishsub(fishsub);
    fm_behavim_main

    for ii = 1:sessionn(fishsub)
        if ~ismember(ii, gh.param.ExcludedSession)

            if gh.param.fishlog.trialdetails.trial(ii,1) == 8
                totalbout_n = find(gh.data.boutmtx(:,1)==ii & gh.data.boutmtx(:,4)<=40);
            else
                totalbout_n = find(gh.data.boutmtx(:,1)==ii & gh.data.boutmtx(:,4)<=20);
            end

            for boutn = 1:length(totalbout_n)
                for jj = 1:length(timeframe)-1
                    if  abs(gh.data.bout_details(totalbout_n(boutn),7)*-1- anglebias_overall(fishsub))>=0
                        if boutn >= -timeframe(jj)+1 && boutn < (length(totalbout_n)-timeframe(jj))+1

                            % bout angle bias at lag jj
                            bout_bias{jj} = [bout_bias{jj}; ...
                                abs(gh.data.bout_details(totalbout_n(boutn+timeframe(jj)),7)*-1 - ...
                                anglebias_overall(fishsub))];

                            % simu vs non-simu flag
                            ifsimu_bout{jj} = [ifsimu_bout{jj}; ...
                                ~isempty(find( ...
                                (gh.data.simuMtx(:,1)==gh.data.bout_details(totalbout_n(boutn),1)) & ...
                                (gh.data.simuMtx(:,9)==gh.data.bout_details(totalbout_n(boutn),2))))];

                            % fish ID for this event
                            fishid_bout{jj} = [fishid_bout{jj}; totalfishsub(fishsub)];

                            % undulation interval (absolute)
                            X = gh.data.anglevect( ...
                                gh.data.bout_details(totalbout_n(boutn+timeframe(jj)),3): ...
                                gh.data.bout_details(totalbout_n(boutn+timeframe(jj)),5));
                            X0 = X - max(X);
                            X0_inv = -X0;
                            pks_X0 = findpeaks(X0);
                            pks_X0_inv = findpeaks(X0_inv);

                            und_int = gh.data.bout_details(totalbout_n(boutn+timeframe(jj)),6) ./ ...
                                ((length(pks_X0)+length(pks_X0_inv))./2);
                            bout_und_int{jj} = [bout_und_int{jj}; und_int];

                        end

                    end
                end
            end
        end
    end
end

% ---------- Figure 3b/c/d: bout_bias time course ----------

mid = []; u = []; l = [];
p_val_boutdynamics = nan(2, length(timeframe)-1);

for ii = 1:length(timeframe)-1

    displaydata = bout_bias{ii};
    displaysimu = logical(ifsimu_bout{ii});


    mid(ii,:) = [nanmedian(displaydata(~displaysimu)), ...
        nanmedian(displaydata(displaysimu))];

    u(ii,:)   = [quantile(displaydata(~displaysimu),0.75), ...
        quantile(displaydata(displaysimu),0.75)];
    l(ii,:)   = [quantile(displaydata(~displaysimu),0.25), ...
        quantile(displaydata(displaysimu),0.25)];

    % original per-bin event-level ranksum (diagnostic only)
    p_val_boutdynamics(1,ii) = ranksum(displaydata(~displaysimu), ...
        displaydata(displaysimu));
end
   
figure(1)
ymax = 8;ymin = 1;
plot((1:length(timeframe)-1)', mid(:,1), 'color', independentColorT, ...
    'linestyle', ':', 'LineWidth', 2.5); hold on
plot((1:length(timeframe)-1)', mid(:,2), 'color', STColor, ...
    'linestyle', ':', 'LineWidth', 2.5); hold on
plot((1:length(timeframe)-1)', smooth(mid(:,1),3), 'color', independentColorT, ...
    'linestyle', '-', 'LineWidth', 4); hold on
plot((1:length(timeframe)-1)', smooth(mid(:,2),3), 'color', STColor, ...
    'linestyle', '-', 'LineWidth', 4); hold on
plot((1:length(timeframe)-1)', max(min(u(:,1),ymax),ymin), 'color', independentColorT, ...
    'linestyle', ':', 'LineWidth', 2.5); hold on
plot((1:length(timeframe)-1)', max(min(u(:,2),ymax),ymin), 'color', STColor, ...
    'linestyle', ':', 'LineWidth', 2.5); hold on
plot((1:length(timeframe)-1)', max(min(l(:,1),ymax),ymin), 'color', independentColorT, ...
    'linestyle', ':', 'LineWidth', 2.5); hold on
plot((1:length(timeframe)-1)', max(min(l(:,2),ymax),ymin), 'color', STColor, ...
    'linestyle', ':', 'LineWidth', 2.5); hold on

xticks(1:length(timeframe)-1); xlim([1, length(timeframe)-1])
xticklabels([]); yticklabels([]); box off
ylim([1 8]); yticks(1:0.5:8)
set(gcf, 'Position', [300 300 300 400])

% ---------- bout_bias: pooled 10–12 and lag 11 LMEs ----------
% pooled lags 10–12
pool_bias_10_12 = [];
pool_simu_10_12 = [];
pool_fish_10_12 = [];

for ii_pool_10_12 = 10:12
    bias_pool_10_12 = bout_bias{ii_pool_10_12};
    simu_pool_10_12 = logical(ifsimu_bout{ii_pool_10_12});
    fish_pool_10_12 = fishid_bout{ii_pool_10_12};

    pool_bias_10_12 = [pool_bias_10_12; double(bias_pool_10_12(:))];
    pool_simu_10_12 = [pool_simu_10_12; double(simu_pool_10_12(:))];
    pool_fish_10_12 = [pool_fish_10_12; fish_pool_10_12(:)];
end

tbl_bias_10_12 = table(pool_bias_10_12(:), pool_simu_10_12(:), ...
    categorical(pool_fish_10_12(:)), ...
    'VariableNames', {'Response','Simu','Fish'});
lme_bias_10_12   = fitlme(tbl_bias_10_12, 'Response ~ Simu + (1|Fish)');
coefTbl_bias_10_12 = lme_bias_10_12.Coefficients;
pval_bias_10_12    = coefTbl_bias_10_12.pValue(strcmp(coefTbl_bias_10_12.Name,'Simu'));
disp(['bout_bias pooled lags 10–12 (clueless case) LME p = ', num2str(pval_bias_10_12)]);

clear ii_pool_10_12 bias_pool_10_12 simu_pool_10_12 fish_pool_10_12 ...
    tbl_bias_10_12 lme_bias_10_12 coefTbl_bias_10_12 ...
    pool_bias_10_12 pool_simu_10_12 pool_fish_10_12

% ---------- bout_bias: pooled lags 2–4 LMEs ----------
% pooled lags 2–4
pool_bias_5_7 = [];
pool_simu_5_7 = [];
pool_fish_5_7 = [];

for ii_pool_5_7 = 5:7
    bias_pool_5_7 = bout_bias{ii_pool_5_7};
    simu_pool_5_7 = logical(ifsimu_bout{ii_pool_5_7});
    fish_pool_5_7 = fishid_bout{ii_pool_5_7};

    pool_bias_5_7 = [pool_bias_5_7; double(bias_pool_5_7(:))];
    pool_simu_5_7 = [pool_simu_5_7; double(simu_pool_5_7(:))];
    pool_fish_5_7 = [pool_fish_5_7; fish_pool_5_7(:)];
end

tbl_bias_5_7 = table(pool_bias_5_7(:), pool_simu_5_7(:), ...
    categorical(pool_fish_5_7(:)), ...
    'VariableNames', {'Response','Simu','Fish'});

lme_bias_5_7 = fitlme(tbl_bias_5_7, 'Response ~ Simu + (1|Fish)');
coefTbl_bias_5_7 = lme_bias_5_7.Coefficients;
pval_bias_5_7 = coefTbl_bias_5_7.pValue(strcmp(coefTbl_bias_5_7.Name, 'Simu'));

disp(['bout_bias pooled lags 5–7 (clueless case) LME p = ', num2str(pval_bias_5_7)]);

clear ii_pool_5_7 bias_pool_5_7 simu_pool_5_7 fish_pool_5_7 ...
    tbl_bias_5_7 lme_bias_5_7 coefTbl_bias_5_7 ...
    pool_bias_5_7 pool_simu_5_7 pool_fish_5_7;


% single lag 11
ii_target = 11;
resp_bias_11 = bout_bias{ii_target};
simu_bias_11 = logical(ifsimu_bout{ii_target});
fish_bias_11 = fishid_bout{ii_target};

tbl_bias_11 = table(double(resp_bias_11(:)), double(simu_bias_11(:)), ...
    categorical(fish_bias_11(:)), ...
    'VariableNames', {'Response','Simu','Fish'});
lme_bias_11   = fitlme(tbl_bias_11, 'Response ~ Simu + (1|Fish)');
coefTbl_bias_11 = lme_bias_11.Coefficients;
pval_bias_11    = coefTbl_bias_11.pValue(strcmp(coefTbl_bias_11.Name,'Simu'));
disp(['bout_bias lag 11 (clueless case) LME p = ', num2str(pval_bias_11)]);

clear ii_target resp_bias_11 simu_bias_11 fish_bias_11 ...
    tbl_bias_11 lme_bias_11 coefTbl_bias_11

% ---------- Figure 4h/i: bout_und_int time course (absolute) ----------
figure(2)
mid2 = []; u2 = []; l2 = [];

for ii = 1:length(timeframe)-1
    displaydata = bout_und_int{ii};
    displaysimu = logical(ifsimu_bout{ii});

    mid2(ii,:) = [nanmedian(displaydata(~displaysimu)), ...
        nanmedian(displaydata(displaysimu))];
    u2(ii,:)   = [quantile(displaydata(~displaysimu),0.75), ...
        quantile(displaydata(displaysimu),0.75)];
    l2(ii,:)   = [quantile(displaydata(~displaysimu),0.25), ...
        quantile(displaydata(displaysimu),0.25)];

    p_val_boutdynamics(2,ii) = ranksum(displaydata(~displaysimu), ...
        displaydata(displaysimu));
end

ymax = 0.048;ymin = 0.036;

plot((1:length(timeframe)-1)', mid2(:,1), 'color', independentColorT, ...
    'linestyle', ':', 'LineWidth', 2.5); hold on
plot((1:length(timeframe)-1)', mid2(:,2), 'color', STColor, ...
    'linestyle', ':', 'LineWidth', 2.5); hold on
plot((1:length(timeframe)-1)', smooth(mid2(:,1),3), 'color', independentColorT, ...
    'linestyle', '-', 'LineWidth', 4); hold on
plot((1:length(timeframe)-1)', smooth(mid2(:,2),3), 'color', STColor, ...
    'linestyle', '-', 'LineWidth', 4); hold on
plot((1:length(timeframe)-1)', max(min(u2(:,1), ymax), ymin), 'color', independentColorT, ...
    'linestyle', ':', 'LineWidth', 2.5); hold on
plot((1:length(timeframe)-1)', max(min(u2(:,2), ymax), ymin), 'color', STColor, ...
    'linestyle', ':', 'LineWidth', 2.5); hold on
plot((1:length(timeframe)-1)', max(min(l2(:,1), ymax), ymin), 'color', independentColorT, ...
    'linestyle', ':', 'LineWidth', 2.5); hold on
plot((1:length(timeframe)-1)', max(min(l2(:,2), ymax), ymin), 'color', STColor, ...
    'linestyle', ':', 'LineWidth', 2.5); hold on

xticks(1:length(timeframe)-1); xlim([1, length(timeframe)-1])
xticklabels([]); yticklabels([]); box off
ylim([0.036 0.048]); yticks(0.036:0.001:0.048);
set(gcf, 'Position', [600 300 300 400])

% ---------- undulation: pooled 10–12 and lag 11 LMEs ----------
% pooled lags 10–12
pool_und_10_12   = [];
pool_simuU_10_12 = [];
pool_fishU_10_12 = [];

for ii_poolU_10_12 = 5:7
    und_pool_10_12   = bout_und_int{ii_poolU_10_12};
    simuU_pool_10_12 = logical(ifsimu_bout{ii_poolU_10_12});
    fishU_pool_10_12 = fishid_bout{ii_poolU_10_12};

    pool_und_10_12   = [pool_und_10_12;   double(und_pool_10_12(:))];
    pool_simuU_10_12 = [pool_simuU_10_12; double(simuU_pool_10_12(:))];
    pool_fishU_10_12 = [pool_fishU_10_12; fishU_pool_10_12(:)];
end

tbl_und_10_12 = table(pool_und_10_12(:), pool_simuU_10_12(:), ...
    categorical(pool_fishU_10_12(:)), ...
    'VariableNames', {'Response','Simu','Fish'});
lme_und_10_12   = fitlme(tbl_und_10_12, 'Response ~ Simu + (1|Fish)');
coefTbl_und_10_12 = lme_und_10_12.Coefficients;
pval_und_10_12    = coefTbl_und_10_12.pValue(strcmp(coefTbl_und_10_12.Name,'Simu'));
disp(['undulation pooled lags 10–12 (clueless case) LME p = ', num2str(pval_und_10_12)]);

clear ii_poolU_10_12 und_pool_10_12 simuU_pool_10_12 fishU_pool_10_12 ...
    pool_und_10_12 pool_simuU_10_12 pool_fishU_10_12 ...
    tbl_und_10_12 lme_und_10_12 coefTbl_und_10_12

% single lag 11
ii_target = 11;
resp_und_11  = bout_und_int{ii_target};
simu_und_11  = logical(ifsimu_bout{ii_target});
fish_und_11  = fishid_bout{ii_target};

tbl_und_11 = table(double(resp_und_11(:)), double(simu_und_11(:)), ...
    categorical(fish_und_11(:)), ...
    'VariableNames', {'Response','Simu','Fish'});
lme_und_11   = fitlme(tbl_und_11, 'Response ~ Simu + (1|Fish)');
coefTbl_und_11 = lme_und_11.Coefficients;
pval_und_11    = coefTbl_und_11.pValue(strcmp(coefTbl_und_11.Name,'Simu'));
disp(['undulation lag 11 (clueless case) LME p = ', num2str(pval_und_11)]);

clear ii_target resp_und_11 simu_und_11 fish_und_11 ...
    tbl_und_11 lme_und_11 coefTbl_und_11