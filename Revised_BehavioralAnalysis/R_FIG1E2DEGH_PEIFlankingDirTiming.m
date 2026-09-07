clear all; close all; LoadFishNColorSel;

global gh

addpath('./stattool'); addpath('./disptool');

% global PEI container
for ii=1:2
    PEInterval_Eventbased{ii} = [];
    PEInterval_Fishbased{ii} = [];
end

timeframe = [-10:11];
timeframe_saccade = [-3:4];

% event-level containers + fish ID trackers for mixed models
for ii=1:length(timeframe)-1
    bout_clock{ii}=[];          ifsimu_bout_clock{ii}=[];      fishid_bout_clock{ii}=[];
    bout_dir{ii}=[];            ifsimu_bout_dir{ii}=[];        fishid_bout_dir{ii}=[];
end

for ii=1:length(timeframe_saccade)-1
    saccade_clock{ii}=[];       ifsimu_saccade_clock{ii}=[];   fishid_saccade_clock{ii}=[];
    saccade_dir{ii}=[];         ifsimu_saccade_dir{ii}=[];     fishid_saccade_dir{ii}=[];
end

for fishsub = 1:length(totalfishsub)
    gh.param.fishid = totalfishsub(fishsub);
    fm_behavim_main
    
    % a local container for storing pre-saccade and pre-tail flip interval for each fish
    for ii=1:2
        PEInterval{ii} = [];
    end

    for ii=1:sessionn(fishsub)
        if ~ismember(ii,gh.param.ExcludedSession)
            if gh.param.fishlog.trialdetails.trial(ii,1)==8
                totalbout_n = find(gh.data.boutmtx(:,1)==ii & gh.data.boutmtx(:,4)<=40);
                totalsaccade_n = find(gh.data.saccademtx(:,1)==ii & gh.data.saccademtx(:,4)<=40);
            else
                totalbout_n = find(gh.data.boutmtx(:,1)==ii & gh.data.boutmtx(:,4)<=20);
                totalsaccade_n = find(gh.data.saccademtx(:,1)==ii & gh.data.saccademtx(:,4)<=20);
            end

            % saccade
            for saccaden = 1:length(totalsaccade_n)
                if saccaden>1
                    PEInterval{1} = [PEInterval{1}; ...
                        (gh.data.saccademtx(totalsaccade_n(saccaden),2)- ...
                         gh.data.saccademtx(totalsaccade_n(saccaden-1),2))];
                end

                for jj=1:length(timeframe_saccade)-1
                    if saccaden>-timeframe_saccade(jj)+1 && ...
                            saccaden <(length(totalsaccade_n)-timeframe_saccade(jj))+1
                        dt = gh.data.saccademtx(totalsaccade_n(saccaden+timeframe_saccade(jj)),2)- ...
                             gh.data.saccademtx(totalsaccade_n(saccaden+timeframe_saccade(jj)-1),2);
                        isSim = ~isempty(find((gh.data.simuMtx(:,1)==gh.data.saccademtx(totalsaccade_n(saccaden),1)) & ...
                                              (gh.data.simuMtx(:,2)==gh.data.saccademtx(totalsaccade_n(saccaden),2))));
                        saccade_clock{jj}        = [saccade_clock{jj}; dt];
                        ifsimu_saccade_clock{jj} = [ifsimu_saccade_clock{jj}; isSim];
                        fishid_saccade_clock{jj} = [fishid_saccade_clock{jj}; totalfishsub(fishsub)];
                    end

                    if saccaden>=-timeframe_saccade(jj)+1 && ...
                            saccaden <(length(totalsaccade_n)-timeframe_saccade(jj))+1
                        same_dir = gh.data.saccademtx(totalsaccade_n(saccaden+timeframe_saccade(jj)),7)== ...
                                   gh.data.saccademtx(totalsaccade_n(saccaden),7);
                        isSim = ~isempty(find((gh.data.simuMtx(:,1)==gh.data.saccademtx(totalsaccade_n(saccaden),1)) & ...
                                              (gh.data.simuMtx(:,2)==gh.data.saccademtx(totalsaccade_n(saccaden),2))));
                        saccade_dir{jj}         = [saccade_dir{jj}; same_dir];
                        ifsimu_saccade_dir{jj}  = [ifsimu_saccade_dir{jj}; isSim];
                        fishid_saccade_dir{jj}  = [fishid_saccade_dir{jj}; totalfishsub(fishsub)];
                    end
                end
            end

            % bout
            for boutn = 1:length(totalbout_n)
                if boutn>1
                    PEInterval{2} = [PEInterval{2}; ...
                        (gh.data.bout_details(totalbout_n(boutn),2)- ...
                         gh.data.bout_details(totalbout_n(boutn-1),2))];
                end

                for jj=1:length(timeframe)-1
                    if boutn>-timeframe(jj)+1 && ...
                            boutn <(length(totalbout_n)-timeframe(jj))+1
                        dt = gh.data.bout_details(totalbout_n(boutn+timeframe(jj)),2)- ...
                             gh.data.bout_details(totalbout_n(boutn+timeframe(jj)-1),2);
                        isSim = ~isempty(find((gh.data.simuMtx(:,1)==gh.data.bout_details(totalbout_n(boutn),1)) & ...
                                              (gh.data.simuMtx(:,9)==gh.data.bout_details(totalbout_n(boutn),2))));
                        bout_clock{jj}        = [bout_clock{jj}; dt];
                        ifsimu_bout_clock{jj} = [ifsimu_bout_clock{jj}; isSim];
                        fishid_bout_clock{jj} = [fishid_bout_clock{jj}; totalfishsub(fishsub)];
                    end

                    if abs(gh.data.bout_details(totalbout_n(boutn),7)*-1-anglebias_overall(fishsub))>=1.5
                        if boutn>=-timeframe(jj)+1 && ...
                                boutn <(length(totalbout_n)-timeframe(jj))+1 && ...
                                abs(gh.data.bout_details(totalbout_n(boutn+timeframe(jj)),7)*-1-anglebias_overall(fishsub))>=1.5
                            same_dir = (sign(gh.data.bout_details(totalbout_n(boutn+timeframe(jj)),7)*-1-anglebias_overall(fishsub))== ...
                                        sign(gh.data.bout_details(totalbout_n(boutn),7)*-1-anglebias_overall(fishsub)));
                            isSim = ~isempty(find((gh.data.simuMtx(:,1)==gh.data.bout_details(totalbout_n(boutn),1)) & ...
                                                  (gh.data.simuMtx(:,9)==gh.data.bout_details(totalbout_n(boutn),2))));
                            bout_dir{jj}        = [bout_dir{jj}; same_dir];
                            ifsimu_bout_dir{jj} = [ifsimu_bout_dir{jj}; isSim];
                            fishid_bout_dir{jj} = [fishid_bout_dir{jj}; totalfishsub(fishsub)];
                        end
                    end
                end
            end
        end
        clear totalbout_n totalsaccade_n
    end

    PEInterval_Eventbased{1} = [PEInterval_Eventbased{1};PEInterval{1}];
    PEInterval_Eventbased{2} = [PEInterval_Eventbased{2};PEInterval{2}];
    PEInterval_Fishbased{1} = [PEInterval_Fishbased{1};median(PEInterval{1})];
    PEInterval_Fishbased{2} = [PEInterval_Fishbased{2};median(PEInterval{2})];
end

%% Figure 1E left  (unchanged)
figure(1)
for ii=1:2
    dotcolor = [0.5 0.5 0.5];
    xpos = ii;
    limiter_vec = 0:0.125:10;
    limiter_vec_indx =[];prop=[];
    tempdata=min(PEInterval_Eventbased{ii},10);

    for jj=1:length(PEInterval_Eventbased{ii})
        [~, limiter_vec_indx(jj)] = min(abs(tempdata(jj)-limiter_vec));
    end

    for jj=1:length(limiter_vec)
        prop(jj)=length(find(limiter_vec_indx==jj))./length(PEInterval_Eventbased{ii});
    end

    weight = prop(limiter_vec_indx)*4.5;
    randx = xpos+rand(1,length(PEInterval_Eventbased{ii})).*weight-weight./2;
    displaydata = PEInterval_Eventbased{ii};

    scatter(randx,min(displaydata,10),10,'MarkerFaceColor',dotcolor,'MarkerEdgeColor','none','MarkerFaceAlpha',.5); hold on
    mid(ii) = nanmedian(displaydata);
    if ii==2
        display(['Median PEI for saccade and tail flip : ',num2str(mid(1)),' and ',num2str(mid(2)) ,' seconds'])
    end
    u(ii) = quantile(displaydata,0.75);
    l(ii) = quantile(displaydata,0.25);

    plot([xpos-0.15.*1.5 xpos+0.15.*1.5],[mid(ii) mid(ii)],'k','linewidth',1.5); hold on
    plot([xpos-0.06.*1.5 xpos+0.06.*1.5],[u(ii) u(ii)],'k','linewidth',1.5); hold on
    plot([xpos-0.06.*1.5 xpos+0.06.*1.5],[l(ii) l(ii)],'k','linewidth',1.5); hold on

    if ii==1, scaling_factor = 0.9; else, scaling_factor = 0.32; end
    violin(xpos,min(PEInterval_Eventbased{ii},10),'facecolor',dotcolor, ...
        'scaling',scaling_factor,'facealpha',0.2,'style',2); hold on;
    ylim([0 10])
end

p_EventbasedPEI = ranksum(PEInterval_Eventbased{1},PEInterval_Eventbased{2});
display(['Ranksum test for PEI Event-based p-value : ',num2str(p_EventbasedPEI)])
set(gcf,'Position',[200 800 300 400]);xticklabels([]);yticklabels([])
xticks([1:2]);xlim([0.5 2.5])
clear mid u l

%% Figure 1e right (fish-based PEI, 1-s density bins and narrower spread for n = 8)
figure(2)
% Connect fish-level saccade and tail-flip median PEIs
PEI_saccade  = min(PEInterval_Fishbased{1}(:), 10);
PEI_tailflip = min(PEInterval_Fishbased{2}(:), 10);
nfish = numel(PEI_saccade);

plot([ones(nfish,1), 2*ones(nfish,1)]', ...
     [PEI_saccade, PEI_tailflip]', ...
     '-', 'Color', [0.75 0.75 0.75], 'LineWidth', 0.8);
hold on
for ii=1:2
    dotcolor = [0.5 0.5 0.5];
    xpos = ii;
    datavec     = PEInterval_Fishbased{ii};
    displaydata = min(datavec,10);

    limiter_vec = 0:1:10;
    limiter_vec_indx = zeros(size(displaydata));
    prop = zeros(size(limiter_vec));

    for jj=1:length(displaydata)
        [~, limiter_vec_indx(jj)] = min(abs(displaydata(jj)-limiter_vec));
    end

    for jj=1:length(limiter_vec)
        prop(jj) = sum(limiter_vec_indx == jj) ./ length(displaydata);
    end

    weight = prop(limiter_vec_indx) * 0.4;
    scatter(xpos, displaydata, 50, ...
        'MarkerFaceColor', dotcolor, ...
        'MarkerEdgeColor', 'none', ...
        'MarkerFaceAlpha', .5); hold on

    mid(ii) = nanmedian(datavec);
    u(ii)   = quantile(datavec,0.75);
    l(ii)   = quantile(datavec,0.25);

    plot([xpos-0.15 xpos+0.15],[mid(ii) mid(ii)],'k','linewidth',1.5); hold on
    plot([xpos-0.06 xpos+0.06],[u(ii) u(ii)],'k','linewidth',1.5); hold on
    plot([xpos-0.06 xpos+0.06],[l(ii) l(ii)],'k','linewidth',1.5); hold on
end

box off
p_FishbasedPEI = signrank(PEInterval_Fishbased{1},PEInterval_Fishbased{2});
display(['Signrank test for PEI Fish-based p-value : ',num2str(p_FishbasedPEI)])
set(gcf,'Position',[600 800 300 400]);xticklabels([]);yticklabels([])
xticks([1:2]); xlim([0.5 2.5]); ylim([0 10]); yticks([0:10])
clear mid u l

%% Figure 2D (event-level plots + mixed model, display capped at 10)

mid_saccade=[];u_saccade=[];l_saccade=[];
pval_saccade_clock_lme = nan(1,length(timeframe_saccade)-1);

for ii=1:length(timeframe_saccade)-1
    displaydata = saccade_clock{ii};
    displaydata_capped = min(displaydata,10);
    displaysimu = logical(ifsimu_saccade_clock{ii});
    fishidvec   = fishid_saccade_clock{ii};

    mid_saccade(ii,:) = [nanmedian(displaydata(~displaysimu)), ...
                         nanmedian(displaydata(displaysimu))];
    u_saccade(ii,:)   = min([quantile(displaydata(~displaysimu),0.75), ...
                         quantile(displaydata(displaysimu),0.75)],10);
    l_saccade(ii,:)   = [quantile(displaydata(~displaysimu),0.25), ...
                         quantile(displaydata(displaysimu),0.25)];

    tbl = table(displaydata(:), double(displaysimu(:)), categorical(fishidvec(:)), ...
                'VariableNames', {'Response','Simu','Fish'});
    lme = fitlme(tbl,'Response ~ Simu + (1|Fish)');
    coefTbl = lme.Coefficients;
    pval_saccade_clock_lme(ii) = coefTbl.pValue(strcmp(coefTbl.Name,'Simu'));

    figure(3)
    if pval_saccade_clock_lme(ii) <0.075
        plot([ii ii],[0 max(u_saccade(ii,:))],'color',[0.5 0.5 0.5],'linestyle','--','LineWidth',1); hold on
    end
end
 
plot(1:(length(timeframe_saccade)-1),mid_saccade(:,1),'color',independentColorS ,'linestyle','-','LineWidth',4); hold on
plot(1:(length(timeframe_saccade)-1),mid_saccade(:,2),'color',STColor,'linestyle','-','LineWidth',4); hold on
plot(1:(length(timeframe_saccade)-1),u_saccade(:,1),'color',independentColorS,'linestyle',':','LineWidth',2.5); hold on
plot(1:(length(timeframe_saccade)-1),u_saccade(:,2),'color',STColor,'linestyle',':','LineWidth',2.5); hold on
plot(1:(length(timeframe_saccade)-1),l_saccade(:,1),'color',independentColorS,'linestyle',':','LineWidth',2.5); hold on
plot(1:(length(timeframe_saccade)-1),l_saccade(:,2),'color',STColor,'linestyle',':','LineWidth',2.5); hold on
xticks(1:(length(timeframe_saccade)-1))
xticklabels(timeframe_saccade(1:end-1));
set(gcf,'Position',[200 300 300 400]); box off
xlim([1 7]);xticklabels([]);yticklabels([]);
ylim([2 10])

%% Figure 2E (event-level plots + mixed model, display capped at 2)

mid_tailflip=[];u_tailflip=[];l_tailflip=[];
pval_bout_clock_lme = nan(1,length(timeframe)-1);

for ii=1:length(timeframe)-1
    displaydata = bout_clock{ii};
    displaydata_capped = min(displaydata,2);
    displaysimu = logical(ifsimu_bout_clock{ii});
    fishidvec   = fishid_bout_clock{ii};

    mid_tailflip(ii,:) = [nanmedian(displaydata(~displaysimu)), ...
                          nanmedian(displaydata(displaysimu))];
    u_tailflip(ii,:)   = min([quantile(displaydata(~displaysimu),0.75), ...
                          quantile(displaydata(displaysimu),0.75)],2);
    l_tailflip(ii,:)   = [quantile(displaydata(~displaysimu),0.25), ...
                          quantile(displaydata(displaysimu),0.25)];

    tbl = table(displaydata(:), double(displaysimu(:)), categorical(fishidvec(:)), ...
                'VariableNames', {'Response','Simu','Fish'});
    lme = fitlme(tbl,'Response ~ Simu + (1|Fish)');
    coefTbl = lme.Coefficients;
    pval_bout_clock_lme(ii) = coefTbl.pValue(strcmp(coefTbl.Name,'Simu'));
end
figure(4)
plot(1:(length(timeframe)-1),smooth(mid_tailflip(:,1),3),'color',independentColorT,'linestyle','-','LineWidth',4); hold on
plot(1:(length(timeframe)-1),smooth(mid_tailflip(:,2),3),'color',STColor,'linestyle','-','LineWidth',4); hold on
plot(1:(length(timeframe)-1),mid_tailflip(:,1),'color',independentColorT,'linestyle',':','LineWidth',2.5); hold on
plot(1:(length(timeframe)-1),mid_tailflip(:,2),'color',STColor,'linestyle',':','LineWidth',2.5); hold on
plot(1:(length(timeframe)-1),u_tailflip(:,1),'color',independentColorT,'linestyle',':','LineWidth',2.5); hold on
plot(1:(length(timeframe)-1),u_tailflip(:,2),'color',STColor,'linestyle',':','LineWidth',2.5); hold on
plot(1:(length(timeframe)-1),l_tailflip(:,1),'color',independentColorT,'linestyle',':','LineWidth',2.5); hold on
plot(1:(length(timeframe)-1),l_tailflip(:,2),'color',STColor,'linestyle',':','LineWidth',2.5); hold on
xticks(1:(length(timeframe)-1));xticklabels(timeframe(1:end-1));
set(gcf,'Position',[600 300 300 400]);box off
xlim([1 21]);xticklabels([]);yticklabels([]);
ylim([0.5 2]);yticks(0.5:.1:2);

%% Pooled 10-12 lag window (Figure 2E)
pool_display_10_12 = [];
pool_simu_10_12    = [];
pool_fish_10_12    = [];

for ii_pool_10_12 = 10:12
    dispdata_pool_10_12 = bout_clock{ii_pool_10_12};
    simuflag_pool_10_12 = logical(ifsimu_bout_clock{ii_pool_10_12});
    fishid_pool_10_12   = fishid_bout_clock{ii_pool_10_12};

    pool_display_10_12 = [pool_display_10_12; dispdata_pool_10_12];
    pool_simu_10_12    = [pool_simu_10_12;    double(simuflag_pool_10_12)];
    pool_fish_10_12    = [pool_fish_10_12;    fishid_pool_10_12];
end

tbl_pool_10_12 = table(pool_display_10_12(:), pool_simu_10_12(:), categorical(pool_fish_10_12(:)), ...
    'VariableNames', {'Response','Simu','Fish'});
lme_pool_10_12 = fitlme(tbl_pool_10_12,'Response ~ Simu + (1|Fish)');
coefTbl_pool_10_12   = lme_pool_10_12.Coefficients;
pval_pool_10_12      = coefTbl_pool_10_12.pValue(strcmp(coefTbl_pool_10_12.Name,'Simu'));
disp(['Figure 2E pooled lags 10-12 (bout clock) LME p = ', num2str(pval_pool_10_12)]);

%% Pooled 17-19 lag window (Figure 2E)
pool_display_17_19 = [];
pool_simu_17_19    = [];
pool_fish_17_19    = [];
for ii_pool_17_19 = 17:19
    dispdata_pool_17_19        = bout_clock{ii_pool_17_19};
    simuflag_pool_17_19        = logical(ifsimu_bout_clock{ii_pool_17_19});
    fishid_pool_17_19          = fishid_bout_clock{ii_pool_17_19};
    pool_display_17_19 = [pool_display_17_19; dispdata_pool_17_19];
    pool_simu_17_19    = [pool_simu_17_19;    double(simuflag_pool_17_19)];
    pool_fish_17_19    = [pool_fish_17_19;    fishid_pool_17_19];
end
tbl_pool_17_19 = table(pool_display_17_19(:), pool_simu_17_19(:), categorical(pool_fish_17_19(:)), ...
    'VariableNames', {'Response','Simu','Fish'});
lme_pool_17_19 = fitlme(tbl_pool_17_19,'Response ~ Simu + (1|Fish)');
coefTbl_pool_17_19 = lme_pool_17_19.Coefficients;
pval_pool_17_19    = coefTbl_pool_17_19.pValue(strcmp(coefTbl_pool_17_19.Name,'Simu'));
disp(['Figure 2E pooled lags 17-19 (bout clock) LME p = ', num2str(pval_pool_17_19)]);

clear ii_pool_10_12 dispdata_pool_10_12 simuflag_pool_10_12 fishid_pool_10_12 ...
      tbl_pool_10_12 lme_pool_10_12 coefTbl_pool_10_12

clear ii_pool_17_19 dispdata_pool_17_19 simuflag_pool_17_19 fishid_pool_17_19 ...
      tbl_pool_17_19 lme_pool_17_19 coefTbl_pool_17_19


%% Figure 2G (event-level plots + mixed model)
figure(5)
mid_saccade = [];
pval_prop_saccade_glme = nan(1,length(timeframe_saccade)-1);
for ii=1:length(timeframe_saccade)-1
    displaydata = saccade_dir{ii};
    displaysimu = logical(ifsimu_saccade_dir{ii});
    fishidvec   = fishid_saccade_dir{ii};

    if ii~=4
        mid_saccade(ii,:) = [sum(displaydata(~displaysimu))./length(displaydata(~displaysimu)), ...
                             sum(displaydata(displaysimu))./length(displaydata(displaysimu))];
        tbl = table(double(displaydata(:)), double(displaysimu(:)), categorical(fishidvec(:)), ...
                    'VariableNames', {'Response','Simu','Fish'});
        glme = fitglme(tbl,'Response ~ Simu + (1|Fish)', ...
            'Distribution','Binomial','Link','logit');
        coefTbl = glme.Coefficients;
        pval_prop_saccade_glme(ii) = coefTbl.pValue(strcmp(coefTbl.Name,'Simu'));
    else
        mid_saccade(ii,1) = nan;
        mid_saccade(ii,2) = nan;
    end
end
plot(1:(length(timeframe_saccade)-1),mid_saccade(:,1),'color',independentColorS,'linestyle','-','LineWidth',4); hold on
plot(1:(length(timeframe_saccade)-1),mid_saccade(:,2),'color',STColor,'linestyle','-','LineWidth',4); hold on
xticks(1:(length(timeframe_saccade)-1));xlim([1 length(timeframe_saccade)-1]);
set(gcf,'Position',[1000 300 300 400])
plot([1 (length(timeframe_saccade)-1)],[.5 .5],'color',[0.5 0.5 0.5],'linewidth',2,'linestyle',':');
yticks([0.1:0.1:0.9])
xticklabels([]);yticklabels([]); ylim([0.1 0.9]);box off
clear n1 N1 n2 N2

%% Figure 2H (event-level plots + mixed model)
figure(6)
mid=[];

pval_prop_bout_glme = nan(1,length(timeframe)-1);

for ii=1:length(timeframe)-1
    if ii~=11
        displaydata = bout_dir{ii};
        displaysimu = logical(ifsimu_bout_dir{ii});
        fishidvec   = fishid_bout_dir{ii};

        mid(ii,:) = [sum(displaydata(~displaysimu))./length(displaydata(~displaysimu)) , ...
                     sum(displaydata(displaysimu))./length(displaydata(displaysimu))];

        tbl = table(double(displaydata(:)), double(displaysimu(:)), categorical(fishidvec(:)), ...
                    'VariableNames', {'Response','Simu','Fish'});
        glme = fitglme(tbl,'Response ~ Simu + (1|Fish)', ...
            'Distribution','Binomial','Link','logit');
        coefTbl = glme.Coefficients;
        pval_prop_bout_glme(ii) = coefTbl.pValue(strcmp(coefTbl.Name,'Simu'));
    else
        mid(ii,1)=nan;   mid(ii,2)=nan;
    end
end

%% Pooled 16-18 lag window (Figure 2H)
pool_resp_16_18 = [];
pool_simu_16_18 = [];
pool_fish_16_18 = [];

for ii_pool_16_18 = 16:18
    resp_pool_16_18 = bout_dir{ii_pool_16_18};
    simu_pool_16_18 = logical(ifsimu_bout_dir{ii_pool_16_18});
    fish_pool_16_18 = fishid_bout_dir{ii_pool_16_18};

    pool_resp_16_18 = [pool_resp_16_18; double(resp_pool_16_18(:))];
    pool_simu_16_18 = [pool_simu_16_18; double(simu_pool_16_18(:))];
    pool_fish_16_18 = [pool_fish_16_18; fish_pool_16_18(:)];
end

tbl_pool_16_18 = table(pool_resp_16_18(:), pool_simu_16_18(:), categorical(pool_fish_16_18(:)), ...
    'VariableNames', {'Response','Simu','Fish'});
glme_pool_16_18 = fitglme(tbl_pool_16_18, 'Response ~ Simu + (1|Fish)', ...
    'Distribution','Binomial','Link','logit');
coefTbl_pool_16_18 = glme_pool_16_18.Coefficients;
pval_pool_16_18 = coefTbl_pool_16_18.pValue(strcmp(coefTbl_pool_16_18.Name,'Simu'));
disp(['Figure 2H pooled lags 16-18 (bout dir) GLME p = ', num2str(pval_pool_16_18)]);

clear ii_pool_16_18 resp_pool_16_18 simu_pool_16_18 fish_pool_16_18 ...
      tbl_pool_16_18 glme_pool_16_18 coefTbl_pool_16_18 ...
      pool_resp_16_18 pool_simu_16_18 pool_fish_16_18

%% Pooled 4-6 lag window (Figure 6)
pool_resp_4_6 = [];
pool_simu_4_6 = [];
pool_fish_4_6 = [];

for ii_pool_4_6 = 4:6
    resp_pool_4_6 = bout_dir{ii_pool_4_6};
    simu_pool_4_6 = logical(ifsimu_bout_dir{ii_pool_4_6});
    fish_pool_4_6 = fishid_bout_dir{ii_pool_4_6};

    pool_resp_4_6 = [pool_resp_4_6; double(resp_pool_4_6(:))];
    pool_simu_4_6 = [pool_simu_4_6; double(simu_pool_4_6(:))];
    pool_fish_4_6 = [pool_fish_4_6; fish_pool_4_6(:)];
end

tbl_pool_4_6 = table(pool_resp_4_6(:), pool_simu_4_6(:), categorical(pool_fish_4_6(:)), ...
    'VariableNames', {'Response','Simu','Fish'});
glme_pool_4_6 = fitglme(tbl_pool_4_6, 'Response ~ Simu + (1|Fish)', ...
    'Distribution','Binomial','Link','logit');
coefTbl_pool_4_6 = glme_pool_4_6.Coefficients;
pval_pool_4_6 = coefTbl_pool_4_6.pValue(strcmp(coefTbl_pool_4_6.Name,'Simu'));
disp(['Figure 2H pooled lags 4-6 (bout dir)   GLME p = ', num2str(pval_pool_4_6)]);

clear ii_pool_4_6 resp_pool_4_6 simu_pool_4_6 fish_pool_4_6 ...
      tbl_pool_4_6 glme_pool_4_6 coefTbl_pool_4_6 ...
      pool_resp_4_6 pool_simu_4_6 pool_fish_4_6

plot(1:(length(timeframe)-1-11),min(max(mid(1:10,1),.35),.65),'color',independentColorT,'linestyle',':','LineWidth',2.5); hold on
plot(1:(length(timeframe)-1-11),min(max(mid(1:10,2),.35),.65),'color',STColor,'linestyle',':','LineWidth',2.5); hold on
plot(12:(length(timeframe)-1),min(max(mid(12:21,1),.35),.65),'color',independentColorT,'linestyle',':','LineWidth',2.5); hold on
plot(12:(length(timeframe)-1),min(max(mid(12:21,2),.35),.65),'color',STColor,'linestyle',':','LineWidth',2.5); hold on
plot(1:(length(timeframe)-1-11),max(smooth(mid(1:10,1),3),.35),'color',independentColorT,'linestyle','-','LineWidth',4); hold on
plot(1:(length(timeframe)-1-11),max(smooth(mid(1:10,2),3),.35),'color',STColor,'linestyle','-','LineWidth',4); hold on
plot(12:(length(timeframe)-1),max(smooth(mid(12:21,1),3),.35),'color',independentColorT,'linestyle','-','LineWidth',4); hold on
plot(12:(length(timeframe)-1),max(smooth(mid(12:21,2),3),.35),'color',STColor,'linestyle','-','LineWidth',4); hold on
xticks(1:(length(timeframe)-1));xlim([1 (length(timeframe)-1)])
set(gcf,'Position',[1400 300 300 400])
plot([1 (length(timeframe)-1)],[.5 .5],'color',[0.5 0.5 0.5],'linewidth',2,'linestyle',':');hold on
xticklabels([]);yticklabels([]);box off
ylim([0.35 0.65]); yticks([0.35:0.05:0.65])

%% Print mixed-model p-values 
disp(' ')
disp('=== Mixed-model p-values (event-level observations; fish random intercept)  ===')
disp(['Figure 2D (saccade clock) LME: ', num2str(pval_saccade_clock_lme)])
disp(['Figure 2E (bout clock)    LME: ', num2str(pval_bout_clock_lme)])
disp(['Figure 2G (saccade dir)   GLME: ', num2str(pval_prop_saccade_glme)])
disp(['Figure 2H (bout dir)       GLME: ', num2str(pval_prop_bout_glme)])