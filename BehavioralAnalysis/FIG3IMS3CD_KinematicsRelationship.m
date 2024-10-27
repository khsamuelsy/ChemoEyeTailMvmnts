%this script plot event delays for the clueless case [figure 3b,c,d, 4h, 4i]
clear all;close all;
global gh
colorsel;
addpath('./stattool');addpath('./disptool');

bout_bias=[];bout_und_int=[];bout_clock=[];bout_mag=[];
bout_bias2=[];bout_mag2=[];

for fishsub = 1:length(totalfishsub)
    temp = [];
    gh.param.fishid = totalfishsub(fishsub);
    fm_behavim_main

    for ii=1:sessionn(fishsub)
        if ~ismember(ii,gh.param.ExcludedSession)
            % find the proportion of ipsi-contra events
            if gh.param.fishlog.trialdetails.trial(ii,1)==8
                totalbout_n = find(gh.data.boutmtx(:,1)==ii & gh.data.boutmtx(:,4)<=40);
            else
                totalbout_n = find(gh.data.boutmtx(:,1)==ii & gh.data.boutmtx(:,4)<=20);
            end
            for boutn = 1:length(totalbout_n)
                temp = [temp;gh.data.bout_details(totalbout_n(boutn),7)];
                if boutn>1
                    bout_clock = [bout_clock; (gh.data.bout_details(totalbout_n(boutn),2)-gh.data.bout_details(totalbout_n(boutn)-1,2))];
                    bout_bias = [bout_bias; abs(gh.data.bout_details(totalbout_n(boutn),7)*-1- anglebias_overall(fishsub))];
                    bout_mag = [bout_mag;abs(gh.data.bout_details(totalbout_n(boutn),8))];

                    X = gh.data.anglevect(gh.data.bout_details(totalbout_n(boutn),3):gh.data.bout_details(totalbout_n(boutn),5));
                    X0 = X-max(X);
                    X0_inv = -X0;
                    pks_X0 = findpeaks(X0);
                    pks_X0_inv = findpeaks(X0_inv);
                    und_int = gh.data.bout_details(totalbout_n(boutn),6)./((length(pks_X0)+length(pks_X0_inv))./2);
                    bout_und_int = [bout_und_int; und_int];
                end
                bout_bias2 = [bout_bias2; abs(gh.data.bout_details(totalbout_n(boutn),7)*-1- anglebias_overall(fishsub))];
                bout_mag2 = [bout_mag2;abs(gh.data.bout_details(totalbout_n(boutn),8)) ];
            end
        end
    end
    [median(temp) mean(temp)]
end

num_bootstrap=1000;
%1st Figure
X=bout_clock;Y=bout_bias;
fixBinEdges = [0,0.5:.25:2.5,max(X)];
[binCounts, binEdges, binIndices] = histcounts(X, 'BinEdges',fixBinEdges);
numBins = length(binCounts);binMedians = zeros(1, numBins);binQuan75 = zeros(1, numBins);binQuan25 = zeros(1, numBins);
for ii = 1:numBins
    binData = Y(binIndices == ii);
    if ~isempty(binData)
        binMedian(ii) = median(binData);
        bootstrap_medians = bootstrp(num_bootstrap, @median, binData);
        ci_median = prctile(bootstrap_medians, [2.5 97.5]);
        binQuan75(ii) = ci_median(2);
        binQuan25(ii) = ci_median(1);
    end
end
figure(1);
binCenters = (binEdges(1:end-1) + binEdges(2:end)) / 2;
binCenters(1) = .375; binCenters(end) = 2.625;
scatter(max(min(X,2.625),.375),min(Y,10),5,'MarkerFaceColor',[0.65 0.65 0.65],'MarkerEdgeColor','none'); hold on
for ii=1:numBins
    scatter(binCenters(ii),binMedian(ii),50,[.5 .5 .5],'MarkerFaceColor',[.5 .5 .5]);hold on
    plot([binCenters(ii) binCenters(ii)],[binQuan25(ii) binQuan75(ii)],'color',[.5 .5 .5],'LineWidth',3); hold on
end
xticks([binCenters]);yticks([0:0.5:10]);ylim([0 10]);xlim([binCenters(1) binCenters(end)]);
set(gcf,'Position',[901 900 300 400]);xticklabels({});yticklabels({});
[rho_bias,pval_bias] = corr(X, Y, 'Type', 'Spearman');
clear X Y binCounts binMedian binEdges binIndices bootstrap_medians ci_median binCenters

%2nd Figure
X=bout_clock;Y=bout_und_int;
fixBinEdges = [0,0.5:.25:2.5,max(X)];
[binCounts, binEdges, binIndices] = histcounts(X, 'BinEdges',fixBinEdges);
numBins = length(binCounts);binMedians = zeros(1, numBins);binQuan75 = zeros(1, numBins);binQuan25 = zeros(1, numBins);
for ii = 1:numBins
    binData = Y(binIndices == ii);
    if ~isempty(binData)
        binMedian(ii) = median(binData);
        bootstrap_medians = bootstrp(num_bootstrap, @median, binData);
        ci_median = prctile(bootstrap_medians, [2.5 97.5]);
        binQuan75(ii) = ci_median(2);        binQuan25(ii) = ci_median(1);
    end
end
figure(2);
binCenters = (binEdges(1:end-1) + binEdges(2:end)) / 2;
binCenters(1) = .375;
binCenters(end) = 2.625;
scatter(max(min(bout_clock,2.625),.375),max(min(Y,.044),.034),5,'MarkerFaceColor',[0.65 0.65 0.65],'MarkerEdgeColor','none'); hold on
for ii=1:numBins
    scatter(binCenters(ii),binMedian(ii),50,[.5 .5 .5],'MarkerFaceColor',[.5 .5 .5]);hold on
    plot([binCenters(ii) binCenters(ii)],[binQuan25(ii) binQuan75(ii)],'color',[.5 .5 .5],'LineWidth',3); hold on
end
xticks([binCenters]); yticks([0.034:0.001:.044]);ylim([0.034 .044]);xlim([binCenters(1) binCenters(end)])
set(gcf,'Position',[1201 900 300 400]),xticklabels({});yticklabels({});
[rho_und_int,pval_und_int] = corr(X, Y, 'Type', 'Spearman');
clear X Y binCounts binMedian binEdges binIndices bootstrap_medians ci_median binCenters

%3rd Figure
X=bout_clock; Y=bout_mag;
fixBinEdges = [0,0.5:.25:2.5,max(X)];
[binCounts, binEdges, binIndices] = histcounts(X, 'BinEdges',fixBinEdges );
numBins = length(binCounts); binMedians = zeros(1, numBins);binQuan75 = zeros(1, numBins);binQuan25 = zeros(1, numBins);
for ii = 1:numBins
    binData = Y(binIndices == ii);
    if ~isempty(binData)
        binMedian(ii) = median(binData);
        bootstrap_medians = bootstrp(num_bootstrap, @median, binData);
        ci_median = prctile(bootstrap_medians, [2.5 97.5]);
        binQuan75(ii) = ci_median(2);        binQuan25(ii) = ci_median(1);
    end
end
figure(3);
binCenters = (binEdges(1:end-1) + binEdges(2:end)) / 2;
binCenters(1) = .375; binCenters(end) = 2.625;
scatter(max(min(bout_clock,2.625),.375),max(min(Y,22),8),5,'MarkerFaceColor',[0.65 0.65 0.65],'MarkerEdgeColor','none'); hold on
for ii=1:numBins
    scatter(binCenters(ii),binMedian(ii),50,[.5 .5 .5],'MarkerFaceColor',[.5 .5 .5]);hold on
    plot([binCenters(ii) binCenters(ii)],[binQuan25(ii) binQuan75(ii)],'color',[.5 .5 .5],'LineWidth',3); hold on
end
xticks([binCenters]);yticks([8:22]);ylim([8 22]);xlim([binCenters(1) binCenters(end)]);
set(gcf,'Position',[1501 900 300 400]);xticklabels({});yticklabels({});
[rho_mag,pval_mag] = corr(X, Y, 'Type', 'Spearman');
clear X Y binCounts binMedian binEdges binIndices bootstrap_medians ci_median binCenters

%4th Figure
X=bout_mag2; Y=bout_bias2;
fixBinEdges = [0,10:4:30,max(X)];
[binCounts, binEdges, binIndices] = histcounts(X, 'BinEdges',fixBinEdges);
numBins = length(binCounts); binMedians = zeros(1, numBins); binQuan75 = zeros(1, numBins); binQuan25 = zeros(1, numBins);
for ii = 1:numBins
    binData = Y(binIndices == ii);
    if ~isempty(binData)
        binMedian(ii) = median(binData);
        bootstrap_medians = bootstrp(num_bootstrap, @median, binData);
        ci_median = prctile(bootstrap_medians, [2.5 97.5]);
        binQuan75(ii) = ci_median(2);        binQuan25(ii) = ci_median(1);
    end
end
figure(4);
binCenters = (binEdges(1:end-1) + binEdges(2:end)) / 2;
binCenters(1) = 8; binCenters(end) = 32;
scatter(max(min(X,32),8),min(Y,20),5,'MarkerFaceColor',[0.65 0.65 0.65],'MarkerEdgeColor','none'); hold on
for ii=1:numBins
    scatter(binCenters(ii),binMedian(ii),50,[.5 .5 .5],'MarkerFaceColor',[.5 .5 .5]);hold on
    plot([binCenters(ii) binCenters(ii)],[binQuan25(ii) binQuan75(ii)],'color',[.5 .5 .5],'LineWidth',3); hold on
end
xticks([binCenters]); yticks([0:2:20]); ylim([0 20]); xlim([binCenters(1) binCenters(end)]);
set(gcf,'Position',[1801 900 300 400]); xticklabels({});yticklabels({});
[rho_magbias,pval_magbias] = corr(X, Y, 'Type', 'Spearman');
clear X Y binCounts binMedian binEdges binIndices bootstrap_medians ci_median binCenters