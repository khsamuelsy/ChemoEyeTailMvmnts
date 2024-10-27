clear all; close all;
addpath('./stattool');addpath('./disptool');
load('previous_StaticArenaData.mat');colorsel;

% figure(1)
% cal.all_PBI =cal.all_PBI./240;
% cal.all_PBI2 =cal.all_PBI2./240;
% plot(smooth(nanmedian(cal.all_PBI),3),'color',independentColorT,'linestyle','-','LineWidth',4); hold on
% plot(smooth(nanmedian(cal.all_PBI2),3),'color',STColor,'linestyle','-','LineWidth',4); hold on; %large angle
% plot(nanmedian(cal.all_PBI),'color',independentColorT,'linestyle',':','LineWidth',2.5); hold on;
% plot(nanmedian(cal.all_PBI2),'color',STColor,'linestyle',':','LineWidth',2.5); hold on;
% black=[];
% orange=[];
% for ii=1:21
%         if ii >=10 & ii<=12
%         black=[black;cal.all_PBI(:,ii)];
%         orange=[orange;cal.all_PBI2(:,ii)];
%         end
%     u(ii,:) = [quantile(cal.all_PBI(:,ii),0.75) quantile(cal.all_PBI2(:,ii),0.75)];
%     l(ii,:) = [quantile(cal.all_PBI(:,ii),0.25) quantile(cal.all_PBI2(:,ii),0.25)];
% end
% plot(min(u(:,1),2),'color',independentColorT,'linestyle',':','LineWidth',2.5); hold on
% plot(min(u(:,2),2),'color',STColor,'linestyle',':','LineWidth',2.5); hold on
% plot(l(:,1),'color',independentColorT,'linestyle',':','LineWidth',2.5); hold on
% plot(l(:,2),'color',STColor,'linestyle',':','LineWidth',2.5);  hold on
% set(gcf,'Position',[100 900 300 400]);box off
% xlim([1 21]);ylim([0.5 2])
% xticks([1:21]);yticks([0.5:0.1:2])
% xticklabels([]);yticklabels([])
% 


% n1 = 0;n2 = 0;N1 = 0;N2 = 0;
% for ii=1:21
% 
%     a1=sum(cal.all_DirSignHist(:,ii)==1);
%     b1=sum(cal.all_DirSignHist(:,ii)==0);
%     c1(ii) = a1./(a1+b1);
%     a2=sum(cal.all_DirSignHist2(:,ii)==1);
%     b2=sum(cal.all_DirSignHist2(:,ii)==0);
%     c2(ii) = a2./(a2+b2);
%     if ii>=7 & ii<=9
%         n1 = n1 + a1;
%         N1 = N1 + a1+b1;
%         n2 = n2 + a2;
%         N2 = N2 + a2+b2;
%     end
% %     if ii==10
% %         w1 = a1+a2;
% %         w2 = a1+a2+b1+b2;
% %         b1+b2
% %     end
% end
% 
% 
% figure(2)
% plot((c1(1:10)),'color',independentColorT,'linestyle',':','LineWidth',2.5);hold on;
% plot(12:21,(c1(12:end)),'color',independentColorT,'linestyle',':','LineWidth',2.5);hold on
% plot(smooth(c1(1:10),3),'color',independentColorT,'linestyle','-','LineWidth',4);hold on;
% plot(12:21,smooth(c1(12:end),3),'color',independentColorT,'linestyle','-','LineWidth',4);hold on
% plot((c2(1:10)),'color',STColor,'linestyle',':','LineWidth',2.5);hold on;
% plot(12:21,(c2(12:end)),'color',STColor,'linestyle',':','LineWidth',2.5);hold on
% plot(smooth(c2(1:10),3),'color',STColor,'linestyle','-','LineWidth',4);hold on;
% plot(12:21,smooth(c2(12:end),3),'color',STColor,'linestyle','-','LineWidth',4);hold on
% plot([1 21],[.5 .5],'color',[0.5 0.5 0.5],'linewidth',2,'linestyle',':');hold on
% xlim([1 21]);ylim([0.35 .65])
% xticks([1:21]);yticks([0.35:0.05:0.65])
% set(gcf,'Position',[450 900 300 400]);box off
% xticklabels([]);yticklabels([])


% figure(3)
% 
% plot(min(nanmean(cal.all_DirHist),80),'color',independentColorT,'linestyle',':','LineWidth',2.5); hold on
% plot(min(nanmean(cal.all_DirHist2),80),'color',STColor,'linestyle',':','LineWidth',2.5); hold on;
% plot(min(smooth(nanmean(cal.all_DirHist),3),80),'color',independentColorT,'linestyle','-','LineWidth',4); hold on
% plot(min(smooth(nanmean(cal.all_DirHist2),3),80),'color',STColor,'linestyle','-','LineWidth',4); hold on;
% black=[];
% orange=[];
% for ii=1:21
%     if ii >=19 & ii<=21
%         black=[black;cal.all_DirHist(:,ii)];
%         orange=[orange;cal.all_DirHist2(:,ii)];
%     end
%     sem1 = nanstd(cal.all_DirHist(:,ii))./sqrt(sum(~isnan(cal.all_DirHist(:,ii))));
%     sem2 = nanstd(cal.all_DirHist2(:,ii))./sqrt(sum(~isnan(cal.all_DirHist2(:,ii))));
%     u(ii,:) = [nanmean(cal.all_DirHist(:,ii))+sem1, nanmean(cal.all_DirHist2(:,ii))+sem2];
%     l(ii,:) = [nanmean(cal.all_DirHist(:,ii))-sem1, nanmean(cal.all_DirHist2(:,ii))-sem2];
% end
% plot(min(u(:,1),80),'color',independentColorT,'linestyle',':','LineWidth',2.5); hold on
% plot(min(u(:,2),80),'color',STColor,'linestyle',':','LineWidth',2.5); hold on
% plot(l(:,1),'color',independentColorT,'linestyle',':','LineWidth',2.5); hold on
% plot(l(:,2),'color',STColor,'linestyle',':','LineWidth',2.5);  hold on
% set(gcf,'Position',[800 900 300 400]);box off
% xlim([1 21]);
% xticks([1:21]);
% yticks([0:10:80])
% ylim([0 80])
% xticklabels([]);yticklabels([])

figure(4)
ax = gca();
pieData = [305 346];
% pieData = [441 501];
h = pie(ax, pieData);
% Define 3 colors, one for each of the 3 wedges
newColors = [...
    .65,       0.65, 0.65;   %hot pink
    0.9,       .9,       .9];  %dark orchid
ax.Colormap = newColors;
set(gcf,'Position',[1150 900 300 300])


% num_bootstrap=1000;
% %1st Figure
% X1 = cal.all_datamtx(:,1)./240;
% fixBinEdges = [0,0.5:.25:2.5,max(X1)];
% [binCounts, binEdges, binIndices] = histcounts(X1, 'BinEdges',fixBinEdges );
% numBins = length(binCounts);
% binMeans = zeros(1, numBins );
% binVariances = zeros(1, numBins );
% y1 = cal.all_datamtx(:,3);
% for ii = 1:numBins
%     binData = y1(binIndices == ii);
%     if ~isempty(binData)
%         binMedian(ii) = median(binData);
%         bootstrap_medians = bootstrp(num_bootstrap, @median, binData);
%         ci_median = prctile(bootstrap_medians, [2.5 97.5]);
%         binQuan75(ii) = ci_median(2);
%         binQuan25(ii) = ci_median(1);
%     end
% end
% figure(5);
% binCenters = (binEdges(1:end-1) + binEdges(2:end)) / 2;
% binCenters(1) = .375;
% binCenters(end) = 2.625;
% scatter(max(min(X1,2.625),.375),min(y1+(rand(size(y1))-.5)*3,80.5),2,'MarkerFaceColor',[0.65 0.65 0.65],'MarkerEdgeColor','none'); hold on
% for ii=1:numBins
%     scatter(binCenters(ii),binMedian(ii),50,[.5 .5 .5],'MarkerFaceColor',[.5 .5 .5]);hold on
%     plot([binCenters(ii) binCenters(ii)],[binQuan25(ii) binQuan75(ii)],'color',[.5 .5 .5],'LineWidth',3); hold on
% end
% xticks([binCenters]);yticks([0:10:80]);ylim([-1.5 81.5])
% xlim([binCenters(1) binCenters(end)]);
% set(gcf,'Position',[1501 900 300 400])
% xticklabels({});yticklabels({});
% [rho_bias,pval_bias] = corr(X1, y1, 'Type', 'Spearman');