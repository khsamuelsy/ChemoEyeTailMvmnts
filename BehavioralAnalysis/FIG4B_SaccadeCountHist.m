clear all;close all;

global gh
colorsel;
addpath('./stattool');addpath('./disptool');
for ii=1:5
    datavec{ii} = [];    
    datavec_fishbased{ii} = [];
end
for fishsub = 1:length(totalfishsub)
    gh.param.fishid = totalfishsub(fishsub);
    fm_behavim_main
    for jj=1:5
        groupdata{jj}=[];
    end
    for ii=1:sessionn(fishsub)
        if ~ismember(ii,gh.param.ExcludedSession)
            stim_saccaden = find(gh.data.saccademtx(:,1)==ii & gh.data.saccademtx(:,2)>20 & gh.data.saccademtx(:,4)<=30);
            prestim_saccaden = find(gh.data.saccademtx(:,1)==ii & gh.data.saccademtx(:,2)>10 & gh.data.saccademtx(:,4)<=20);
            
            if gh.param.fishlog.trialdetails.trial(ii,1)==8
                group(ii) = 1;
            elseif gh.param.fishlog.trialdetails.trial(ii,1)==1
                group(ii) = 2;
            elseif gh.param.fishlog.trialdetails.trial(ii,1)==2
                group(ii) = 3;
            elseif gh.param.fishlog.trialdetails.trial(ii,1)==4
                group(ii) = 4;
            elseif gh.param.fishlog.trialdetails.trial(ii,1)==5
                group(ii) = 5;
            end

            datacmp = length(stim_saccaden)-length(prestim_saccaden);
            datavec{group(ii)} = [datavec{group(ii)};datacmp];
            groupdata{group(ii)} = [groupdata{group(ii)};datacmp];
        end
    end
    for jj=1:5
        datavec_fishbased{jj} = [datavec_fishbased{jj}; mean(groupdata{jj})];
    end
    clear datacmp groupdata
end

specifiedcolor{1}=blankColor;
specifiedcolor{2}=appeColor;specifiedcolor{3}=appeColor;
specifiedcolor{4}=averColor;specifiedcolor{5}=averColor;
for ii=1:5
    figure(ii)
    plot([0 0],[0 0.8],'color',[0.5 0.5 0.5],'linestyle',':','linewidth',2); hold on;
    histogram(min(datavec{ii},3.5),[-3.5:1:3.5],'normalization','probability','displaystyle','stairs','LineWidth',4,'EdgeColor',specifiedcolor{ii}); hold on;
    set(gcf,'Position',[300 1400-ii*220 450 120]); ylim([0 0.8]);
    scatter(sort(datavec_fishbased{ii}),[0.1:0.1:0.5],60,'MarkerFaceColor',specifiedcolor{ii},'MarkerEdgeColor','none','MarkerFaceAlpha',1); hold on;
    xticklabels({});yticklabels([]);box off;yticks([0.:0.2:0.8]);xlim([-3.5 3.5])
end
    
datamtx=[];
datamtx(:,1) = [datavec{1};nan(500-length(datavec{1}),1)];
datamtx(:,2) = [datavec{2};datavec{3};nan(500-length(datavec{2})-length(datavec{3}),1)];
datamtx(:,3) = [datavec{4};datavec{5};nan(500-length(datavec{4})-length(datavec{5}),1)];
[p,tbl,stats] = kruskalwallis(datamtx);
c = multcompare(stats,"CType","dunn-sidak");