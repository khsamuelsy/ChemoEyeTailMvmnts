clear all; close all; global gh
colorsel;
addpath('./stattool');addpath('./disptool');
for ii=1:3
    datamtx{ii} = [];
end

nbin = 15;
for fishsub = 1:length(totalfishsub)
    gh.param.fishid = totalfishsub(fishsub);
    fm_behavim_main;
    for ii=1:sessionn(fishsub)
        tailn = [];  simu = []; saccaden = []; ration = [];
        if ~ismember(ii,gh.param.ExcludedSession)
            for jj=1:nbin
                tailn(jj) = length(find(gh.data.bout_details(:,1)==ii & gh.data.bout_details(:,2)>(10+(jj-1)*30/nbin) & gh.data.bout_details(:,2)<=(10+jj*30/nbin)));
                simu(jj) = length(find(gh.data.simuMtx(:,1)==ii & gh.data.simuMtx(:,9)>(10+(jj-1)*30/nbin) & gh.data.simuMtx(:,9)<=(10+jj*30/nbin)));
                saccaden(jj) = length(find(gh.data.saccademtx(:,1)==ii & gh.data.saccademtx(:,2)>(10+(jj-1)*30/nbin) & gh.data.saccademtx(:,2)<=(10+jj*30/nbin)));
                ration(jj) = simu(jj)-(tailn(jj)-simu(jj));
            end
            
            n_ration = ration-mean(ration(1:(nbin./3)));
            
            if gh.param.fishlog.trialdetails.trial(ii,1)==8
                group(ii) = 1;
            elseif gh.param.fishlog.trialdetails.trial(ii,1)==1 || gh.param.fishlog.trialdetails.trial(ii,1)==2
                group(ii) = 2;
            elseif gh.param.fishlog.trialdetails.trial(ii,1)==4 || gh.param.fishlog.trialdetails.trial(ii,1)==5
                group(ii) = 3;
            end
            
            trialboutn = find(gh.data.bout_details(:,1)==ii);
            
            datamtx{group(ii)} = [datamtx{group(ii)};n_ration/2];
            clear n_ration
        end
        clear tailn simu saccaden ration
    end
end

linecolor{1} = blankColor; linecolor{2} = appeColor ; linecolor{3} = averColor;
linestyle{1} = '--'; linestyle{2} = ':'; linestyle{3} = '-';

datamtx_kw =[];

for ii=1:3
    figure(1)
    plot([1:nbin], nanmean(datamtx{ii}),'color',linecolor{ii},'linewidth',1.5,'linestyle',linestyle{ii}); hold on
    tempmtx = datamtx{ii};
    size(tempmtx)
    datamtx_kw(:,ii) = [tempmtx(:,6);tempmtx(:,7);tempmtx(:,8);tempmtx(:,9);tempmtx(:,10);nan(500-5*length(tempmtx(:,1)),1)];
    notnancount = sum(isnan(tempmtx(:,1))==0);
    pgon = polyshape([1:nbin,fliplr([1:nbin])],[nanmean(datamtx{ii})+nanstd(datamtx{ii})./sqrt(notnancount),fliplr(nanmean(datamtx{ii})-nanstd(datamtx{ii})./sqrt(notnancount))]);
    plot(pgon,'FaceColor',linecolor{ii},'FaceAlpha',0.3,'EdgeColor','none');
    tempmtx = datamtx{ii};
    xlim([1 nbin])
    xticklabels([1:nbin])
    box off
    pause(1)
end

ylim([-.3 .3]);yticks([-.3:0.05:.3])
set(gcf,'Position',[800 800 100 250])
plot([.5 nbin+.5],[0 0],'k','linestyle','--')
xticklabels({});yticklabels({});yl = ylim;
pgon = polyshape([nbin./3+0.5 nbin.*2/3+0.5 nbin.*2/3+0.5 nbin./3+0.5],[yl(1) yl(1) yl(2) yl(2)]);
plot(pgon,'FaceColor','none','EdgeColor','k','linestyle','--');
xlim([.5 nbin+.5]);h = gca;h.XAxis.Visible = 'off';
[p,tbl,stats] = kruskalwallis(datamtx_kw);
c = multcompare(stats,"CType","dunn-sidak");

ylim([-.3 0.3]);

