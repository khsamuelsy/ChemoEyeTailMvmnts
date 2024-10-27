clear all; close all; 
global gh
colorsel;
addpath('./stattool');addpath('./disptool');
for ii=1:3
    datamtx{ii} = [];
end

nbin = 36;
for fishsub = 1:length(totalfishsub)
    gh.param.fishid = totalfishsub(fishsub);
    fm_behavim_main;
    for ii=1:sessionn(fishsub)
        tailn = [];          
        if ~ismember(ii,gh.param.ExcludedSession)
            for jj=1:nbin
                tailn(jj) = length(find(gh.data.bout_details(:,1)==ii&gh.data.bout_details(:,2)>10+(jj-1)*30/nbin & gh.data.bout_details(:,2)<=10+jj*30/nbin));
                simu(jj) = length(find(gh.data.simuMtx(:,1)==ii&gh.data.simuMtx(:,2)>10+(jj-1)*30/nbin & gh.data.simuMtx(:,2)<=10+jj*30/nbin));
            end
            varval = std(tailn(1:(nbin./3)));
            n1_tailn = tailn./varval;
            n2_tailn = n1_tailn-mean(n1_tailn(1:(nbin./3)));
            
            if gh.param.fishlog.trialdetails.trial(ii,1)==8
                group(ii) = 1;
            elseif gh.param.fishlog.trialdetails.trial(ii,1)==1 || gh.param.fishlog.trialdetails.trial(ii,1)==2
                group(ii) = 2;
            elseif gh.param.fishlog.trialdetails.trial(ii,1)==4 || gh.param.fishlog.trialdetails.trial(ii,1)==5
                group(ii) = 3;
            end
            if varval~=0 
                datamtx{group(ii)} = [datamtx{group(ii)};n2_tailn];
            end
            clear n1_tailn n2_tailn
        end
        clear tailn
    end
end

linecolor{1} = blankColor; linecolor{2} = appeColor; linecolor{3} = averColor;
linestyle{1} = '--'; linestyle{2} = ':'; linestyle{3} = '-';

datamtx_kw =[];

for ii=1:3
    figure(1)
    plot([1:nbin], smooth(mean(datamtx{ii}),3),'color',linecolor{ii},'linewidth',1.5,'linestyle',linestyle{ii}); hold on
    plot([1:nbin], mean(datamtx{ii}),'color',linecolor{ii},'linewidth',1,'linestyle',linestyle{ii}); hold on
    tempmtx = datamtx{ii};
    datamtx_kw(:,ii) = [...
        tempmtx(:,13);tempmtx(:,14);tempmtx(:,15);...
        tempmtx(:,16); tempmtx(:,17);tempmtx(:,18);...
        tempmtx(:,19);tempmtx(:,20);tempmtx(:,21);tempmtx(:,22);tempmtx(:,23);tempmtx(:,24);...
        nan(1000-12*length(tempmtx(:,1)),1)];
    notnancount = sum(isnan(tempmtx(:,1))==0)
    %sum(isnan(tempmtx(:,1))==1)
    pgon = polyshape([1:nbin,fliplr([1:nbin])],[nanmean(datamtx{ii})+nanstd(datamtx{ii})./sqrt(notnancount),fliplr(nanmean(datamtx{ii})-nanstd(datamtx{ii})./sqrt(notnancount))]);
    plot(pgon,'FaceColor',linecolor{ii},'FaceAlpha',0.3,'EdgeColor','none');
    tempmtx = datamtx{ii};
    xlim([1 nbin]);xticklabels([1:nbin]);    box off
    pause(1)
end

ylim([-1 1]);yticks([-1:0.1:1])
set(gcf,'Position',[800 800 100 250])
plot([.5 nbin+.5],[0 0],'k','linestyle','--')
xticklabels({});yticklabels({});yl = ylim;
pgon = polyshape([nbin./3+0.5 nbin.*2/3+0.5 nbin.*2/3+0.5 nbin./3+0.5],[yl(1) yl(1) yl(2) yl(2)]);
plot(pgon,'FaceColor','none','EdgeColor','k','linestyle','--');
xlim([.5 nbin+.5]);h = gca;h.XAxis.Visible = 'off';

[p,tbl,stats] = kruskalwallis(datamtx_kw);
c = multcompare(stats,"CType","dunn-sidak");



