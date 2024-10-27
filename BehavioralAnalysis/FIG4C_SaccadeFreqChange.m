clear all; close all; global gh
colorsel;
addpath('./stattool');addpath('./disptool');

for ii=1:3
    datamtx{ii} = [];
end

for fishsub = 1:length(totalfishsub)
    gh.param.fishid = totalfishsub(fishsub);

    fm_behavim_main;
    for ii=1:sessionn(fishsub)
        saccaden=[];
        if ~ismember(ii,gh.param.ExcludedSession)
            for jj=1:6
                saccaden(jj) = length(find(gh.data.saccademtx(:,1)==ii&gh.data.saccademtx(:,2)>(10+(jj-1)*30/6) & gh.data.saccademtx(:,2)<=(10+jj*30/6)));
            end
            n_saccaden = saccaden-mean(saccaden(1:2));
            
            if gh.param.fishlog.trialdetails.trial(ii,1)==8
                group(ii) = 1;
            elseif gh.param.fishlog.trialdetails.trial(ii,1)==1 || gh.param.fishlog.trialdetails.trial(ii,1)==2
                group(ii) = 2;
            elseif gh.param.fishlog.trialdetails.trial(ii,1)==4 || gh.param.fishlog.trialdetails.trial(ii,1)==5
                group(ii) = 3;
            end
            datamtx{group(ii)} = [datamtx{group(ii)};n_saccaden./5];
            clear n_saccaden
        end
        clear saccaden
    end
end

linecolor{1} = blankColor; linecolor{2} = appeColor; linecolor{3} = averColor;
linestyle{1} = '--'; linestyle{2} = ':'; linestyle{3} = '-';

datamtx_kw =[];

for ii=1:3
    figure(1)
    plot([1:6], mean(datamtx{ii}),'color',linecolor{ii},'linewidth',1.5,'linestyle',linestyle{ii}); hold on
    tempmtx = datamtx{ii};
    datamtx_kw(:,ii) = [tempmtx(:,3);tempmtx(:,4);nan(500-2*length(tempmtx(:,1)),1)];
    notnancount = sum(isnan(tempmtx(:,1))==0);
    pgon = polyshape([1:6,fliplr([1:6])],[nanmean(datamtx{ii})+(nanstd(datamtx{ii})./sqrt(notnancount)),fliplr(nanmean(datamtx{ii})-(nanstd(datamtx{ii})./sqrt(notnancount)))]);
    plot(pgon,'FaceColor',linecolor{ii},'FaceAlpha',0.3,'EdgeColor','none');
    tempmtx = datamtx{ii};
    xlim([1 6]);xticklabels([1:6]);
    box off;
end


ylim([-0.08 0.08]);yticks([-0.08:0.01:0.08])
set(gcf,'Position',[800 800 100 250])
plot([.5 6.5],[0 0],'k','linestyle','--')
xticklabels({});yticklabels({});yl = ylim;
pgon = polyshape([2.5 4.5 4.5 2.5],[yl(1) yl(1) yl(2) yl(2)]);
plot(pgon,'FaceColor','none','EdgeColor','k','linestyle','--');
xlim([.5 6.5])
h = gca;h.XAxis.Visible = 'off';
[p,tbl,stats] = kruskalwallis(datamtx_kw);
c = multcompare(stats,"CType","dunn-sidak");

