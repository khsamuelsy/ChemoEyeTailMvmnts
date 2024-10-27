clear all;close all;
global gh
colorsel; addpath('./stattool'); addpath('./disptool');
datapair_sp=[]; datapair_ap=[]; datapair_av=[];

for fishsub = 1:length(totalfishsub)
    gh.param.fishid = totalfishsub(fishsub);
    fm_behavim_main
    for ii=1:sessionn(fishsub)
        stimsimu_n=[];
        if ~ismember(ii,gh.param.ExcludedSession)
            stimsimu_n = find(gh.data.simuMtx(:,1)==ii & gh.data.simuMtx(:,2)>20 & gh.data.simuMtx(:,4)<=30);
            if gh.param.fishlog.trialdetails.trial(ii,1)==8
                group(ii) = 1;
            elseif gh.param.fishlog.trialdetails.trial(ii,1)==1 || gh.param.fishlog.trialdetails.trial(ii,1)==2
                group(ii) = 2;
            elseif gh.param.fishlog.trialdetails.trial(ii,1)==4 || gh.param.fishlog.trialdetails.trial(ii,1)==5
                group(ii) = 3;
            end
            for jj=1:length(stimsimu_n)
                if group(ii)==1 & gh.data.simuMtx(stimsimu_n(jj),9)-gh.data.simuMtx(stimsimu_n(jj),2)>=-0.0625
                    datapair_sp = [datapair_sp;   ...
                        gh.data.simuMtx(stimsimu_n(jj),9)-20, ...
                        abs(gh.data.simuMtx(stimsimu_n(jj),14)*-1-anglebias_overall(fishsub))];
                elseif group(ii)==2 & gh.data.simuMtx(stimsimu_n(jj),9)-gh.data.simuMtx(stimsimu_n(jj),2)>=-0.0625
                    datapair_ap = [datapair_ap;   ...
                        gh.data.simuMtx(stimsimu_n(jj),9)-20, ...
                        abs(gh.data.simuMtx(stimsimu_n(jj),14)*-1-anglebias_overall(fishsub))];
                elseif group(ii)==3 & gh.data.simuMtx(stimsimu_n(jj),9)-gh.data.simuMtx(stimsimu_n(jj),2)>=-0.0625
                    datapair_av = [datapair_av;   ...
                        gh.data.simuMtx(stimsimu_n(jj),9)-20, ...
                        abs(gh.data.simuMtx(stimsimu_n(jj),14)*-1-anglebias_overall(fishsub))];
                end
            end
        end
    end
end

figure
scatter(datapair_sp(:,1),min(datapair_sp(:,2),20),30,'square','MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor','none','MarkerFaceAlpha',0.4); hold on
scatter(datapair_ap(:,1),min(datapair_ap(:,2),20),30,'diamond','MarkerFaceColor',[119 136 172]./255,'MarkerEdgeColor','none','MarkerFaceAlpha',0.4); hold on
scatter(datapair_av(:,1),min(datapair_av(:,2),20),30,'MarkerFaceColor',[238 129 114]./255,'MarkerEdgeColor','none','MarkerFaceAlpha',0.4); hold on
xrange=[0 10];drawlinewidth = 2;

[R,P] = corrcoef(datapair_sp(:,1),datapair_sp(:,2));
p = polyfit(datapair_sp(:,1),datapair_sp(:,2),1);
y = polyval(p,xrange); plot(xrange,y,'linewidth',drawlinewidth,'linestyle','--','color',blankColor); hold on
clear R P p y

[R,P] = corrcoef(datapair_ap(:,1),datapair_ap(:,2));
p = polyfit(datapair_ap(:,1),datapair_ap(:,2),1);
y = polyval(p,xrange); plot(xrange,y,'linewidth',drawlinewidth,'linestyle','--','color',appeColor); hold on
clear R P p y

[R,P] = corrcoef(datapair_av(:,1),datapair_av(:,2));
p = polyfit(datapair_av(:,1),datapair_av(:,2),1);
y = polyval(p,xrange); plot(xrange,y,'linewidth',drawlinewidth,'linestyle',':','color',averColor); hold on
clear R P p y

[ p, yhat, ci] = polypredci(datapair_av(:,1),datapair_av(:,2), 1, 0.95, xrange);
pgonci = polyshape([xrange,fliplr(xrange)],[yhat+ci',fliplr(yhat-ci')]);
pgci= plot(pgonci);
pgci.FaceAlpha = .075;
pgci.EdgeColor = 'none';
pgci.FaceColor = averColor;
ylim([0 20]);yticks([0:2:20]);xticks([0:10]);xlim(xrange);
xticklabels([]);yticklabels([]);
set(gcf,'Position',[300 800 300 300])