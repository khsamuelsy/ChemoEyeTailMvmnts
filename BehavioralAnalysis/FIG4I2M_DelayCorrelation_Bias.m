clear all;close all;
global gh
colorsel;
addpath('./stattool');addpath('./disptool');
datapair_sp=[];datapair_ap=[];datapair_av=[];

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
                newentry = [gh.data.simuMtx(stimsimu_n(jj),9)-gh.data.simuMtx(stimsimu_n(jj),2) , ...
                    gh.data.simuMtx(stimsimu_n(jj),11)-gh.data.simuMtx(stimsimu_n(jj),4) , ...
                    gh.data.simuMtx(stimsimu_n(jj),4)-gh.data.simuMtx(stimsimu_n(jj),9) , ...
                    -gh.data.simuMtx(stimsimu_n(jj),2)+gh.data.simuMtx(stimsimu_n(jj),4),...
                    -gh.data.simuMtx(stimsimu_n(jj),9)+gh.data.simuMtx(stimsimu_n(jj),11),...
                    abs(gh.data.simuMtx(stimsimu_n(jj),14)*-1-anglebias_overall(fishsub))];
                if group(ii)==1
                    datapair_sp = [datapair_sp; newentry];
                elseif group(ii)==2
                    datapair_ap = [datapair_ap; newentry];
                elseif group(ii)==3
                    datapair_av = [datapair_av; newentry];
                end
                clear newentry
            end
        end
    end
    clear group
end

rown = find(datapair_sp(:,1)<-25/400);datapair_sp(rown,:)=[];
[length(rown) size(datapair_sp,1)]; rown =[];
rown = find(datapair_ap(:,1)<-25/400);datapair_ap(rown,:)=[];
[length(rown) size(datapair_ap,1)]; rown =[];
rown = find(datapair_av(:,1)<-25/400);datapair_av(rown,:)=[];
[length(rown) size(datapair_av,1)]; rown =[];

% figure
for ii=1:5
    if ii<=3
        xrange{ii}=[min([datapair_sp(:,ii);datapair_ap(:,ii);datapair_av(:,ii)]),max([datapair_sp(:,ii);datapair_ap(:,ii);datapair_av(:,ii)])];
    elseif ii==4
        xrange{ii}=[0 0.625];
    elseif ii==5
        xrange{ii}=[0.1 0.3];
    end
end
xticksdef{1}=[0 .2 .4];xticksdef{2}=[-0.2 0 0.2];xticksdef{3}=[0 0.2 0.4];xticksdef{4}=[0:.125:.625];xticksdef{5}=[0.1 0.2 0.3];
drawlinewidth = 2; R_mtx =[];
for ii=1:5
    figure(ii)
    % Scatter Plot
    if ii<=3
        scatter(datapair_sp(:,ii),min(datapair_sp(:,6),20),30,'square','MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor','none','MarkerFaceAlpha',0.4); hold on
        scatter(datapair_ap(:,ii),min(datapair_ap(:,6),20),30,'diamond','MarkerFaceColor',[119 136 172]./255,'MarkerEdgeColor','none','MarkerFaceAlpha',0.4); hold on
        scatter(datapair_av(:,ii),min(datapair_av(:,6),20),30,'MarkerFaceColor',[238 129 114]./255,'MarkerEdgeColor','none','MarkerFaceAlpha',0.4); hold on
    elseif ii==4
        scatter(min(datapair_sp(:,ii)+(rand(length(datapair_sp(:,ii)),1)-0.5)*0.05,0.625),min(datapair_sp(:,6),20),30,'square','MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor','none','MarkerFaceAlpha',0.4); hold on
        scatter(min(datapair_ap(:,ii)+(rand(length(datapair_ap(:,ii)),1)-0.5)*0.05,0.625),min(datapair_ap(:,6),20),30,'diamond','MarkerFaceColor',[119 136 172]./255,'MarkerEdgeColor','none','MarkerFaceAlpha',0.4); hold on
        scatter(min(datapair_av(:,ii)+(rand(length(datapair_av(:,ii)),1)-0.5)*0.05,0.625),min(datapair_av(:,6),20),30,'MarkerFaceColor',[238 129 114]./255,'MarkerEdgeColor','none','MarkerFaceAlpha',0.4); hold on
    elseif ii==5
        scatter(max(min(datapair_sp(:,ii),0.3),0.1),min(datapair_sp(:,6),20),30,'square','MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor','none','MarkerFaceAlpha',0.4); hold on
        scatter(max(min(datapair_ap(:,ii),0.3),0.1),min(datapair_ap(:,6),20),30,'diamond','MarkerFaceColor',[119 136 172]./255,'MarkerEdgeColor','none','MarkerFaceAlpha',0.4); hold on
        scatter(max(min(datapair_av(:,ii),0.3),0.1),min(datapair_av(:,6),20),30,'MarkerFaceColor',[238 129 114]./255,'MarkerEdgeColor','none','MarkerFaceAlpha',0.4); hold on
    end

% spontaneous 
    [R,P] = corrcoef(datapair_sp(:,ii),datapair_sp(:,6));
    p = polyfit(datapair_sp(:,ii),datapair_sp(:,6),1);
    y = polyval(p,xrange{ii});
    if P(1,2)<0.05
        plot(xrange{ii},y,'linewidth',drawlinewidth,'linestyle','-','color',blankColor); hold on
    else
        plot(xrange{ii},y,'linewidth',drawlinewidth,'linestyle','--','color',blankColor); hold on
    end
    clear R P p y

    [R,P] = corrcoef(datapair_ap(:,ii),datapair_ap(:,6));
    p = polyfit(datapair_ap(:,ii),datapair_ap(:,6),1);
    y = polyval(p,xrange{ii});
    if P(1,2)<0.05
        plot(xrange{ii},y,'linewidth',drawlinewidth,'linestyle','-','color',appeColor ); hold on
    else
        plot(xrange{ii},y,'linewidth',drawlinewidth,'linestyle','--','color',appeColor); hold on
    end
    clear R P p y

    [R,P] = corrcoef(datapair_av(:,ii),datapair_av(:,6));
    p = polyfit(datapair_av(:,ii),datapair_av(:,6),1);
    y = polyval(p,xrange{ii});

    if P(1,2)<0.05
        plot(xrange{ii},y,'linewidth',drawlinewidth,'linestyle',':','color',averColor); hold on
        R_mtx(ii,:) = [P(1,2) R(1,2)];
        [p, yhat, ci] = polypredci(datapair_av(:,ii),datapair_av(:,6), 1, 0.95, xrange{ii});
        pgonci = polyshape([xrange{ii},fliplr(xrange{ii})],[yhat+ci',fliplr(yhat-ci')]);
        pgci= plot(pgonci);        pgci.FaceAlpha = .075;
        pgci.EdgeColor = 'none';        pgci.FaceColor = [238 129 114]./255;
    else
        plot(xrange{ii},y,'linewidth',drawlinewidth,'linestyle',':','color',averColor); hold on
    end
    clear R P p y
    ylim([0 20]);    xlim(xrange{ii});
    xticks(xticksdef{ii});xticklabels([]);yticklabels([]);
    set(gcf,'Position',[300*ii 800 300 300])
end