clear all;close all;colorsel;
global gh
addpath('./stattool');addpath('./disptool');
for ii=1:2
    PEInterval_Eventbased{ii} = [];
    PEInterval_Fishbased{ii}=[];
end
timeframe = [-10:11];timeframe_saccade = [-3:4];
for ii=1:length(timeframe)-1
    bout_clock{ii}=[];    ifsimu_bout_clock{ii}=[];
    bout_dir{ii}=[];    ifsimu_bout_dir{ii}=[];
end
for ii=1:length(timeframe_saccade)-1
    saccade_clock{ii}=[];    ifsimu_saccade_clock{ii}=[];
    saccade_dir{ii}=[];    ifsimu_saccade_dir{ii}=[];
end
for fishsub = 1:length(totalfishsub)
    gh.param.fishid = totalfishsub(fishsub);
    fm_behavim_main
    for ii=1:2
        PEInterval{ii} = [];
    end
    for ii=1:sessionn(fishsub)
        if ~ismember(ii,gh.param.ExcludedSession)
            % find the proportion of ipsi-contra events
            if gh.param.fishlog.trialdetails.trial(ii,1)==8
                totalbout_n = find(gh.data.boutmtx(:,1)==ii & gh.data.boutmtx(:,4)<=40);
                totalsaccade_n = find(gh.data.saccademtx(:,1)==ii & gh.data.saccademtx(:,4)<=40);
            else
                totalbout_n = find(gh.data.boutmtx(:,1)==ii & gh.data.boutmtx(:,4)<=20);
                totalsaccade_n = find(gh.data.saccademtx(:,1)==ii & gh.data.saccademtx(:,4)<=20);
            end
           
            % calculation for saccade
            for saccaden = 1:length(totalsaccade_n)
                if saccaden>1
                    PEInterval{1} = [PEInterval{1}; ...
                        (gh.data.saccademtx(totalsaccade_n(saccaden),2)-gh.data.saccademtx(totalsaccade_n(saccaden-1),2))];
                end
                % PEI and dir extraction in sequence of saccade 
                for jj=1:length(timeframe_saccade)-1
                    if saccaden>-timeframe_saccade(jj)+1 && saccaden <(length(totalsaccade_n)-timeframe_saccade(jj))+1
                        saccade_clock{jj} = [saccade_clock{jj}; (gh.data.saccademtx(totalsaccade_n(saccaden+timeframe_saccade(jj)),2)-gh.data.saccademtx(totalsaccade_n(saccaden+timeframe_saccade(jj)-1),2))];
                        ifsimu_saccade_clock{jj} = [ifsimu_saccade_clock{jj};~isempty(find((gh.data.simuMtx(:,1)==gh.data.saccademtx(totalsaccade_n(saccaden),1)) & ...
                            (gh.data.simuMtx(:,2)==gh.data.saccademtx(totalsaccade_n(saccaden),2))))];
                    end
                    if saccaden>=-timeframe_saccade(jj)+1 && saccaden <(length(totalsaccade_n)-timeframe_saccade(jj))+1
                        saccade_dir{jj} = [saccade_dir{jj}; gh.data.saccademtx(totalsaccade_n(saccaden+timeframe_saccade(jj)),7)==gh.data.saccademtx(totalsaccade_n(saccaden),7) ];
                        ifsimu_saccade_dir{jj} = [ifsimu_saccade_dir{jj};~isempty(find((gh.data.simuMtx(:,1)==gh.data.saccademtx(totalsaccade_n(saccaden),1)) & ...
                            (gh.data.simuMtx(:,2)==gh.data.saccademtx(totalsaccade_n(saccaden),2))))];
                    end
                end
            end
            
            % calculation for bout
            for boutn = 1:length(totalbout_n)
                if boutn>1
                    PEInterval{2} = [PEInterval{2}; (gh.data.bout_details(totalbout_n(boutn),2)-gh.data.bout_details(totalbout_n(boutn-1),2))];
                end
                % PEI and dir extraction in sequence of tail flip
                for jj=1:length(timeframe)-1
                    if boutn>-timeframe(jj)+1 && boutn <(length(totalbout_n)-timeframe(jj))+1
                        bout_clock{jj} = [bout_clock{jj}; (gh.data.bout_details(totalbout_n(boutn+timeframe(jj)),2)-gh.data.bout_details(totalbout_n(boutn+timeframe(jj)-1),2))];
                        ifsimu_bout_clock{jj} = [ifsimu_bout_clock{jj};~isempty(find((gh.data.simuMtx(:,1)==gh.data.bout_details(totalbout_n(boutn),1)) & ...
                            (gh.data.simuMtx(:,9)==gh.data.bout_details(totalbout_n(boutn),2))))];
                    end
                    if abs(gh.data.bout_details(totalbout_n(boutn),7)*-1-anglebias_overall(fishsub))>=1.5
                        if boutn>=-timeframe(jj)+1 && boutn <(length(totalbout_n)-timeframe(jj))+1 ...
                            && abs(gh.data.bout_details(totalbout_n(boutn+timeframe(jj)),7)*-1-anglebias_overall(fishsub))>=1.5
                            bout_dir{jj} = [bout_dir{jj}; (sign(gh.data.bout_details(totalbout_n(boutn+timeframe(jj)),7)*-1-anglebias_overall(fishsub))==sign(gh.data.bout_details(totalbout_n(boutn),7)*-1-anglebias_overall(fishsub)))];
                            ifsimu_bout_dir{jj} = [ifsimu_bout_dir{jj};~isempty(find((gh.data.simuMtx(:,1)==gh.data.bout_details(totalbout_n(boutn),1)) & ...
                                (gh.data.simuMtx(:,9)==gh.data.bout_details(totalbout_n(boutn),2))))];
                        end
                    end
                end
            end
        end
        clear totalbout_n totalsaccade_n
    end
    PEInterval_Eventbased{1} = [PEInterval_Eventbased{1};PEInterval{1}];    PEInterval_Eventbased{2} = [PEInterval_Eventbased{2};PEInterval{2}];
    PEInterval_Fishbased{1} = [PEInterval_Fishbased{1};median(PEInterval{1})];    PEInterval_Fishbased{2} = [PEInterval_Fishbased{2};median(PEInterval{2})];
end

%Figure 2c
figure(1)
for ii=1:2
    dotcolor = [0.5 0.5 0.5];
    xpos = ii;
    % setting weights from data values ranging from 0 to 10
    limiter_vec = [0:0.125:10];
    limiter_vec_indx =[];prop=[];
    tempdata=min(PEInterval_Eventbased{ii},10);
    for jj=1:length(PEInterval_Eventbased{ii})
        [val limiter_vec_indx(jj)] = min(abs(tempdata(jj)-limiter_vec));
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
    if ii==1
        scaling_factor = 0.9;
    else
        scaling_factor = 0.32;
    end
    violin(xpos,min(PEInterval_Eventbased{ii},10),'facecolor',dotcolor, ...
        'scaling',scaling_factor,'facealpha',0.2,'style',2); hold on;
    ylim([0 10])
end
p_EventbasedPEI = ranksum(PEInterval_Eventbased{1},PEInterval_Eventbased{2});
display(['Ranksum test for PEI Event-based p-value : ',num2str(p_EventbasedPEI)])
set(gcf,'Position',[200 800 300 400]);xticklabels([]);yticklabels([])
xticks([1:2]);xlim([0.5 2.5])
clear mid u l

%Figure 2d
figure(2)
for ii=1:2
    dotcolor = [0.5 0.5 0.5];
    xpos = ii;
    % setting weights from data values ranging from 0 to 10
    weight = 0.5;
    datavec = PEInterval_Fishbased{ii};
    randdispvec=[0.033 0.067 0.1];
    countrand=0;
    for jj=1:length(datavec)
        if (datavec(jj)~=max(datavec) && datavec(jj)~=min(datavec)) && ii==1
            countrand = countrand+1;
            randx = (ii-0.05+randdispvec(countrand));
        else
            randx=ii;
        end
        displaydata = datavec(jj);
        scatter(randx,min(displaydata,10),50,'MarkerFaceColor',dotcolor,'MarkerEdgeColor','none','MarkerFaceAlpha',.5); hold on
    end
    mid(ii) = nanmedian(datavec);
    u(ii) = quantile(datavec,0.75);
    l(ii) = quantile(datavec,0.25);
    plot([xpos-0.15.*1.5 xpos+0.15.*1.5],[mid(ii) mid(ii)],'k','linewidth',1.5); hold on
    plot([xpos-0.06.*1.5 xpos+0.06.*1.5],[u(ii) u(ii)],'k','linewidth',1.5); hold on
    plot([xpos-0.06.*1.5 xpos+0.06.*1.5],[l(ii) l(ii)],'k','linewidth',1.5); hold on
end
p_FishbasedPEI = ranksum(PEInterval_Fishbased{1},PEInterval_Fishbased{2});
display(['Ranksum test for PEI Fish-based p-value : ',num2str(p_FishbasedPEI)])
set(gcf,'Position',[600 800 300 400]);xticklabels([]);yticklabels([])
xticks([1:2]);xlim([0.5 2.5]);ylim([0 10]);yticks([0:10])
clear mid u l

%Figure 3e
mid_saccade=[];u_saccade=[];l_saccade=[];
pval_saccade_clock=[];
for ii=1:length(timeframe_saccade)-1
    displaydata=saccade_clock{ii};
    displaysimu = ifsimu_saccade_clock{ii};
    mid_saccade(ii,:) = [nanmedian(displaydata(find(~displaysimu))),nanmedian(displaydata(find(displaysimu)))];
    u_saccade(ii,:) = [quantile(displaydata(find(~displaysimu)),0.75),quantile(displaydata(find(displaysimu)),0.75)];
    l_saccade(ii,:) = [quantile(displaydata(find(~displaysimu)),0.25),quantile(displaydata(find(displaysimu)),0.25)];
    pval_saccade_clock(ii) = ranksum(displaydata(find(~displaysimu)),displaydata(find(displaysimu)));
    figure(3) %Figure 3M
    if pval_saccade_clock(ii) <0.075
        plot([ii ii],[0 max(u_saccade(ii,:))],'color',[0.5 0.5 0.5],'linestyle','--','LineWidth',1); hold on
    end
    clear displaydata displaysimu
end
plot([1:(length(timeframe_saccade)-1)]',mid_saccade(:,1),'color',independentColorS ,'linestyle','-','LineWidth',4); hold on
plot([1:(length(timeframe_saccade)-1)]',mid_saccade(:,2),'color',STColor,'linestyle','-','LineWidth',4); hold on
plot([1:(length(timeframe_saccade)-1)]',u_saccade(:,1),'color',independentColorS,'linestyle',':','LineWidth',2.5); hold on
plot([1:(length(timeframe_saccade)-1)]',u_saccade(:,2),'color',STColor,'linestyle',':','LineWidth',2.5); hold on
plot([1:(length(timeframe_saccade)-1)]',l_saccade(:,1),'color',independentColorS,'linestyle',':','LineWidth',2.5); hold on
plot([1:(length(timeframe_saccade)-1)]',l_saccade(:,2),'color',STColor,'linestyle',':','LineWidth',2.5); hold on
xticks([1:(length(timeframe_saccade)-1)])
xticklabels(timeframe_saccade(1:end-1));
set(gcf,'Position',[200 300 300 400]); box off
xlim([1 7]);xticklabels([]);yticklabels([]);ylim([2 10])

%Figure 3f
black = [];
orange= [];
mid_tailflip=[];u_tailflip=[];l_tailflip=[];
for ii=1:length(timeframe)-1
    displaydata = bout_clock{ii};
    displaysimu = ifsimu_bout_clock{ii};
    mid_tailflip(ii,:) = [nanmedian(displaydata(find(~displaysimu))), nanmedian(displaydata(find(displaysimu)))];
    u_tailflip(ii,:) = [quantile(displaydata(find(~displaysimu)),0.75),quantile(displaydata(find(displaysimu)),0.75)];
    l_tailflip(ii,:) = [quantile(displaydata(find(~displaysimu)),0.25),quantile(displaydata(find(displaysimu)),0.25)];
    pval_bout_clock(ii) = ranksum(displaydata(find(~displaysimu)),displaydata(find(displaysimu)));
    figure(4) 
    if ii>=10 & ii<=12
        black = [black; displaydata(find(~displaysimu))];
        orange = [orange; displaydata(find(displaysimu))];
    end
    clear displaydata displaysimu
end
% plot smoothed and raw median
plot([1:(length(timeframe)-1)]',smooth(mid_tailflip(:,1),3),'color',independentColorT,'linestyle','-','LineWidth',4); hold on
plot([1:(length(timeframe)-1)]',smooth(mid_tailflip(:,2),3),'color',STColor,'linestyle','-','LineWidth',4); hold on
plot([1:(length(timeframe)-1)]',mid_tailflip(:,1),'color',independentColorT,'linestyle',':','LineWidth',2.5); hold on
plot([1:(length(timeframe)-1)]',mid_tailflip(:,2),'color',STColor,'linestyle',':','LineWidth',2.5); hold on
% plot 75 25 quantile
plot([1:(length(timeframe)-1)]',min(u_tailflip(:,1),2),'color',independentColorT,'linestyle',':','LineWidth',2.5); hold on
plot([1:(length(timeframe)-1)]',min(u_tailflip(:,2),2),'color',STColor,'linestyle',':','LineWidth',2.5); hold on
plot([1:(length(timeframe)-1)]',l_tailflip(:,1),'color',independentColorT,'linestyle',':','LineWidth',2.5); hold on
plot([1:(length(timeframe)-1)]',l_tailflip(:,2),'color',STColor,'linestyle',':','LineWidth',2.5); hold on
xticks([1:(length(timeframe)-1)]);xticklabels(timeframe(1:end-1));
set(gcf,'Position',[600 300 300 400]);box off
xlim([1 21]);xticklabels([]);yticklabels([]);
ylim([0.5 2]);yticks([0.5:0.1:2]);

%Figure 3k
figure(5)
mid_saccade = [];
for ii=1:length(timeframe_saccade)-1
    displaydata = saccade_dir{ii};
    displaysimu = ifsimu_saccade_dir{ii};
    if ii~=4
        mid_saccade(ii,:) = [length(find(displaydata(find(~displaysimu))))./length(displaydata(find(~displaysimu))), ...
            length(find(displaydata(find(displaysimu))))./length(displaydata(find(displaysimu)))];
        
        n1 = length(find(displaydata(find(~displaysimu))));    N1 = length(displaydata(find(~displaysimu)));
        n2 = length(find(displaydata(find(displaysimu))));    N2 = length(displaydata(find(displaysimu)));
        x1 = [repmat('a',N1,1); repmat('b',N2,1)];
        x2 = [repmat(1,n1,1); repmat(2,N1-n1,1); repmat(1,n2,1); repmat(2,N2-n2,1)];
        [tbl,chi2stat,pval] = crosstab(x1,x2);
        pval_prop_saccade(ii) = pval;
    else
        mid_saccade(ii,1) = nan;
        mid_saccade(ii,2) = nan;
    end
end
plot([1:(length(timeframe_saccade)-1)]',mid_saccade(:,1),'color',independentColorS,'linestyle','-','LineWidth',4); hold on
plot([1:(length(timeframe_saccade)-1)]',mid_saccade(:,2),'color',STColor,'linestyle','-','LineWidth',4); hold on
xticks([1:(length(timeframe_saccade)-1)]);xlim([1 length(timeframe_saccade)-1]);
set(gcf,'Position',[1000 300 300 400])
plot([1 (length(timeframe)-1)],[.5 .5],'color',[0.5 0.5 0.5],'linewidth',2,'linestyle',':');
yticks([0.1:0.1:0.9])
xticklabels([]);yticklabels([]); ylim([0.1 0.9]);box off
clear n1 N1 n2 N2

figure(6)
mid=[];
cum_prop = zeros(1,4);

for ii=1:length(timeframe)-1
    if ii~=11
        displaydata = bout_dir{ii};
        displaysimu = ifsimu_bout_dir{ii};
        mid(ii,:) = [length(find(displaydata(find(~displaysimu))))./length(displaydata(find(~displaysimu))) , ...
            length(find(displaydata(find(displaysimu))))./length(displaydata(find(displaysimu)))];
        n1 = length(find(displaydata(find(~displaysimu))));        N1 = length(displaydata(find(~displaysimu)));
        n2 = length(find(displaydata(find(displaysimu))));        N2 = length(displaydata(find(displaysimu)));
        
        if ii>=16 & ii<=18
            cum_prop = cum_prop + [n1,N1,n2,N2];
        end
        x1 = [repmat('a',N1,1); repmat('b',N2,1)];
        x2 = [repmat(1,n1,1); repmat(2,N1-n1,1); repmat(1,n2,1); repmat(2,N2-n2,1)];
        [tbl,chi2stat,p_value] = crosstab(x1,x2);
        pval_prop_bout(ii) = p_value;
    else
        mid(ii,1)=nan;   mid(ii,2)=nan;
    end
end
% plot([1:(length(timeframe)-1)]',(mid(:,1)),'color',independentColorT,'linestyle',':','LineWidth',2.5); hold on
% plot([1:(length(timeframe)-1)]',(mid(:,2)),'color',STColor,'linestyle',':','LineWidth',2.5); hold on

plot([1:(length(timeframe)-1-11)]',min(max(mid(1:10,1),.35),.65),'color',independentColorT,'linestyle',':','LineWidth',2.5); hold on
plot([1:(length(timeframe)-1-11)]',min(max(mid(1:10,2),.35),.65),'color',STColor,'linestyle',':','LineWidth',2.5); hold on
plot([12:(length(timeframe)-1)]',min(max(mid(12:21,1),.35),.65),'color',independentColorT,'linestyle',':','LineWidth',2.5); hold on
plot([12:(length(timeframe)-1)]',min(max(mid(12:21,2),.35),.65),'color',STColor,'linestyle',':','LineWidth',2.5); hold on

plot([1:(length(timeframe)-1-11)]',max(smooth(mid(1:10,1),3),.35),'color',independentColorT,'linestyle','-','LineWidth',4); hold on
plot([1:(length(timeframe)-1-11)]',max(smooth(mid(1:10,2),3),.35),'color',STColor,'linestyle','-','LineWidth',4); hold on
plot([12:(length(timeframe)-1)]',max(smooth(mid(12:21,1),3),.35),'color',independentColorT,'linestyle','-','LineWidth',4); hold on
plot([12:(length(timeframe)-1)]',max(smooth(mid(12:21,2),3),.35),'color',STColor,'linestyle','-','LineWidth',4); hold on
xticks([1:(length(timeframe)-1)]);xlim([1 (length(timeframe)-1)])
set(gcf,'Position',[1400 300 300 400])
plot([1 (length(timeframe)-1)],[.5 .5],'color',[0.5 0.5 0.5],'linewidth',2,'linestyle',':');hold on
xticklabels([]);yticklabels([]);box off
ylim([0.35 0.65]);yticks([0.35:0.05:0.65])