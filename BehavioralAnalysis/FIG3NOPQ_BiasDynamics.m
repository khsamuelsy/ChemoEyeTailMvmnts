%this script plot event delays for the clueless case [figure 3b,c,d, 4h, 4i]
clear all;close all;
global gh
colorsel;
timeframe = [-10:11];
addpath('./stattool');addpath('./disptool');
for ii=1:length(timeframe)-1 
    bout_bias{ii}=[];   ifsimu_bout{ii}=[];    bout_und_int{ii}=[];
end

for fishsub = 1:length(totalfishsub)
    gh.param.fishid = totalfishsub(fishsub);
    fm_behavim_main

    for ii=1:sessionn(fishsub)
        if ~ismember(ii,gh.param.ExcludedSession)
            
            if gh.param.fishlog.trialdetails.trial(ii,1)==8
                totalbout_n = find(gh.data.boutmtx(:,1)==ii & gh.data.boutmtx(:,4)<=40);
            else
                totalbout_n = find(gh.data.boutmtx(:,1)==ii & gh.data.boutmtx(:,4)<=20);
            end
            
            for boutn = 1:length(totalbout_n)
 
                for jj=1:length(timeframe)-1
                    if  abs(gh.data.bout_details(totalbout_n(boutn),7)*-1- anglebias_overall(fishsub))>=1.5
                        
                        if boutn>=-timeframe(jj)+1 && boutn <(length(totalbout_n)-timeframe(jj))+1
                            %                         if abs(gh.data.bout_details(totalbout_n(boutn+timeframe(jj)),7)*-1-anglebias_overall(fishsub))<1.5
                            
                            bout_bias{jj} = [bout_bias{jj};abs(gh.data.bout_details(totalbout_n(boutn+timeframe(jj)),7)*-1- anglebias_overall(fishsub)) ];
                            ifsimu_bout{jj} = [ifsimu_bout{jj};~isempty(find((gh.data.simuMtx(:,1)==gh.data.bout_details(totalbout_n(boutn),1)) & ...
                                (gh.data.simuMtx(:,9)==gh.data.bout_details(totalbout_n(boutn),2))))];
                            
                            X = gh.data.anglevect(gh.data.bout_details(totalbout_n(boutn+timeframe(jj)),3):gh.data.bout_details(totalbout_n(boutn+timeframe(jj)),5));
                            X0 = X-max(X);
                            X0_inv = -X0;
                            pks_X0 = findpeaks(X0);
                            pks_X0_inv = findpeaks(X0_inv);
                            
                            und_int = gh.data.bout_details(totalbout_n(boutn+timeframe(jj)),6)./((length(pks_X0)+length(pks_X0_inv))./2);
                            bout_und_int{jj} = [bout_und_int{jj}; und_int];
                            %                         end
                            
                        end
                    end
                end
            end
        end
    end
end


black = [];
orange= [];

mid=[];u=[];l=[];
for ii=1:length(timeframe)-1
    
    displaydata = bout_bias{ii};
    displaysimu = ifsimu_bout{ii};
    %6 11 16%
    
    mid(ii,:) = [nanmedian(displaydata(find(~displaysimu))), nanmedian(displaydata(find(displaysimu)))];
    u(ii,:) = [quantile(displaydata(find(~displaysimu)),0.75), quantile(displaydata(find(displaysimu)),0.75)];
    l(ii,:) = [quantile(displaydata(find(~displaysimu)),0.25), quantile(displaydata(find(displaysimu)),0.25)];
    p_val_boutdynamics(1,ii) = ranksum(displaydata(find(~displaysimu)),displaydata(find(displaysimu)));
    
    if ii>=3 & ii<=5
        black = [black; displaydata(find(~displaysimu))];
        orange = [orange; displaydata(find(displaysimu))];
    end
end
figure(1)
ymax=8;

plot([1:(length(timeframe)-1)]',mid(:,1),'color',independentColorT,'linestyle',':','LineWidth',2.5); hold on
plot([1:(length(timeframe)-1)]',mid(:,2),'color',STColor,'linestyle',':','LineWidth',2.5); hold on
plot([1:(length(timeframe)-1)]',smooth(mid(:,1),3),'color',independentColorT,'linestyle','-','LineWidth',4); hold on
plot([1:(length(timeframe)-1)]',smooth(mid(:,2),3),'color',STColor,'linestyle','-','LineWidth',4); hold on
plot([1:(length(timeframe)-1)]',min(u(:,1),ymax),'color',independentColorT,'linestyle',':','LineWidth',2.5); hold on
plot([1:(length(timeframe)-1)]',min(u(:,2),ymax),'color',STColor,'linestyle',':','LineWidth',2.5); hold on
plot([1:(length(timeframe)-1)]',min(l(:,1),ymax),'color',independentColorT,'linestyle',':','LineWidth',2.5); hold on
plot([1:(length(timeframe)-1)]',min(l(:,2),ymax),'color',STColor,'linestyle',':','LineWidth',2.5); hold on

xticks([1:(length(timeframe)-1)]);xlim([1,(length(timeframe)-1)])
xticklabels([]);yticklabels([]);box off
ylim([1 8]);yticks([1:0.5:8])
set(gcf,'Position',[300 300 300 400])



figure(2)
mid2=[];u2=[];l2=[];
% black = [];
% orange= [];
for ii=1:length(timeframe)-1
    displaydata = bout_und_int{ii};
    displaysimu = ifsimu_bout{ii};

    mid2(ii,:) = [nanmedian(displaydata(find(~displaysimu))), nanmedian(displaydata(find(displaysimu)))];
    u2(ii,:) = [quantile(displaydata(find(~displaysimu)),0.75), quantile(displaydata(find(displaysimu)),0.75)];
    l2(ii,:) =[quantile(displaydata(find(~displaysimu)),0.25), quantile(displaydata(find(displaysimu)),0.25)];
    
    p_val_boutdynamics(2,ii) = ranksum(displaydata(find(~displaysimu)),displaydata(find(displaysimu)));
    % if ii>=18 & ii<=20
    %     black = [black; displaydata(find(~displaysimu))];
    %     orange = [orange; displaydata(find(displaysimu))];
    % end
end
plot([1:(length(timeframe)-1)]',mid2(:,1),'color',independentColorT,'linestyle',':','LineWidth',2.5); hold on
plot([1:(length(timeframe)-1)]',mid2(:,2),'color',STColor,'linestyle',':','LineWidth',2.5); hold on
plot([1:(length(timeframe)-1)]',smooth(mid2(:,1),3),'color',independentColorT,'linestyle','-','LineWidth',4); hold on
plot([1:(length(timeframe)-1)]',smooth(mid2(:,2),3),'color',STColor,'linestyle','-','LineWidth',4); hold on
plot([1:(length(timeframe)-1)]',max(min(u2(:,1),.046),.034),'color',independentColorT,'linestyle',':','LineWidth',2.5); hold on
plot([1:(length(timeframe)-1)]',max(min(u2(:,2),.046),.034),'color',STColor,'linestyle',':','LineWidth',2.5); hold on
plot([1:(length(timeframe)-1)]',max(min(l2(:,1),.046),.034),'color',independentColorT,'linestyle',':','LineWidth',2.5); hold on
plot([1:(length(timeframe)-1)]',max(min(l2(:,2),.046),.034),'color',STColor,'linestyle',':','LineWidth',2.5); hold on
xticks([1:(length(timeframe)-1)]);xlim([1,(length(timeframe)-1)])
xticklabels([]);yticklabels([]);box off
ylim([.034 .046]); yticks([0.034:0.001:0.046]);
set(gcf,'Position',[600 300 300 400])


