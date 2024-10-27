clear all;close all;

global gh
colorsel;
addpath('./stattool');addpath('./disptool');
for ii=1:5
    datavec{ii} = [];
end
for fishsub = 1:length(totalfishsub)
    gh.param.fishid = totalfishsub(fishsub);
    
    fm_behavim_main
   
    for ii=1:sessionn(fishsub)
        if ~ismember(ii,gh.param.ExcludedSession)
            stim_boutn = find(gh.data.boutmtx(:,1)==ii&gh.data.boutmtx(:,2)>20 & gh.data.boutmtx(:,4)<=30);
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
            for jj=2:length(stim_boutn)
                datavec{group(ii)} = [datavec{group(ii)}; (gh.data.bout_details(stim_boutn(jj),2)-gh.data.bout_details(stim_boutn(jj-1),2))];
            end
            
        end
    end
end

for ii=1:5
    if ii==2
        xpos = 2.2;
    elseif ii==3
        xpos = 2.8;
    elseif ii==4
        xpos = 3.7; %3.8 for multiodor cmp
    elseif ii==5
        xpos = 4.3;
    else
        xpos = ii;
    end
    figure(1)
    dotcolor{1} = blankColor;
    dotcolor{2} = appeColor ;    dotcolor{3} = appeColor ;
    dotcolor{4} = averColor;    dotcolor{5} = averColor;
    prop=[];
    limiter_vec = [0.5:0.05:1.5];
    limiter_vec_indx =[];prop=[];
    tempdata=max(min(datavec{ii},2),0.5);
    for jj=1:length(datavec{ii})
        [val limiter_vec_indx(jj)] = min(abs(tempdata(jj)-limiter_vec));
    end
    for jj=1:length(limiter_vec)
        prop(jj)=length(find(limiter_vec_indx==jj))./length(datavec{ii});
    end
    
    weight = prop(limiter_vec_indx).*2;
    randx = xpos+rand(1,length(datavec{ii})).*weight-weight./2;
    scatter(randx,max(min(datavec{ii},1.5),0.5),20,'MarkerFaceColor',dotcolor{ii},'MarkerEdgeColor','none','MarkerFaceAlpha',.1); hold on
    mid(ii) = nanmedian(datavec{ii});
    u(ii) = quantile(datavec{ii},0.75);
    l(ii) = quantile(datavec{ii},0.25);

    plot([xpos-0.15 xpos+0.15],[mid(ii) mid(ii)],'color', dotcolor{ii},'linewidth',1.5); hold on
    plot([xpos-0.06 xpos+0.06],[u(ii) u(ii)],'color',dotcolor{ii},'linewidth',1.5); hold on
    plot([xpos-0.06 xpos+0.06],[l(ii) l(ii)],'color',dotcolor{ii},'linewidth',1.5); hold on
    violin(xpos,max(min(datavec{ii},2),0.5),'facecolor',dotcolor{ii}, ...
        'scaling',0.13,'facealpha',0.2,'style',2); hold on;
    set(gcf,'Position',[800 800 300 400])
    weight = [];limiter_vec_indx=[];
end
yticklabels([]);xticklabels({})
h = gca; h.XAxis.Visible = 'off';
xticks([1 2.5 4]); xlim([0 5]);ylim([.5 1.5]);
yticks([0.5:0.1:2]);
datamtx=[];
datamtx(:,1) = [datavec{1};nan(5000-length(datavec{1}),1)];
datamtx(:,2) = [datavec{2};datavec{3};nan(5000-length(datavec{2})-length(datavec{3}),1)];
datamtx(:,3) = [datavec{4};datavec{5};nan(5000-length(datavec{4})-length(datavec{5}),1)];
[p,tbl,stats] = kruskalwallis(datamtx);
c = multcompare(stats,"CType","dunn-sidak");