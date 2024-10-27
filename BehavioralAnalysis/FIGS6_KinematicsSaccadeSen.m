clear all;close all;global gh
colorsel;
addpath('./stattool');addpath('./disptool');
for ii=1:3
    for jj=1:3
        DataGroup{ii,jj}=[];
    end
end

for fishsub = 1:length(totalfishsub)
    gh.param.fishid = totalfishsub(fishsub);
    fm_behavim_main
    for ii=1:sessionn(fishsub)
        if ~ismember(ii,gh.param.ExcludedSession)
            stim_boutn = find(gh.data.boutmtx(:,1)==ii&gh.data.boutmtx(:,2)>20 & gh.data.boutmtx(:,4)<=30);
            if gh.param.fishlog.trialdetails.trial(ii,1)==8
                group(ii) = 1;
            elseif gh.param.fishlog.trialdetails.trial(ii,1)==1 || gh.param.fishlog.trialdetails.trial(ii,1)==2
                group(ii) = 2;
            elseif gh.param.fishlog.trialdetails.trial(ii,1)==4 || gh.param.fishlog.trialdetails.trial(ii,1)==5
                group(ii) = 3;
            end

            for boutn=1:length(stim_boutn)
                boutbias = abs(gh.data.bout_details(stim_boutn(boutn),7)*-1- anglebias_overall(fishsub));

                X = gh.data.anglevect(gh.data.bout_details(stim_boutn(boutn),3):gh.data.bout_details(stim_boutn(boutn),5));
                X0 = X-max(X);
                X0_inv = -X0;
                pks_X0 = findpeaks(X0);
                pks_X0_inv = findpeaks(X0_inv);
                dur = gh.data.bout_details(stim_boutn(boutn),6);

                avgundulationinterval = dur./((length(pks_X0)+length(pks_X0_inv))./2);

                if boutn>=2
                    pretailflip_interval = [gh.data.bout_details(stim_boutn(boutn),2)-gh.data.bout_details(stim_boutn(boutn-1),2)];
                else
                    pretailflip_interval = nan;
                end
                ifchemosaccadecoupled = find(gh.data.simuMtx(:,9)==gh.data.boutmtx(stim_boutn(boutn),2) & ...
                    gh.data.simuMtx(:,1)==gh.data.boutmtx(stim_boutn(boutn),1));
                ifchemosaccadefollowing = find(gh.data.simuMtx(:,1)==gh.data.boutmtx(stim_boutn(boutn),1) & ...
                    gh.data.simuMtx(:,2)<gh.data.boutmtx(stim_boutn(boutn),2) & ...
                    gh.data.simuMtx(:,2)>20);
                if ~isempty(ifchemosaccadecoupled)
                    subgroup = 2;
                elseif ~isempty(ifchemosaccadefollowing)
                    subgroup = 3;
                else
                    subgroup = 1;
                end
                DataGroup{subgroup,group(ii)} = [DataGroup{subgroup,group(ii)}; boutbias, ...
                    avgundulationinterval,...
                    pretailflip_interval];
                clear ifchemosaccadecoupled ifchemosaccadefollowing boutbias avgundulationinterval pretailflip_interval
            end
        end
    end
end

min_val = [0, 0.032, 0.5];
max_val = [10, 0.048, 2];
int_val = [2, 0.002, 0.15];
weight_val = [.4, .5, .6];  
scaling_val = [.5, .0008, .06];
yticks_val{1} = [0:10];
yticks_val{2} = [.032:.002:.048];
yticks_val{3} = [.5:.1:2];
mkrsz = 25;
specifiedFaceColor = 'none';
specifiedColor{1} = blankColor ;
specifiedColor{2} = appeColor;
specifiedColor{3} = averColor;
specifiedalpha = 0.05;

for ii=1:length(min_val)
    ii
    count_observations(ii) = 0;
    figure(ii)
    allstatmtx{ii} =[];varstatmtx=[];varstatmtx2=[];
    for jj=1:3 %stimuli group
        for kk=1:3 %subgroup
            datamtx = DataGroup{kk,jj};
            datavec = datamtx(:,ii);
            count_observations(ii) = count_observations(ii)+sum(isnan(datavec)==0);
            varstatmtx(:,(jj-1)*3+kk) = [datavec;nan(1000-length(datavec),1)];
             
            xpos = jj+(kk-1)*0.25-0.25;
            
            limiter_vec = [min_val(ii):int_val(ii):max_val(ii)];
            limiter_vec_indx =[];prop=[];
            
            tempdata=max(min(datavec,max_val(ii)),min_val(ii));
            for mm=1:length(datavec)
                [val limiter_vec_indx(mm)] = min(abs(tempdata(mm)-limiter_vec));
            end
            for mm=1:length(limiter_vec)
                prop(mm)=length(find(limiter_vec_indx==mm))./length(datavec);
            end
            
            
            weight = prop(limiter_vec_indx)*weight_val(ii);
            randx = xpos+rand(1,length(datavec)).*weight-weight./2;
            
            scatter(randx,max(min(datavec,max_val(ii)),min_val(ii)),mkrsz,'MarkerEdgeColor',specifiedColor{jj},'MarkerFaceColor',specifiedColor{jj}, ...
                'MarkerFaceAlpha',specifiedalpha,'MarkerEdgeAlpha',specifiedalpha); hold on
            mid = nanmedian(datavec);
            u = quantile(datavec,0.75);
            l = quantile(datavec,0.25);
            plot([xpos-0.15./2 xpos+0.15./2],[mid mid],'color',specifiedColor{jj},'linewidth',1.5); hold on
            plot([xpos-0.06./2 xpos+0.06./2],[min(u,max_val(ii)) min(u,max_val(ii))],'color',specifiedColor{jj},'linewidth',1.5); hold on
            plot([xpos-0.06./2 xpos+0.06./2],[max(l,min_val(ii)) max(l,min_val(ii))],'color',specifiedColor{jj},'linewidth',1.5); hold on
                
            violin(xpos,max(min(datavec,max_val(ii)),min_val(ii)),'facecolor',specifiedColor{jj}, ...
                'scaling',scaling_val(ii),'facealpha',0.2,'style',2); hold on;
            
            mid =[];u=[];l=[];

        end
    end
    allstatmtx{ii} = varstatmtx;
    [p_kw{ii},tbl{ii},stats] =kruskalwallis(varstatmtx);
    c{ii} = multcompare(stats,"CType","dunn-sidak");
    stats=[];
    xticklabels([]);yticklabels([]);yticks(yticks_val{ii})
    xticks([0.75,1,1.25,1.75,2,2.25,2.75,3,3.25]);
    set(gcf,'Position',[(ii-1)*500+1 800 500 400])
    ylim([min_val(ii) max_val(ii)]);xlim([.5 3.5])
end




