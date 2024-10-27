close all; clear all;

figure
figuren = 'elementalMerged';
chemicalColor = [[252 229 214];...
    [248 203 173];...
    [244 177 131]]./255;

switch figuren
    case 'elementalMerged'
        datasummary=FAA_Quadrant_extractdata_WaterDMSOMerged;
        xvec=[-.275 0 0.275];
        plot([0.5 10.5],[0 0],'color','k','linewidth',2); hold on
        for ii=1:10
            for kk=1:size(datasummary{ii}.ts,1)
                plot(ii+xvec,(datasummary{ii}.ts(kk,2:4)-datasummary{ii}.ts(kk,1)),'color',[0.85 0.85 0.85],'linewidth',1,'linestyle','-','Marker','o', ...
                    'MarkerFaceColor',[0.85 0.85 0.85],'MarkerEdgeColor',[0.7 0.7 0.7],'MarkerSize',2); hold on
            end
            plot(ii+xvec,median(datasummary{ii}.ts(:,2:4)-datasummary{ii}.ts(:,1)),'k','linewidth',1,'linestyle','-','Marker','none'); hold on
            for jj=1:3
                scatter(ii+xvec(jj),median(datasummary{ii}.ts(:,jj+1)-datasummary{ii}.ts(:,1)),'MarkerFaceColor',chemicalColor(jj,:),'MarkerEdgeColor','k');hold on
            end
        end

        yticks([-0.8:0.1:0.8]);yticklabels([]);        ylim(gca,[-.8 .8]);xlim(gca,[0.5 10.5]);box off
        h = gca;        h.XAxis.Visible = 'off';        set(gcf,'Position',[100 100 650 400])
        set(h, 'LineWidth', 2);

    case 'DMSOCtrl'
        datasummary=FAA_Quadrant_extractdata_DMSOCtrl;
        xvec=[-.275 0 0.275];
        plot([0.5 3.5],[0 0],'color','k','linewidth',2); hold on
        for ii=1:3
            for kk=1:size(datasummary{ii}.ts,1)
                plot(ii+xvec,(datasummary{ii}.ts(kk,2:4)-datasummary{ii}.ts(kk,1)),'color',[0.85 0.85 0.85],'linewidth',1,'linestyle','-','Marker','o', ...
                    'MarkerFaceColor',[0.85 0.85 0.85],'MarkerEdgeColor',[0.7 0.7 0.7],'MarkerSize',2); hold on
            end
            plot(ii+xvec,median(datasummary{ii}.ts(:,2:4)-datasummary{ii}.ts(:,1)),'k','linewidth',1,'linestyle','-','Marker','none'); hold on
            for jj=1:3
                scatter(ii+xvec(jj),median(datasummary{ii}.ts(:,jj+1)-datasummary{ii}.ts(:,1)),'MarkerFaceColor',chemicalColor(jj,:),'MarkerEdgeColor','k');hold on
            end
        end
        yticks([-0.8:0.1:0.8]);yticklabels([]);     ylim(gca,[-.8 .8]);xlim(gca,[0.5 3.5]);box off
        h = gca;
        h.XAxis.Visible = 'off';
        set(gcf,'Position',[100 100 300 400]);        set(h, 'LineWidth', 2);

    case 'WaterCtrl'
        datasummary=FAA_Quadrant_extractdata_WaterCtrl;
        xvec=[-.275 0 0.275];
        plot([0.5 3.5],[0 0],'color','k','linewidth',2); hold on
        for ii=1:3
            for kk=1:size(datasummary{ii}.ts,1)
                plot(ii+xvec,(datasummary{ii}.ts(kk,2:4)-datasummary{ii}.ts(kk,1)),'color',[0.85 0.85 0.85],'linewidth',1,'linestyle','-','Marker','o', ...
                    'MarkerFaceColor',[0.85 0.85 0.85],'MarkerEdgeColor',[0.7 0.7 0.7],'MarkerSize',2); hold on
            end
            plot(ii+xvec,median(datasummary{ii}.ts(:,2:4)-datasummary{ii}.ts(:,1)),'k','linewidth',1,'linestyle','-','Marker','none'); hold on
            for jj=1:3
                scatter(ii+xvec(jj),median(datasummary{ii}.ts(:,jj+1)-datasummary{ii}.ts(:,1)),'MarkerFaceColor',chemicalColor(jj,:),'MarkerEdgeColor','k');hold on
            end
        end
        yticks([-0.8:0.1:0.8]);yticklabels([]);        ylim(gca,[-.8 .8]);xlim(gca,[0.5 3.5]);box off
        h = gca;        h.XAxis.Visible = 'off';        set(gcf,'Position',[100 100 300 400])
        set(h, 'LineWidth', 2);
end

