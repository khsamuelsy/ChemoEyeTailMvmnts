function  FIG2_Example
% clear all; close all
global gh
close all
drivename='E:\';

% trialn=24;fishchoice=4;
% trialn=29;fishchoice=4;
trialn=7;fishchoice=5; % for figure 2 and supp fig 3a

behavframepersession = 8000;
totalfishsub = [325 330 365 366 367];
fishn = totalfishsub(fishchoice);
anglebias_overall = [2.7117, 1.9589, 1.5048, 4.6064, 2.2587];

addpath([drivename,'TransitFMBehavior\fm_behavim\analysis'])
load([drivename,'TransitFMBehavior\fm',num2str(fishn),'\Registered\',num2str(fishn),'_cleaned.mat'])
load([drivename,'FM_IntegratedAnalysis\regressors\fm',num2str(fishn),'_ROI_XBlur.mat'])

gh.param.rootdir = [drivename,'TransitFMBehavior\'];
gh.param.framepersession = 8000;
gh.param.framepersecond = 200;
gh.param.fishid = num2str(fishn);
gh.param.ExcludedSession =[]; %can set excluded session to zero here
gh.data.anglevect = behavim.cleaneddata.angle;
gh.data.saccadevect = behavim.cleaneddata.angle_LE;
gh.data.saccadevect_R = behavim.cleaneddata.angle_RE;
gh.data.boutmtx = fm_behavim_boutid;
gh.data.bout_details = fm_behavim_boutcharacterization;
gh.data.saccademtx = fm_behavim_saccadeid;
gh.data.simuMtx = fm_behavim_simuBoutSaccadeID;


figuren = 'b';

%case a; saccade, a1: supp fig 3a
%case b: tail flip, b1: supp fig 3a

switch figuren

    case 'a'

        %initialize saccade vector
        conj_saccade = zeros(size(gh.data.saccadevect));

        trans_setting = 0.25; %0.25 for supp fig 3a, 0.75 for fig 2

        %create the values for the saccade vector
        rown = find(gh.data.saccademtx(:,1)==trialn);
        for ii=1:length(rown)
            if gh.data.saccademtx(rown(ii),7)==1
                conj_saccade(round(gh.data.saccademtx(rown(ii),3):gh.data.saccademtx(rown(ii),5)))=1;
            else
                conj_saccade(round(gh.data.saccademtx(rown(ii),3):gh.data.saccademtx(rown(ii),5)))=-1;
            end
        end
        clear rown

        %draw S-T events shadow
        rown = find(gh.data.simuMtx(:,1)==trialn);
        for ii=1:length(rown)
            startt = min(gh.data.simuMtx(rown(ii),3),gh.data.simuMtx(rown(ii),10));
            endt = max(gh.data.simuMtx(rown(ii),5),gh.data.simuMtx(rown(ii),12));
            pgon = polyshape([startt,endt,endt,startt]-((trialn-1)*behavframepersession +200*5),[1 1 -1 -1]);
            plot(pgon,'FaceColor',[191 152 149]./255,'FaceAlpha',trans_setting,'EdgeColor','none'); hold on
        end
        clear rown

        %drarw te saccade-onlz shadow
        rown = find(gh.data.saccademtx(:,1)==trialn);
        for ii=1:length(rown)
            if isempty(find(gh.data.simuMtx(:,1)==gh.data.saccademtx(rown(ii),1) & gh.data.simuMtx(:,2)==gh.data.saccademtx(rown(ii),2) ))
                startt = gh.data.saccademtx(rown(ii),3);
                endt = gh.data.saccademtx(rown(ii),5);
                pgon = polyshape([startt,endt,endt,startt]-((trialn-1)*behavframepersession +200*5),[20 20 -10 -10]);
                plot(pgon,'FaceColor',[170 179 202]./255,'FaceAlpha',trans_setting,'EdgeColor','none'); hold on
            end
        end
        clear rown

        plot(conj_saccade(((trialn-1)*behavframepersession +200*5+1):trialn*behavframepersession),'k','linewidth',1); hold on
        ylim([-1 1]); 
        xlim([2214-10 2292+10]); set(gcf,'Position',[100 100 300 100])
        % xlim([0 7000]); set(gcf,'Position',[100 100 500 100])
        ax1 = gca;                   % gca = get current axis
        ax1.YAxis.Visible = 'off';   % remove y-axis
        ax1.XAxis.Visible = 'off';   % remove x-axis
        box off

        clear conj_saccade

    case 'b'

        trans_setting = 0.25; %0.25 for supp fig 3a, 0.75 for fig 2

        rown = find(gh.data.simuMtx(:,1)==trialn);
        for ii=1:length(rown)
            startt = min(gh.data.simuMtx(rown(ii),3),gh.data.simuMtx(rown(ii),10));
            endt = max(gh.data.simuMtx(rown(ii),5),gh.data.simuMtx(rown(ii),12));
            pgon = polyshape([startt,endt,endt,startt]-((trialn-1)*behavframepersession +200*5),[20 20 -10 -10]);
            plot(pgon,'FaceColor',[191 152 149]./255,'FaceAlpha',trans_setting,'EdgeColor','none'); hold on
        end

        rown = find(gh.data.boutmtx(:,1)==trialn);
        for ii=1:length(rown)
            if isempty(find(gh.data.simuMtx(:,1)==gh.data.boutmtx(rown(ii),1) & gh.data.simuMtx(:,9)==gh.data.boutmtx(rown(ii),2) ))
                startt = gh.data.boutmtx(rown(ii),3);
                endt = gh.data.boutmtx(rown(ii),5);
                pgon = polyshape([startt,endt,endt,startt]-((trialn-1)*behavframepersession +200*5),[20 20 -10 -10]);
                plot(pgon,'FaceColor',[100 164 177]./255,'FaceAlpha',trans_setting,'EdgeColor','none'); hold on
            end
        end
        plot(behavim.cleaneddata.angle(((trialn-1)*behavframepersession +200*5+1):trialn*behavframepersession)*-1-anglebias_overall(fishchoice) ,'k','linewidth',1); hold on
        ylim([-10 20])

        % xlim([0 7000]); set(gcf,'Position',[100 100 500 100])
        xlim([2214-10 2292+10]); set(gcf,'Position',[100 100 300 100])
        ax1 = gca;                   % gca = get current axis
        ax1.YAxis.Visible = 'off';   % remove y-axis
        ax1.XAxis.Visible = 'off';   % remove x-axis
        box off

    case 'b1'
        plot(behavim.cleaneddata.angle(51251-30:51293+30)*-1-anglebias_overall(fishchoice) ,'k','linewidth',1);
        set(gcf,'Position',[100 100 500 100])
        ax1 = gca;                   % gca = get current axis
        ax1.YAxis.Visible = 'off';   % remove y-axis
        ax1.XAxis.Visible = 'off';   % remove x-axis
        box off


    case 'c'
        ROI_heatmap = zeros(35,size(ROI_Regressor_XBlur.dfof_session,2));
        count = 0;
        list = find(ROI_Regressor_XBlur.region.coarse==275);length(find(ROI_Regressor_XBlur.region.coarse==275))
        list = [list,find(ROI_Regressor_XBlur.region.coarse==1)]; length(find(ROI_Regressor_XBlur.region.coarse==1))
        list = [list,find(ROI_Regressor_XBlur.region.coarse==94)]; length(find(ROI_Regressor_XBlur.region.coarse==94))
        list = [list,find(ROI_Regressor_XBlur.region.coarse==114)]; length(find(ROI_Regressor_XBlur.region.coarse==114))
        list = [list,find(ROI_Regressor_XBlur.region.coarse==0)]; length(find(ROI_Regressor_XBlur.region.coarse==0))


        ROI_heatmap = ROI_Regressor_XBlur.dfof_session((trialn-1)*80+11:trialn*80,list);
        %ROI_heatmap = ROI_Regressor_XBlur.dfof_session((trialn-1)*40+6:trialn*40,list);
        imshow(ROI_heatmap',[0 0.1])


        set(gcf,'Position',[100 100 500 500])
        daspect([1 10 1])

end
end