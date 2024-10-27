
switch figuren

    case 'elementalMerged'
        for ii=1:10
            if ii==1
                datamtx = mean(datasummary{ii}.ts(:,2:4)')' - datasummary{ii}.ts(:,1);
                size(datamtx)
                std_water = (std(datamtx)).*sqrt(length(datamtx))./sqrt(length(datamtx)-1);

            else
                datamtx = datasummary{ii}.ts(:,2:4) - datasummary{ii}.ts(:,1);
                for jj=1:3
                    [~,p_ztest_quadrant(ii,jj)]=ztest(datamtx(:,jj),0,std_water);
                end
                [~,p_ztest_quadrant_groupwise(ii,1)]=ztest([datamtx(:,1);datamtx(:,2);datamtx(:,3)],0,std_water);
                [~,p_ztest_quadrant_groupwise(ii,2)]=ztest([datamtx(:,1);datamtx(:,2);],0,std_water);
                [~,p_ztest_quadrant_groupwise(ii,3)]=ztest([datamtx(:,1);datamtx(:,3);],0,std_water);
                [~,p_ztest_quadrant_groupwise(ii,4)]=ztest([datamtx(:,2);datamtx(:,3);],0,std_water);
            end
        end


    case 'DMSOCtrl'

        %         std_water = 0.147283137598374;
        std_water = 0.159932295654260;

        for ii=1:3
            datamtx = datasummary{ii}.ts(:,2:4) - datasummary{ii}.ts(:,1);
            for jj=1:3
                [~,p_ztest_quadrant(ii,jj)]=ztest(datamtx(:,jj),0,std_water);
            end
            [~,p_ztest_quadrant_groupwise(ii,1)]=ztest([datamtx(:,1);datamtx(:,2);datamtx(:,3)],0,std_water);
            [~,p_ztest_quadrant_groupwise(ii,2)]=ztest([datamtx(:,1);datamtx(:,2);],0,std_water);
            [~,p_ztest_quadrant_groupwise(ii,3)]=ztest([datamtx(:,1);datamtx(:,3);],0,std_water);
            [~,p_ztest_quadrant_groupwise(ii,4)]=ztest([datamtx(:,2);datamtx(:,3);],0,std_water);
        end

    case 'WaterCtrl'

        %         std_water = 0.147283137598374;
        std_water = 0.159932295654260;
        for ii=1:3
            datamtx = datasummary{ii}.ts(:,2:4) - datasummary{ii}.ts(:,1);
            for jj=1:3
                [~,p_ztest_quadrant(ii,jj)]=ztest(datamtx(:,jj),0,std_water);
            end
            [~,p_ztest_quadrant_groupwise(ii,1)]=ztest([datamtx(:,1);datamtx(:,2);datamtx(:,3)],0,std_water);
            [~,p_ztest_quadrant_groupwise(ii,2)]=ztest([datamtx(:,1);datamtx(:,2);],0,std_water);
            [~,p_ztest_quadrant_groupwise(ii,3)]=ztest([datamtx(:,1);datamtx(:,3);],0,std_water);
            [~,p_ztest_quadrant_groupwise(ii,4)]=ztest([datamtx(:,2);datamtx(:,3);],0,std_water);
        end

end