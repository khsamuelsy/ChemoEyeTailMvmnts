function boutmtx = fm_behavim_boutid

%% this function identifies discrete bout


global gh

load([gh.param.rootdir,'fm',num2str(gh.param.fishid),'\Registered\',num2str(gh.param.fishid),'_cleaned.mat']);
gh.data.anglevect = behavim.cleaneddata.angle;
boutmtx = [];
framepersession = 8000;
framepersecond = 200;
detectingbout = false;
boutn = 0;
for ii=5:length(gh.data.anglevect)
    % when there is difference in angle, raise flag and store bout start
    % info
    if ((mod(ii-1,framepersession)+1)>=5) && (gh.data.anglevect(ii)~=gh.data.anglevect(ii-1) && ~detectingbout)
        sessionn = floor((ii-1)./framepersession)+1;
        ifvalid = (gh.data.anglevect(ii-1)==gh.data.anglevect(ii-2)) ...
            && ((gh.data.anglevect(ii-2)==gh.data.anglevect(ii-3))...
         && (gh.data.anglevect(ii-3)==gh.data.anglevect(ii-4)));
        
        if ~ismember(sessionn,gh.param.ExcludedSession) && ifvalid
            detectingbout = true;
            starttimept = (mod((ii-2),framepersession)+1)./framepersecond;
            boutn = boutn+1;
            boutmtx(boutn,:) = [sessionn,starttimept,ii-1,nan,nan,nan];
        end
    end
    % when there is no more difference in angle, raise start bout flag and
    % store info
    if gh.data.anglevect(ii)==gh.data.anglevect(ii-1) && detectingbout
        if ii<(length(gh.data.anglevect)-2)
            % check for 3 consecutive no difference angle frames, raise end
            % bout flag and store info
            if gh.data.anglevect(ii+1)==gh.data.anglevect(ii) && gh.data.anglevect(ii+2)==gh.data.anglevect(ii)
                endtimept = (mod((ii-2),framepersession)+1)./framepersecond;
                boutmtx(boutn,4:6) = [endtimept,ii-1,endtimept-starttimept];
                detectingbout=false;
            end
        end
    end
    if mod(ii,framepersession)==0
        detectingbout=false;
    end
end
% filter bout that is <= 2 frames
boutmtx(find((boutmtx(:,5)-boutmtx(:,3)) <= 2),:)=[];
end