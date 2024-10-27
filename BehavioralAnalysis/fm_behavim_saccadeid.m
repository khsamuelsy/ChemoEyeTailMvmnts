function saccademtx = fm_behavim_saccadeid

%% this function identifies discrete saccade
%% default avg frame 25


global gh

saccademtx = [];
gh.param.framepersession = 8000;
gh.param.framepersecond = 200;

if isfield(gh.data,'saccadevect')
    saccadevect = gh.data.saccadevect;
    saccadevect_R = gh.data.saccadevect_R;
else
    saccadevect = gh.PulledData{gh.param.fishseq}.saccadevect;
    saccadevect_R = gh.PulledData{gh.param.fishseq}.saccadevect_R;
end

detectingsaccade = false;
saccaden = 0;



for ii=26:length(saccadevect)
    
    if saccadevect(ii)~=saccadevect(ii-25) && ~detectingsaccade
        
        sessionn = floor((ii-1)./gh.param.framepersession)+1;
        
        if ~ismember(sessionn,gh.param.ExcludedSession) && (mod((ii-1),gh.param.framepersession)+1)>25
            detectingsaccade = true;
            starttimept = (mod((ii-12.5-1),gh.param.framepersession)+1)./gh.param.framepersecond;
            saccaden = saccaden+1;
            saccademtx(saccaden,:) = [sessionn,starttimept,ii-12.5,nan,nan,nan,nan];
        end
    end
    
    if saccadevect(ii)==saccadevect(ii-25) && detectingsaccade
        %         endtimept = (mod((ii+12.5-2),gh.param.framepersession)+2)./gh.param.framepersecond;
        %pending confirm
        endtimept = (mod((ii-12.5-1),gh.param.framepersession)+1)./gh.param.framepersecond;
        if ii<(length(saccadevect)-2)
            if saccadevect(ii+1)==saccadevect(ii) && saccadevect(ii+2)==saccadevect(ii)
                saccademtx(saccaden,4:6) = [endtimept,ii-12.5,endtimept-starttimept-1./gh.param.framepersecond];
                saccademtx(saccaden,7) = sign(saccadevect(ii)-saccadevect(round(saccademtx(saccaden,3)-25+12.5)));
                detectingsaccade=false;
            end
            if sign(saccadevect_R(ii)-saccadevect_R(round(saccademtx(saccaden,3)-25+12.5))) ~= saccademtx(saccaden,7)
                % signal L R incoherent direction
                                % [gh.param.fishid saccademtx(saccaden,1) (mod((ii-1),gh.param.framepersession)+1)]
                saccademtx(end,:) =[];
                saccaden = saccaden - 1;
            end
        end
    end
    if mod(ii,gh.param.framepersession)==0
        detectingsaccade=false;
    end
end
end