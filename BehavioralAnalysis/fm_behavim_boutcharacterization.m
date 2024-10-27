function bout_details=fm_behavim_boutcharacterization

global gh

anglevect = gh.data.anglevect;
if isfield(gh.data,'boutmtx')
    boutmtx = gh.data.boutmtx;
else
    boutmtx = gh.PulledData{gh.param.fishseq}.boutmtx;
end
bout_details = boutmtx;

for ii=1:size(boutmtx,1)
    
    if ~isnan(boutmtx(ii,5))
        
        L = boutmtx(ii,5)-boutmtx(ii,3)+1;
        Fs = 200;
        T = 1/Fs;
        
        X = anglevect(boutmtx(ii,3):boutmtx(ii,5));
        
        X0 = X-max(X);
        X0_inv = -X0;
        pks_X0 = findpeaks(X0);
        pks_X0_inv = findpeaks(X0_inv);
        if ~isempty(pks_X0) && ~isempty(pks_X0_inv)
            bias = (mean(pks_X0)-mean([X0(1,end)]))-(mean(pks_X0_inv)-mean([X0_inv(1,end)]));
            mag = (mean(pks_X0)-mean([X0(1,end)]))+(mean(pks_X0_inv)-mean([X0_inv(1,end)]));
            
        elseif ~isempty(pks_X0)
            bias = (mean(pks_X0)-mean([X0(1,end)]));
            mag = (mean(pks_X0)-mean([X0(1,end)]));
            
        elseif ~isempty(pks_X0_inv)
            bias = -(mean(pks_X0_inv)-mean([X0_inv(1,end)]));
            mag = (mean(pks_X0_inv)-mean([X0_inv(1,end)]));
        end
        
        bout_details(ii,7) = bias;
        bout_details(ii,8) = mag;
    end
end

end