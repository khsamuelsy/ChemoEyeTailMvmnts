function simuMtx = fm_behavim_simuBoutSaccadeID

global gh

simuMtx = [];
simuN = 0;


if isfield(gh.data,'saccademtx')
    saccademtx = gh.data.saccademtx;
else
    saccademtx = gh.PulledData{gh.param.fishseq}.saccademtx;
end

if isfield(gh.data,'bout_details')
    bout_details = gh.data.bout_details;
else
    bout_details = gh.PulledData{gh.param.fishseq}.bout_details;
end


for ii=1:size(saccademtx,1)
    eyen = ii;
    sessionn = saccademtx(ii,1);
    startt = saccademtx(ii,2);
    
    rown = find(bout_details(:,1)==sessionn & abs(bout_details(:,2)-startt)<=.5);
    if length(rown)~=0
        simuN = simuN+1;
        simuMtx(simuN,1:7) = saccademtx(ii,:);
        if length(rown)>1
            [value idx] = min(abs(bout_details(rown,2)-startt));
            simuMtx(simuN,8:15) = bout_details(rown(idx),:);
        elseif length(rown)==1
            simuMtx(simuN,8:15) = bout_details(rown,:);
        end
    else
        eyemtx(ii,7:12) = nan;
    end
    
end
end