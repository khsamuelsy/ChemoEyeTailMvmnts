function outputim=fm_funcim_dft(ImRaw)
for jj=1:5
    mip=mean(ImRaw,3);
    for ii=2:size(ImRaw,3)
        drifty = [];
        driftx = [];
        for kk=2:14
            refframe = fft2(mip(round(end*kk/16):round(end*(kk+1)/16),round(end*3/8):round(end*5/8)));
            movingframe = fft2(ImRaw(round(end*kk/16):round(end*(kk+1)/16),round(end*3/8):round(end*5/8),ii));
            output = dftregistration(refframe,movingframe,1);
            if abs(output(3))<=30
                drifty = [drifty;output(3)];
            end
            if abs(output(4))<=30
                driftx = [driftx;output(4)];
            end
        end
        output(3) = round(mean(drifty));
        output(4) = round(mean(driftx));
        if ~((output(3)==0) && (output(4)==0))
            if ((output(3)>=0) && (output(4)>=0))
                ImRaw(output(3)+1:end,output(4)+1:end,ii) = ...
                    ImRaw(1:end-output(3),1:end-output(4),ii);
            elseif ((output(3)>=0) && (output(4)<=0))
                ImRaw(output(3)+1:end,1:end+output(4),ii) = ...
                    ImRaw(1:end-output(3),-output(4)+1:end,ii);
            elseif ((output(3)<=0) && (output(4)>=0))
                ImRaw(1:end+output(3),output(4)+1:end,ii) = ...
                    ImRaw(-output(3)+1:end,1:end-output(4),ii);
            elseif ((output(3)<=0) && (output(4)<=0))
                ImRaw(1:end+output(3),1:end+output(4),ii) = ...
                    ImRaw(-output(3)+1:end,-output(4)+1:end,ii);
            end
        end
    end
end
outputim=ImRaw;
end