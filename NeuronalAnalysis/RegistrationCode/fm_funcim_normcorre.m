function output_stack = fm_funcim_normcorre(ImRaw)

%input is single ImRaw
%output is single ImRaw
padframe=1000;
ImRaw_avg = mean(ImRaw,3);
for ii=1:padframe
    ImRaw_padded(:,:,ii)=ImRaw_avg;
end
ImRaw_padded(:,:,(padframe+1):(size(ImRaw,3)+padframe))=ImRaw;

options_rigid = NoRMCorreSetParms('d1',size(ImRaw_padded,1),'d2',size(ImRaw_padded,2),'bin_width',20,'max_shift',15,'us_fac',50,'init_batch',200,'iter',1);
%50 (from 200) for max shift
%10 (from 20) for bin_width
[M1,shifts1,template1,options_rigid] = normcorre(ImRaw_padded,options_rigid);
options_nonrigid = NoRMCorreSetParms('d1',size(M1,1),'d2',size(M1,2),'grid_size',[32,32],'mot_uf',4,'bin_width',20,'max_shift',15,'max_dev',3,'us_fac',50,'init_batch',200,'iter',1);
%grid size 16, 16(from [32, 32])
[ImRaw_reg,shifts2,template2,options_nonrigid] = normcorre_batch(M1,options_nonrigid);
output_stack=ImRaw_reg(:,:,padframe+1:end);
end