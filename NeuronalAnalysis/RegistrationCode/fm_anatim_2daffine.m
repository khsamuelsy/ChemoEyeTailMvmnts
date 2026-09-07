function [tform2d,newIm,ref2d_sr,im2d_sr]=fm_anatim_2daffine(Im,Ref)

Im_avg=mean(Im,3);
Ref_avg=mean(Ref,3);

ref2d_sr = imref2d(size(Ref_avg));
im2d_sr = imref2d(size(Im_avg));

[optimizer, metric] = imregconfig('multimodal');
optimizer.InitialRadius = 0.0009;
optimizer.Epsilon = 1.5e-4;
optimizer.GrowthFactor = 1.01;
optimizer.MaximumIterations = 300;

tform2d = imregtform(Im_avg, Ref_avg, "affine", optimizer, metric);
Im_avg_reg = imwarp(Im_avg, tform2d, 'OutputView', imref2d(size(Ref_avg)));

newIm = zeros(size(Ref));

for i = 1:size(newIm, 3)
    newIm(:, :,i) = imwarp(Im(:, :, i), tform2d, 'OutputView', imref2d(size(Ref_avg)));
end

end