function [tform3d,newIm,ref3d_sr,im3d_sr]=fm_anatim_3daffine(Im,Ref)


param.maxit = 1200; %1200

param.dper = 0.02;
param.initradius = 9e-4; %9e-4
param.epilson = 1.5e-4; %1.5e-4
param.growthfactor = 1.001;%1.001

[optimizer, metric] = imregconfig('multimodal');
optimizer.InitialRadius = param.initradius; %9e-4
optimizer.Epsilon = param.epilson; %1.5e-4
optimizer.GrowthFactor = param.growthfactor;%1.001
optimizer.MaximumIterations = param.maxit;

tform3d = imregtform(Im, Ref,'affine', optimizer, metric);

ref3d_sr = imref3d(size(Ref));
im3d_sr = imref3d(size(Im));

newIm = imwarp(Im, tform3d, 'OutputView', imref3d(size(Ref)));

end