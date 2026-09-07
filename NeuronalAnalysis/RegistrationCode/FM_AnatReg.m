function FM_AnatReg(fishn)

addpath('./saveastiff_4.5')

ImLs = double(loadtiff(['fm',num2str(fishn),'_anatstack.tif']));
ImRef= double(loadtiff('refstack_elavl3_h2brfp.tif'));

[regdata.tform2d1,im.ls2daffine1,~,~] = fm_anatim_2daffine(ImLs,ImRef);
saveastiff(uint16(im.ls2daffine1),['fm',num2str(fishn),'_2dReg1 .tif']);

[regdata.tform3d1,im.ls3daffine1,~,~] = fm_anatim_3daffine(im.ls2daffine1,ImRef);
saveastiff(uint16(im.ls3daffine1),['fm',num2str(fishn),'_3dReg1.tif']);

[regdata.tform2d2,im.ls2daffine2,~,~] = fm_anatim_2daffine(im.ls3daffine1,ImRef);
saveastiff(uint16(im.ls2daffine2),['fm',num2str(fishn),'_2dReg2.tif']);

[regdata.tform3d2,im.ls3daffine2,~,~] = fm_anatim_3daffine(im.ls2daffine2,ImRef);
saveastiff(uint16(im.ls3daffine2),['fm',num2str(fishn),'_3dReg2.tif']);

save(['fm',num2str(fishn),'_regdata2.mat'],'regdata')

end

