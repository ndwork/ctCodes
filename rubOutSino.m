function [imSinoRubbed, sinoMask]=rubOutSino(sinogram, metalThresh, thetas, nDetectors,dSize, cx,cy,Nx,Ny, dx, dy, window);

% Inputs: 
%         im         -- image to apply RubOut on
%         metalTresh -- threshold level (all values above is considered metal)
% Outputs:
%         imOut      -- rubbed out image
                      

% Steps:
% Threshold im based on metalThresh. (Form binary metal mask)
% Dilate mask to close regions
% Forward project both metal mask and im
% Set imSino to 0 wherever MmaskSino=1
% Interpolate values for each projection to fill in black sinusoids
% Filtered back project rubbed imSino

delta=dx;
recon=ctIRadon(sinogram, thetas, dSize, cx, cy, Nx, Ny, delta, delta, window);

metalMask=findMetal(recon,metalThresh);

imSino=sinogram;

imSinoRubbed=imSino;
[L, num]=bwlabel(metalMask);
sinoMask=zeros(size(sinogram));

for n=1:num
    MmaskSino = ctRadon( (L==n), delta, nDetectors, dSize, thetas );
    sinoMask=(sinoMask + MmaskSino) > 0;

    for i=1:length(thetas)
        metalIndx=find(MmaskSino(i,:));
        LofM=metalIndx(1)-1;
        RofM=metalIndx(end)+1;
        width=RofM-LofM;
        Rweights=(metalIndx-LofM)./width;
        Lweights=(RofM-metalIndx)./width;
        imSinoRubbed(i,metalIndx)=imSino(i,LofM).*Lweights + imSino(i,RofM).*Rweights;
    end    
end



