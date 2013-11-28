function imOut=rubOut(im, metalThresh)

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

% clear
% clc
% close all
% load imMetalArt    %Siemens data
% im=imgctF2PiRadon;
% % load imMetalArtGE      %GE data
% % im=imMATLAB;
% 
% 
% metalThresh=0.18; %%Siemens data
% %  metalThresh=2;  %%GE data
metalMask=findMetal(im,metalThresh);

% imshow(metalMask,[]),title('metalMask')

thetas=0:0.5:360;
imSino=(radon(im,thetas))';

imSinoRubbed=imSino;
[L, num]=bwlabel(metalMask);
for n=1:num
    MmaskSino=radon((L==n),thetas);    
    MmaskSino=MmaskSino';
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

% figure,imshow(imSinoRubbed,[]),title('imSinoRubbed')
% figure,imshow(imSino,[]),title('imSino')

sizeIm=size(im);
imOut=iradon(imSinoRubbed',thetas,sizeIm(1));

% figure,imshow(imOut,[min(imOut(:)) max(imOut(:))]),title('imRubbedOut')
% figure,imshow(im, [min(imOut(:)) max(imOut(:))]),title('im')


