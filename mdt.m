function outIm = mdt(im, metalThresh)
% Metal Deletion Technique
% Inputs:
%         im          -- input image
%         metalThresh -- metal threshold value
% Outputs:
%         outIm -- output image


thetas=0:0.5:360;
% 1. 'original' projection data
imSino=radon(im,thetas);

sizeIm=size(im);

% 3. Linear interpolation
imLI = rubOut(im,metalThresh);

metalMask=findMetal(im,metalThresh);
metalSino=radon(metalMask,thetas);

imCache=im.*metalMask + imLI.*(~metalMask);
%TODO: change above so that transition between im and imLI is smooth

for iter=1:4
    disp(['MDT: iteration ' num2str(iter)]);
    % 4. edge preserving blur filter
    imCache=bFilter(imCache,4,2,0.1);
    
    % 5. forward project 4
    imCacheSino=radon(imCache, thetas);
    
    %6. replace metal data from 1 with values from 5
    imCacheSino=imSino.*(~(metalSino>0)) + imCacheSino.*(metalSino>0);
    
    %7. filtered back projection
    imCache=iradon(imCacheSino, thetas, sizeIm(1));
     
end

% Add back metal pixels from original image
outIm = im.*metalMask + imCache.*(~metalMask);
