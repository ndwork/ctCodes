function outIm = mdt(im, metalThresh)
% Metal Deletion Technique
% Inputs:
%         im          -- input image
%         metalThresh -- metal threshold value
% Outputs:
%         outIm -- output image
% MDT patent: US8233586


thetas=0:0.5:180;
% 1. 'original' projection data
imSino=radon(im,thetas);

sizeIm=size(im);

% 3. Linear interpolation
imLI = rubOut(im,metalThresh);

metalMask=findMetal(im,metalThresh);
metalSino=radon(metalMask,thetas);

% imCache=im.*metalMask + imLI.*(~metalMask);
masks=createLIweights(metalMask,10);
sinoMasks=createLIweights(metalSino,20);

imCache=im.*metalMask + imLI.*masks(:,:,3) + im.*masks(:,:,2) + imLI.*masks(:,:,1);

for iter=1:4
    disp(['MDT: iteration ' num2str(iter)]);
    % 4. edge preserving blur filter
    imCache=bFilter(imCache,4,2,0.1);
    
    % 5. forward project 4
    imCacheSino=radon(imCache, thetas);
    
    % 6. replace metal data from 1 with values from 5
%     imCacheSino=imSino.*(~(metalSino>0)) + imCacheSino.*(metalSino>0);
    imCacheSino = imCacheSino.*(metalSino>0) + imSino.*sinoMasks(:,:,3) + imCacheSino.*sinoMasks(:,:,2) + imSino.*sinoMasks(:,:,1);
    
    % 7. filtered back projection
    imCache=iradon(imCacheSino, thetas, sizeIm(1));
     
end

% Add back metal pixels from original image
outIm = im.*metalMask + imCache.*(~metalMask);
