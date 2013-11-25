function metalMask=findMetal(im, metalThresh)

metalMask=zeros(size(im));
metalMask(im(:)>metalThresh)=1;

metalMask=imdilate(metalMask,strel('disk',1,8));
metalMask=imerode(metalMask,strel('disk',1,8));