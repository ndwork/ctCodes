function masks = createLIweights(metalMask, extent)
%Create linear interpolation weighting masks
% Inputs:
%        metalMask -- mask where metal is identified with values >0
%        extent    -- furthest away in number of pixels to linearly
%                      interpolate (defines width of interpolation region)
% Outputs:
%         masks    -- masks(:,:,1) is the weighting image where the
%                      interpolation region ramps from 0 (away from metal)
%                      to 1 (next to metal) and is 0 everywhere else
%                     masks(:,:,2) is the same as masks(:,:,1) except the
%                      ramp goes from 0 to 1 in the opposite direction
%                     masks(:,:,3) is 1 outside the interpolation region
%                      and 0 inside the interpolation region and metal

[R,C] = size(metalMask);

[n,m]=find(metalMask>0);
masks=zeros(R,C,3);

for i=1:R
    for j=1:C
        d=sqrt((i-n).^2 + (j-m).^2);
        minD=min(d);
        if(minD<=extent)
            masks(i,j,1)=minD/extent;
            masks(i,j,2)=1-(minD/extent);       
        else
            masks(i,j,3)=1;
        end        
    end
end
masks(:,:,2)=(~metalMask).*masks(:,:,2);
