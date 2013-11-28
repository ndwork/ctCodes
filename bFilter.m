function outIm = bFilter(im, w, sigmDist, sigmInt)
%Bilateral filtering for grayscale images
% % Inputs:
%         im       -- input image
%         w        -- window size (window will be -w:w, -w:w)
%         sigmDist -- spatial domain standard deviation
%         sigmInt  -- intensity domain standard deviation
% % Outputs:
%         outIm    -- output image

% Normalize image to [0,1]
minIm=min(im(:));
maxIm=max(im(:));
im = im - minIm;
im = im./maxIm;    

%Calculate gaussian distance weights
[X,Y]=meshgrid(-w:w, -w:w);
G = exp(-(X.^2 + Y.^2)/(2*sigmDist^2));

outIm=zeros(size(im));

[R,C]=size(im);
for i=1:R
    for j=1:C        
        % Extract local region
        iMin=max(i-w, 1);
        iMax=min(i+w, R);
        jMin=max(j-w, 1);
        jMax=min(j+w, C);
        localReg = im(iMin:iMax, jMin:jMax);
        
        % Compute Gaussian intensity weights
        H=exp(-(localReg-im(i,j)).^2./(2*sigmInt^2));
        
        % Calculate bilateral filter response
        F=H.*G((iMin:iMax)-i+w+1, (jMin:jMax)-j+w+1);
        outIm(i,j) = sum(F(:).*localReg(:))/sum(F(:));
    end
end

outIm=outIm.*maxIm;
outIm=outIm+minIm;

