function  dsInterped   = dsInterp(dsIn, threshold)
% Nearest neighbor interpolation of metal region



[R, C]=size(dsIn);

[X,Y]=meshgrid(1:C, 1:R);

metalMask = dsIn >= threshold;

metalMask=imdilate(metalMask,strel('disk',1,8));
metalMask=imerode(metalMask,strel('disk',1,8));

metalIndxs = find( metalMask );
nonmetalIndxs = find( ~metalMask );

dsInterped=dsIn;


minIndxs = zeros(1,length(metalIndxs));
parfor mIndx=1:length(metalIndxs)
    if(mod(mIndx,10)==0)
        disp(['Interping ' num2str(mIndx) ' of ' num2str(length(metalIndxs)) ' points']);
    end
   
    dists = sqrt( ...
    ( X(metalIndxs(mIndx)) - X(nonmetalIndxs) ).^2 + ...
    ( Y(metalIndxs(mIndx)) - Y(nonmetalIndxs) ).^2 );

    [~,minIndx]=min(dists(:));
    minIndxs(mIndx) = nonmetalIndxs(minIndx);
end

dsInterped(metalIndxs)= dsIn(minIndxs);

end