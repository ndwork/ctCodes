
function out = dsInterpLUT( dsIn, threshold, lut )

  metalMask = dsIn >= threshold;
  metalMask=imdilate(metalMask,strel('disk',1,8));
  metalMask=imdilate(metalMask,strel('disk',1,8));
  metalIndxs = find( metalMask );

  [nA,nPhi]=size(dsIn);
  [xs,ys]=meshgrid(1:nPhi, 1:nA);

  dsMasked = dsIn .* (1-metalMask);

  out = zeros(size(dsIn));
  for i=1:numel(metalIndxs)
    disp(['Working on metal index ', num2str(i), ' of ', ...
      num2str(numel(metalIndxs)) ]);
    x = xs(metalIndxs(i));
    y = ys(metalIndxs(i));

    dsInIndxs = find( lut(:,1)==y & lut(:,2)==x );
    if numel( dsInIndxs ) == 0 continue;

    aIndxs = lut(dsInIndxs,3);
    phiIndxs = lut(dsInIndxs,4);

    out(y,x) = sum( dsIn(aIndxs,phiIndxs) );

    %data = dsIn(aIndxs,phiIndxs);
    %data = data( find(data>0) );
    %out(y,x) = mean(data);
  end

  out = out + dsIn .* ( 1-metalMask );

end

