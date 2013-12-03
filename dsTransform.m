
function out = dsTransform( img, A, phi, x, y )

  nA = size( img, 1 );
  nPhi = size( img, 2 );
  Nx = numel(x);
  Ny = numel(y);

  minX = min(x);
  maxX = max(x);

  out = zeros( Ny, Nx );

  for aIndx = 1:nA
    for phiIndx = 1:nPhi
      T = img( aIndx, phiIndx );

      for yIndx = 1:Ny
        traj = A * sin( y(yIndx) + phi );
        if traj < minX || traj > maxX continue; end;
        dist = abs( traj - x );
        [junk minIndx] = min( dist );
        out(yIndx,minIndx) = out(yIndx,minIndx) + T;
      end

    end
  end

end
