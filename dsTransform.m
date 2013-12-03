
function out = dsTransform( img, A, phi, x, y )

  nA = size( img, 1 );
  nPhi = size( img, 2 );
  Nx = numel(x);
  Ny = numel(y);

  minX = min(x);
  maxX = max(x);

  out = zeros( Ny, Nx );

  for aIndx = 1:nA  
    disp(['dsTransform: A=' num2str(A(aIndx))]);

    parfor phiIndx = 1:nPhi
      T = img( aIndx, phiIndx );
      if T==0 continue; end;

      thisA = A;
      thisY = y;
      tmp = zeros(Ny,Nx);
      for yIndx = 1:Ny
        traj = thisA(aIndx) * sin( thisY(yIndx) + phi(phiIndx) );
        if traj < minX || traj > maxX continue; end;
        dist = abs( traj - x );
        [~, minIndx] = min( dist );
        %out(yIndx,minIndx) = out(yIndx,minIndx) + T;
        tmp(yIndx,minIndx) = T;
      end
      out = out + tmp;

    end
    
    imshow( out, [] );
    drawnow
  end

end
