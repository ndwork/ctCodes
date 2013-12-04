
function dsMatrix = makeDsMatrix( A, phi, x, y )

  Nx = numel(x);
  Ny = numel(y);
  nA = numel(A);
  nPhi = numel(phi);

  dsMatrix = sparse(nA*nPhi,Nx*Ny)

  parfor i=1:Nx*Ny
    dsIn = zeros(Ny,Nx);
    dsIn(i) = 1;
    dsOut = dsTransform( dsIn, A, phi, x, y );
    dsMatrix(:,i) = dsOut(:);
  end

end
