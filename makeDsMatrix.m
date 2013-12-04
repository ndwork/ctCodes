
function dsMatrix = makeDsMatrix( A, phi, x, y )

  Nx = numel(x);
  Ny = numel(y);
  nA = numel(A);
  nPhi = numel(phi);

  nNonzero = Ny * nA * nPhi;
  colIndxs = zeros( nNonzero, 1 );
  rowIndxs = zeros( nNonzero, 1 );
  values   = zeros( nNonzero, 1 );

  rowImg = (1:Ny)' * ones(1,Nx);
  colImg = ones(Ny,1) * (1:Nx);

  for i=1:nA*nPhi
    if mod(i,100)==0
      disp(['dsMatrix ', num2str(i), ' of ', num2str(nA*nPhi)]);
    end
    dsIn = zeros(nA,nPhi);
    dsIn(i) = 1;
    dsOut = dsTransform( dsIn, A, phi, x, y );
    nonzeroIndxs = find( dsOut ~= 0 );
    colIndxs( (i-1)*Ny + 1 : i*Ny ) = colImg( nonzeroIndxs );
    rowIndxs( (i-1)*Ny + 1 : i*Ny ) = rowImg( nonzeroIndxs );
    values( (i-1)*Ny + 1 : i*Ny ) = dsOut( nonzeroIndxs );
  end

  dsMatrix = sparse( rowIndxs, colIndxs, values );
end


