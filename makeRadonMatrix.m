
function radonMatrix = makeRadonMatrix( Ny, Nx, delta, ...
  nDetectors, dSize, thetas )

  nThetas = numel( thetas );

  nPossible = 100000;
  rowIndxs = zeros(1,nPossible);
  colIndxs = zeros(1,nPossible);
  values = zeros(1,nPossible);

  nValues = 0;
  nIncrease = 100000;
  for i=1:Ny*Nx
    disp(['makeRadonMatrix ', num2str(i), ' of ', num2str(Ny*Nx)]);
    in = zeros( Ny, Nx );
    in(i) = 1;
    sino = ctRadon( in, delta, nDetectors, dSize, thetas );

    nonzeroIndxs = find( sino > 0 );
    nNonzero = numel( nonzeroIndxs );
    if nNonzero == 0 continue; end;

    if nValues + nNonzero > nPossible
      while nValues + nNonzero > nPossible
        newRowIndxs = zeros(1,nPossible+nIncrease);
        newColIndxs = zeros(1,nPossible+nIncrease);
        newValues = zeros(1,nPossible+nIncrease);
        nPossible = nPossible + nIncrease;
      end
      
      if nValues > 0
        newRowIndxs(1:nValues) = rowIndxs;
        newColIndxs(1:nValues) = colIndxs;
        newValues(1:nValues) = values;
      end
      rowIndxs = newRowIndxs;
      colIndxs = newColIndxs;
      values = newValues;

    end

    rowIndxs(nValues+1:nValues+nNonzero) = nonzeroIndxs;
    colIndxs(nValues+1:nValues+nNonzero) = i;
    values(nValues+1:nValues+nNonzero) = sino(nonzeroIndxs);
  end

save( 'rowIndxs.mat', 'rowIndxs' );
save( 'colIndxs.mat', 'colIndxs' );
save( 'values.mat', 'values' );
  
  radonMatrix = sparse( rowIndxs, colIndxs, values, Ny*Nx, ...
    nDetectors*nThetas );
  
save( 'radonMatrix.mat', 'R' );
end
