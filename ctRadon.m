
function sinogram = ctRadon( img, delta, nDetectors, dSize, thetas )
  dOffset=0;  % center channel offset

  nTheta = numel(thetas);

  dLocs = ( [0:nDetectors-1] - floor(0.5*nDetectors) ) * dSize - dOffset;

  thetas_deg = thetas * 180/pi;

  Nx = size( img, 2 );
  if mod( Nx, 2 )==0
    locs = ( ([0:Nx-1]) - 0.5*Nx + 0.5 ) * delta;
  else
    locs = ( ([0:Nx-1]) - floor(0.5*Nx) ) * delta;
  end

  sinogram = zeros( nTheta, nDetectors );
  parfor th=1:numel(thetas)
    theta = thetas_deg(th);
    rotImg = imrotate( img, theta, 'bilinear','crop' );
    sumResult = sum( rotImg, 1 ) * delta;

    interped = interp1( locs, sumResult, dLocs,'linear',0 );

    sinogram(th,:) = interped;
    if mod(th,5)==0 disp(['ctRadon Theta: ', num2str(th) ]); end;
  end

end
