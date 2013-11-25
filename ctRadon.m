
function sinogram = ctRadon( object, delta, nDetectors, dSize, thetas )
  dOffset=0;  % center channel offset

  sObject = size( object );
  nTheta = numel(thetas);

  sinogram = zeros( nTheta, nDetectors );

  dLocs = ( [0:nDetectors-1] - floor(0.5*nDetectors) ) * dSize - dOffset;

  thetas_deg = thetas * 180/pi;
  
  parfor th=1:numel(thetas)
    theta = thetas_deg(th);
    rotImg = imrotate( object, theta, 'bilinear','crop' );
    sumResult = sum( rotImg, 1 ) * delta;
    sRotImg = size( rotImg );
    locs = ( [0:sRotImg(2)-1] - floor(0.5*sRotImg(2)) ) * delta;
    interped = interp1( locs, sumResult, dLocs,'linear',0 );
    
    sinogram(th,:) = interped;
    if mod(th,5)==0 disp(['ctRadon Theta: ', num2str(th) ]); end;
  end

end
