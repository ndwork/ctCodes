
function sinogram = ctRadon( img, delta, nDetectors, dSize, thetas, ...
  varargin )
  % img:  2D array - will take the Radon transform of this image
  % delta: horizontal and vertical size of pixel (assumed square)
  % nDetectors: the number of detectors
  % thetas: a 1D array, each element is the angle that corresponds to row
  %    radon domain
  % type (optional): 'iso' or 'fast'
  %    'iso' (default) uses a rotation that's an isometry
  %    'fast' faster implementation

  defaultType = 'fast';
  expectedTypes = { 'iso', 'fast' };

  p = inputParser;
  p.addRequired('img', @(x) ndims(x)==2);
  p.addRequired('delta',@isnumeric);
  p.addRequired('nDetectors',@isnumeric);
  p.addRequired('dSize',@isnumeric);
  p.addRequired('thetas');
  p.addOptional('type',defaultType, ...
    @(x) any(validatestring(x,expectedTypes)) );

  p.parse( img, delta, nDetectors, dSize, thetas, ...
    varargin{:} );
  inputs = p.Results;
  type = inputs.type;

  if strcmp( type, 'fast' )
    sinogram = ctRadonFast( img, delta, nDetectors, dSize, thetas );
  elseif strcmp( type, 'iso' )
    sinogram = ctRadonIso( img, delta, nDetectors, dSize, thetas );
  end
end


function sinogram = ctRadonFast( img, delta, nDetectors, dSize, thetas )
  dOffset=0;  % center channel offset

  nTheta = numel(thetas);

  dLocs = ( [0:nDetectors-1] - floor(0.5*nDetectors) ) * dSize - dOffset;

  thetas_deg = thetas * 180/pi;

  Ny = size( img, 1 );  halfY = Ny/2;
  Nx = size( img, 2 );  halfX = Nx/2;
  xs = ones(Ny,1) * (1:Nx);
  ys = (1:Ny)' * ones(1,Nx);
  xs = xs - halfX;
  ys = ys - halfY;
  radiusMask = sqrt( xs.*xs + ys.*ys ) < min(Nx/2,Ny/2);
  radiusImg = img .* radiusMask;

  if mod( Nx, 2 )==0
    locs = ( ([0:Nx-1]) - 0.5*Nx + 0.5 ) * delta;
  else
    locs = ( ([0:Nx-1]) - floor(0.5*Nx) ) * delta;
  end

  sinogram = zeros( nTheta, nDetectors );
  parfor th=1:numel(thetas)
    theta = thetas_deg(th);
    rotImg = imrotate( radiusImg, theta, 'bilinear','crop' );
    sumResult = sum( rotImg, 1 ) * delta;

    interped = interp1( locs, sumResult, dLocs,'linear',0 );

    sinogram(th,:) = interped;
    if mod(th,10)==0 disp(['ctRadon Theta: ', num2str(th), ' of ', ...
        num2str(numel(thetas)) ]); end;
  end

end

function [sinogram B] = ctRadonIso( img, delta, nDetectors, dSize, thetas )
  dOffset=0;  % center channel offset

  nTheta = numel(thetas);

  dLocs = ( [0:nDetectors-1] - floor(0.5*nDetectors) ) * dSize - dOffset;

  Ny = size( img, 1 );  halfY = Ny/2;
  Nx = size( img, 2 );  halfX = Nx/2;
  xs = ones(Ny,1) * (1:Nx);
  ys = (1:Ny)' * ones(1,Nx);
  xs = xs - halfX;
  ys = ys - halfY;
  radiusMask = sqrt( xs.*xs + ys.*ys ) < min(Nx/2,Ny/2);
  radiusImg = img .* radiusMask;

  if mod( Nx, 2 )==0
    locs = ( ([0:Nx-1]) - 0.5*Nx + 0.5 ) * delta;
  else
    locs = ( ([0:Nx-1]) - floor(0.5*Nx) ) * delta;
  end

  % create the interpolation matrix
  B = zeros(nDetectors, Nx);
  for i=1:Nx
    tmp = zeros(1,Nx);
    tmp(i) = 1;
    interped = interp1( locs, tmp, dLocs,'linear',0 );
    B(:,i) = interped;
  end


  sinogram = zeros( nTheta, nDetectors );
  parfor th=1:numel(thetas)
    theta = thetas(th);
    rotImg = isoRot( radiusImg, theta );
    sumResult = sum( rotImg, 1 ) * delta;

    %interped = interp1( locs, sumResult, dLocs,'linear',0 );
    interped = B * sumResult';

    sinogram(th,:) = interped;
    if mod(th,10)==0 disp(['ctRadon Theta: ', num2str(th), ' of ', ...
        num2str(numel(thetas)) ]); end;
  end

end


