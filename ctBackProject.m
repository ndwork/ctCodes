
function bp = ctBackProject( sinogram, thetas, dSize, cx, cy, Nx, Ny, ...
  dx, dy, varargin )
  % Written by Nicholas Dwork
  % sinogram: 2D array, each row corresponds to a theta, each column 
  %    to a detector
  % thetas: a 1D array of thetas corresponding to each sinogram row
  % dSize: the size of each detector in meters
  % cx: center x position of reconstructed region
  % cy: center y position of reconstructed region
  % xFOV: horizontal field of view in meters
  % yFOV: vertical field of view in meters
  % Nx: number of pixels horizontally in reconstruction
  % Ny: number of pixels vertically in reconstruction
  % dx: horizontal size of pixel in meters
  % dy: vertical size of pixel in meters
  % type (optional): 'iso' or 'fast'
  %    'iso' uses a rotation that's an isometry
  %    'fast' faster implementation

  defaultType = 'fast';
  expectedTypes = { 'iso', 'fast' };

  p = inputParser;
  p.addRequired('sinogram', @(x) ndims(x)==2);
  p.addRequired('thetas', @(x) ndims(x)==2);
  p.addRequired('dSize',@isnumeric);
  p.addRequired('cx',@isnumeric);
  p.addRequired('cy',@isnumeric);
  p.addRequired('Nx',@isnumeric);
  p.addRequired('Ny',@isnumeric);
  p.addRequired('dx',@isnumeric);
  p.addRequired('dy',@isnumeric);
  p.addOptional('type',defaultType, ...
    @(x) any(validatestring(x,expectedTypes)) );

  p.parse( sinogram, thetas, dSize, cx, cy, Nx, Ny, ...
    dx, dy, varargin{:} );
  inputs = p.Results;
  type = inputs.type;

  if strcmp(type, 'fast')
    bp = ctBackProjectFast( sinogram, thetas, dSize, cx, cy, Nx, Ny, ...
      dx, dy );
  elseif strcmp(type, 'iso')
    bp = ctBackProjectIso( sinogram, thetas, dSize, cx, cy, Nx, Ny, ...
      dx, dy );
  end

end

function bp = ctBackProjectIso( sinogram, thetas, dSize, cx, cy, Nx, Ny, ...
  dx, dy )

  if dx ~= dy
    err = MException('ctCodes:ctBackProjectIso', ...
        'dx and dy must be equal for this version of back projection');
    throw(err);
  end
  
  [nThetas, nDetectors] = size(sinogram);

  dOffset = 0;   % detector center offset
  dLocs = ( [0:nDetectors-1] - floor(0.5*nDetectors) ) * dSize - dOffset;

  % create the adjoint of the interpolation matrix
  if mod( Nx, 2 )==0
    locs = ( ([0:Nx-1]) - 0.5*Nx + 0.5 ) * dx;
  else
    locs = ( ([0:Nx-1]) - floor(0.5*Nx) ) * dx;
  end
  B = zeros(nDetectors, Nx);
  for i=1:Nx
    tmp = zeros(1,Nx);
    tmp(i) = 1;
    interped = interp1( locs, tmp, dLocs,'linear',0 );
    B(:,i) = interped;
  end
  Bt = B';

  bp = zeros(Ny,Nx);
  parfor thIndx = 1:nThetas
    theta = thetas( thIndx );
    line = sinogram(thIndx,:);
    tmp = Bt * line';
    smeared = ones(Ny,1) * tmp';
    bp = bp + isoRot( smeared, -theta );
    if mod(thIndx,5)==0
      disp([ 'ctBackProject Theta Indx / Theta: ', ...
        num2str(thIndx), '/', num2str(theta) ]);
    end
  end

end


function bp = ctBackProjectFast( sinogram, thetas, dSize, cx, cy, Nx, Ny, ...
  dx, dy )

  [nThetas, nDetectors] = size(sinogram);

  dOffset = 0;   % detector center offset
  dLocs = ( [0:nDetectors-1] - floor(0.5*nDetectors) ) * dSize - dOffset;

  % Make arrays of x and y positions of each pixel
  if mod( Nx, 2 )==0
    lineXs = ( ([0:Nx-1]) - 0.5*Nx + 0.5 ) * dx + cx;
  else
    lineXs = ( ([0:Nx-1]) - floor(0.5*Nx) ) * dx + cx;
  end
  if mod( Ny, 2 )==0
    lineYs = ( ([0:Ny-1]) - 0.5*Ny + 0.5 ) * dy + cy;
  else
    lineYs = ( ([0:Ny-1]) - floor(0.5*Ny) ) * dy + cy;
  end
  xs = ones(Ny,1) * lineXs;
  ys = lineYs' * ones(1,Nx);
  xs=xs(:) + cx;
  ys=ys(:) + cy;

  angles = atan2(ys,xs);
  pixDs = sqrt( xs.*xs + ys.*ys );

  bp = zeros(Ny,Nx);
  parfor thIndx = 1:nThetas
    theta = thetas( thIndx );
    projections = pixDs .* cos( angles - theta );
    interped = interp1( dLocs, sinogram(thIndx,:), projections, ...
      'linear', 0);
    bp = bp + reshape(interped,Ny,Nx);
    if mod(thIndx,10)==0
      disp([ 'ctBackProject Theta Indx / Theta: ', ...
        num2str(thIndx), '/', num2str(theta) ]);
    end
  end

end
