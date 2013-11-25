
function bp = ctBackProject( sinogram, thetas, dSize, cx, cy, Nx, Ny, dx, dy )
  % Written by Nicholas Dwork
  % sinogram: each row corresponds to a theta, each column to a detector
  % thetas: a 1D array of thetas corresponding to each sinogram row
  % dSize: the size of each detector in meters
  % cx: center x position of reconstructed region
  % cy: center y position of reconstructed region
  % xFOV: horizontal field of view in meters
  % yFOV: vertical field of view in meters
  % Nx: number of pixels horizontally in reconstruction
  % Ny: number of pixels vertically in reconstruction

  sSino = size(sinogram);
  nThetas = sSino(1);
  nDetectors = sSino(2);

  dOffset = 0;   % detector center offset
  dLocs = ( [0:nDetectors-1] - floor(0.5*nDetectors) ) * dSize - dOffset;

  % Make arrays of x and y positions of each pixel
  lineXs = ( [0:Nx-1] - floor(0.5*Nx) ) * dx;
  lineYs = ( [0:Ny-1] - floor(0.5*Ny) ) * dy;
  NLxs = numel(lineXs);
  NLys = numel(lineYs);
  xs = ones(NLys,1) * lineXs;
  ys = lineYs' * ones(1,NLxs);
  xs=xs(:);
  ys=ys(:);

  angles = atan2(ys,xs);
  pixDs = sqrt( xs.*xs + ys.*ys );

  bp = zeros(Ny,Nx);
  parfor thIndx = 1:nThetas
    theta = thetas( thIndx );
    projections = pixDs .* cos( angles - theta );
    interped = interp1( dLocs, sinogram(thIndx,:), projections, ...
      'linear', 0);
    bp = bp + reshape(interped,Ny,Nx);
    if mod(thIndx,5)==0
      disp([ 'ctBackProject Theta Indx / Theta: ', ...
        num2str(thIndx), '/', num2str(theta) ]);
    end
  end

end

