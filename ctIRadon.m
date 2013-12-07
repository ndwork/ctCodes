
function out = ctIRadon( sinogram, thetas, dSize, cx, cy, Nx, Ny, ...
  dx, dy, window, type )
  % Written by Nicholas Dwork
  % sinogram: each row corresponds to a theta, each column to a detector
  % thetas: a 1D array of thetas corresponding to each sinogram row
  % dSize: the size of each detector in meters
  % cx: center x position of reconstructed region
  % cy: center y position of reconstructed region
  % Nx: number of pixels horizontally in reconstruction
  % Ny: number of pixels vertically in reconstruction
  % dx: the horizontal size of each pixel in meters
  % dy: the vertical size of each pixel in meters
  % window: the window to be applied to the ramp filter
  % type: 'fast' or 'iso'

  defaultType = 'fast';
  expectedTypes = { 'iso', 'fast' };
  defaultWin = 'none';
  expectedWins = { 'Hanning', 'none' };
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
  p.addOptional('window',defaultWin, ...
    @(x) any(validatestring(x,expectedWins)) );
  p.addOptional('type',defaultType, ...
    @(x) any(validatestring(x,expectedTypes)) );

  p.parse( sinogram, thetas, dSize, cx, cy, Nx, Ny, ...
  dx, dy, window, type );
  inputs = p.Results;
  type = inputs.type;
  window = inputs.window;

  [nThetas nDetectors] = size(sinogram);

  %% Filter the sinogram
  cn = floor( 0.5 * nDetectors );
  n = ( [0:nDetectors-1] - cn );
  %dLocs = n * dSize;

  h = zeros( 1, numel(n) );
  oddNs = find( mod(n,2)~=0 );
  h(oddNs) = -1 ./ ( n(oddNs).*n(oddNs) * pi*pi * dSize*dSize );
  zeroIndx = find( n==0 );
  h(zeroIndx) = 1/(4*dSize*dSize);
  nPadded = 2*nDetectors;
  hZP = zeros(1,nPadded);
  hZP(1:nDetectors) = h;

  if nargin > 9
    if strcmp(window,'Hanning')

      df = 1 / ( dSize * nPadded );
      halfP = floor(nPadded / 2) + 1;
      fIndxs = ([1:nPadded] - halfP);
      f = df .* fIndxs;
      absF = abs(f);

      fN = 1./(2*dSize);
      hannWin = 0.5.*(1+cos(pi .* absF ./ fN));
      hannWin = hannWin .* ( absF < fN );
      fftHZP = fftshift( fft(hZP) );
      filtFftHZP = fftHZP .* hannWin;
      hZP =  ifft( ifftshift( filtFftHZP) );
    end
  end

  hZP = ones(nThetas,1) * hZP;
  sinoZP = zeros(nThetas, nPadded);
  sinoZP(:,1:nDetectors) = sinogram;

  fftH = fft( hZP, [], 2 );
  fftSino = fft( sinoZP, [], 2 );

  fftFiltSino = fftH .* fftSino;
  filtSino = dSize * ifft( fftFiltSino, [], 2 );
  filtSino = circshift( filtSino, [0 -zeroIndx+1] );
  filtSino = filtSino(:,1:nDetectors);


  %% Perform backprojection
  out = ctBackProject( filtSino, thetas, dSize, cx, cy, Nx, Ny, dx, dy, type );
  out = out * pi / nThetas;

end
