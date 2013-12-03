
function out = ctIRadon( sinogram, thetas, dSize, cx, cy, Nx, Ny, ...
  dx, dy, filterWin )
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
  % filterWin: the window to be applied to the ramp filter


  nThetas = size(sinogram,1);
  nDetectors = size(sinogram,2);

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
    if strcmp(filterWin,'Hanning')

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
  out = ctBackProject( filtSino, thetas, dSize, cx, cy, Nx, Ny, dx, dy );
  out = out * pi / nThetas;

end
