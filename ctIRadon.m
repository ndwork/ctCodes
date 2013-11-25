
function out = ctIRadon( sinogram, thetas, dSize, cx, cy, Nx, Ny, dx, dy )
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


  nThetas = size(sinogram,1);
  nDetectors = size(sinogram,2);


  %% Filter the sinogram
  n = ( [0:nDetectors-1] - floor(0.5*nDetectors) );
  %dLocs = n * dSize;

  h = zeros( 1, numel(n) );
  oddNs = find( mod(n,2)~=0 );
  h(oddNs) = -1 ./ ( n(oddNs).*n(oddNs) * pi*pi * dSize*dSize );
  zeroIndx = find( n==0 );
  h(zeroIndx) = 1/(4*dSize*dSize);
  
  nPadded = 2*nDetectors;
  hZP = zeros(1,nPadded);
  hZP(1:nDetectors) = h;
  hZP = ones(nThetas,1) * hZP;

  sinoZP = zeros(nThetas, nPadded);
  sinoZP(:,1:nDetectors) = sinogram;

  fftH = fft( hZP, [], 2 );
  fftSino = fft( sinoZP, [], 2 );

  fftFiltSino = fftH .* fftSino;
  filtSino = dSize * ifft( fftFiltSino, [], 2 );
  filtSino = circShift( filtSino, [0 -zeroIndx] );
  filtSino = filtSino(:,1:nDetectors);


  %% Perform backprojection
  out = ctBackProject( filtSino, thetas, dSize, cx, cy, Nx, Ny, dx, dy );
  out = out * pi / nThetas;
  
end
