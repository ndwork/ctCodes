
function out = ctIRadonMetal( sinogram, thetas, dSize, cx, cy, Nx, Ny, ...
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
  
%   [rubbedSino,sinoMask]=rubOutSino(sinogram,50,thetas,nDetectors,dSize,cx,cy,Nx,Ny,dx,dy,filterWin);
load rubbedSino.mat;
load sinoMask.mat;


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
  sinoZP(:,1:nDetectors) = rubbedSino;

  fftH = fft( hZP, [], 2 );
  fftSino = fft( sinoZP, [], 2 );

  fftFiltSino = fftH .* fftSino;
  filtSino = dSize * ifft( fftFiltSino, [], 2 );
  filtSino = circshift( filtSino, [0 -zeroIndx+1] );
  filtSino = filtSino(:,1:nDetectors);

  
  

  %% Perform backprojection

  nonmetalSino = filtSino .* (1-sinoMask);
  
  sSino = size(filtSino);
  nThetas = sSino(1);
  nDetectors = sSino(2);

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
  nContributions = zeros(Ny,Nx);
  parfor thIndx = 1:nThetas
    theta = thetas( thIndx );
    projections = pixDs .* cos( angles - theta );
    interped = interp1( dLocs, nonmetalSino(thIndx,:), projections, ...
      'linear', 0);
    metalInterped = interp1( dLocs, sinoMask(thIndx,:), projections, ...
      'linear', 0);
    metalInterped = reshape( metalInterped, Ny, Nx );
    nContributions = nContributions + ( metalInterped == 0 );
    interped = reshape( interped, Ny, Nx );
    interped = interped .* ( metalInterped == 0 );
    bp = bp + interped;
    if mod(thIndx,10)==0
      disp([ 'ctIRadonMetal Theta Indx / Theta: ', ...
        num2str(thIndx), '/', num2str(theta) ]);
    end
  end

  nonmetalIndxs = find( nContributions == 0 );
  out = bp;
  out(nonmetalIndxs) = bp(nonmetalIndxs) * pi ./ ...
      nContributions(nonmetalIndxs);

end
