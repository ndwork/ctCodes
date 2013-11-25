
function recon = ctIRadon2DFT( sinogram, thetas, dSize, cx, cy, Nx, Ny, dx, dy, ...
  window )
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
  % window can be 'none', 'ramp', or 'Hanning'

  bp = ctBackProject( sinogram, thetas, dSize, cx, cy, 2*Nx, 2*Ny, dx, dy );


  %% Now filter the back projected image
  sBP = size(bp);
  dky = 1 / ( dSize * sBP(1) );
  dkx = 1 / ( dSize * sBP(2) );

  halfP = floor(sBP / 2) + 1;
  kxIndxs = ones(sBP(2),1) * ([1:sBP(1)] - halfP(1));
  kyIndxs = ([1:sBP(2)] - halfP(2))' * ones(1,sBP(1));

  kx = dkx .* kxIndxs;
  ky = dky .* kyIndxs;

  absK = sqrt( kx.*kx + ky.*ky );
  filter = absK;    % This is the ramp filter

  if strcmp(window,'Hanning')
    fN = 1./(2*dSize);
    hannWin = 0.5.*(1+cos(pi .* absK ./ fN));
    hannWin = hannWin .* ( absK < fN );
    filter = filter .* hannWin;
  end

  bpFFT = fftshift(fft2(bp));
  filtBpFFT = bpFFT .* filter;
  reconBig = real( ifft2( ifftshift(filtBpFFT) ) );

  half = sBP/2;
  mid = floor( half ) + 1;
  left = floor(Nx/2);
  right = floor((Nx-1)/2);
  bottom = floor(Ny/2);
  top = floor((Ny-1)/2);
  recon = reconBig( mid(1)-left:mid(1)+right, ...
                    mid(2)-bottom:mid(2)+top  );

  recon = recon*pi/numel(thetas);
end
