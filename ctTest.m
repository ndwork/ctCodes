
function ctTest
  clear; close all;

  im = phantom();
  delta = 0.001;

  figure( 'name', 'Original Image' );
  imshow( im, [] );

  type = 'fast';

  nDetectors = 500;
  dSize = 0.001;
  dTheta = 1 * pi/180;
  thetas = 0:dTheta:pi-dTheta;
  sinogram = ctRadon( im, delta, nDetectors, dSize, thetas, type );

  figure( 'name', 'sinogram' );
  imshow( sinogram, [] );

  % Reconstruction parameters
  cx = 0;   Nx=256;   dx=delta;
  cy = 0;   Ny=256;   dy=delta;

  bp = ctBackProject( sinogram, thetas, dSize, cx, cy, Nx, Ny, ...
    dx, dy, type );
  figure( 'name', 'Back Projection' );
  imshow( bp, [] );

  recon = ctIRadon( sinogram, thetas, dSize, cx, cy, Nx, Ny, dx, dy, 'Hanning' );
  figure( 'name', 'Reconstruction' );
  imshow( recon, [] );

  recon2DFT = ctIRadon2DFT( sinogram, thetas, dSize, cx, cy, Nx, Ny, dx, dy, 'none' );
  figure( 'name', '2DFT Reconstruction' );
  imshow( recon2DFT, [] );

  sinogram2 = ctRadon( recon, delta, nDetectors, dSize, thetas );
  figure( 'name', 'Sinogram 2' );
  imshow( sinogram2, [] );
  
end
