
function ctTest
  clear; close all;

  delta = 0.001;
  object = createPhantom(delta, [0,0], 0.01);

  figure( 'name', 'Original Image' );
  imshow( object.im, [] );

  nDetectors = 500;
  dSize = 0.001;
  dTheta = 1 * pi/180;
  thetas = 0:dTheta:pi-dTheta;
  %sinogram = ctRadon( object.im, object.delta, nDetectors, dSize, thetas );
load 'sinogram.mat';

  figure( 'name', 'sinogram' );
  imshow( sinogram, [] );

  % Reconstruction parameters
  cx = 0;   Nx=256;   dx=delta;
  cy = 0;   Ny=256;   dy=delta;

  recon = ctIRadon( sinogram, thetas, dSize, cx, cy, Nx, Ny, dx, dy );
  figure( 'name', 'Reconstruction' );
  imshow( recon, [] );

  recon2DFT = ctIRadon2DFT( sinogram, thetas, dSize, cx, cy, Nx, Ny, dx, dy, 'none' );
  figure( 'name', '2DFT Reconstruction' );
  imshow( recon2DFT, [] );
%end
