
function recon = marShearRad( sinogram, thetas,nDetectors, dSize, cx, cy, ...
  Nx, Ny, dx, dy, window )

  mkdir('tmp_marShearRad');

  % determine the sinogram mask
  %recon = ctIRadon( sinogram, thetas, dSize, cx, cy, Nx, Ny, dx, dy, 'Hanning' );
  %metalMask=findMetal(recon,50);
  %sinoMask = ctRadon( metalMask, dx, nDetectors, dSize, thetas );
  %sinoMask = ( sinoMask == 0 );
load( 'recon.mat' );
load( 'sinoMask.mat' );

  % scale sinogram so that max is 1
  sino = sinogram;
  sino = sino - min(sino(:));
  sino = sino / max(sino(:));

  sino = imresize( sino, 0.5, 'bilinear' );
  sinoMask = imresize( sinoMask, 0.5, 'nearest' );

  sinoY = size( sino, 1 );
  sinoX = size( sino, 2 );

  useGPU = 0;
  %shearletSystem = SLgetShearletSystem2D( useGPU, sinoY, sinoX, 4 );
load shearletSystem.mat

  tolerance = 1d-3;
  diff = tolerance + 1;
  nIter = 0;
  origSino = sino;
  while( diff > tolerance && nIter < 1000 )
    oldSino = sino;

    shearSino = SLsheardec2D( sino, shearletSystem );
    st = softThresh( shearSino, 0.0005 );
    stSino = SLshearrec2D( st, shearletSystem );

    % Enforce the data constraint
    sino = ( origSino .* sinoMask ) + (stSino .* (1-sinoMask));

    imwrite( sino, ['tmp_marShearRad/sino',num2str(nIter,'%4.4i'),'.jpg'], 'jpeg' );
    %save( ['tmp_marShearRad/wSino',num2str(nIter,'%4.4i'),'.mat'], 'sino' );

    diff = norm( sino - oldSino );
    disp(['Difference Value: ', num2str(diff) ]);

    nIter = nIter + 1;
  end

  recon = ctIRadon( sino, thetas, dSize, cx, cy, Nx, Ny, dx, dy, 'Hanning' );
  
end
