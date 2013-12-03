
function recon = marWavRad( sinogram, thetas,nDetectors, dSize, cx, cy, ...
  Nx, Ny, dx, dy, window )

  W = Wavelet;
  mkdir('tmp');

  % determine the sinogram mask
  recon = ctIRadon( sinogram, thetas, dSize, cx, cy, Nx, Ny, dx, dy, 'Hanning' );
  metalMask=findMetal(recon,50);
  sinoMask = ctRadon( metalMask, dx, nDetectors, dSize, thetas );
  sinoMask = ( sinoMask == 0 );
% load( 'recon.mat' );
% load( 'sinoMask.mat' );

  % scale sinogram so that max is 1
  sino = sinogram;
  sino = sino / max(sino(:));

  % pad sinogram so that of appropriate size for wavelet transform
  sinoY = size( sino, 1 );
  sinoX = size( sino, 2 );
  M = max( sinoY, sinoX );
  N = 2^ceil(log(M)/log(2));
  padX = max( N-sinoX, 0 );
  padY = max( N-sinoY, 0 );
  padSino = padarray( sino, [padY,padX], 'replicate', 'post'  );
  padMask = padarray( sinoMask, [padY,padX], 'replicate', 'post' );
  origSino = padSino;

  tolerance = 1d-3;
  diff = tolerance + 1;
  nIter = 0;
  while( diff > tolerance && nIter < 1000 )
    oldSino = padSino;

    wSino = W * padSino;
    st = softThresh( wSino, 0.05 );
    stSino = W' * st;

    padSino = ( padSino .* padMask ) + (stSino .* (1-padMask));

    sino = padSino(1:sinoY,1:sinoX);
    imwrite( sino, ['tmp/sino',num2str(nIter,'%4.4i'),'.jpg'], 'jpeg' );
    save( ['tmp/wSino',num2str(nIter,'%4.4i'),'.mat'], 'wSino' );

    diff = norm( padSino - oldSino );
    disp(['Difference Value: ', num2str(diff) ]);

    nIter = nIter + 1;
  end

  recon = ctIRadon( reconSino, thetas, dSize, cx, cy, Nx, Ny, dx, dy, 'Hanning' );
  
end
