
function recon = marWav( sinogram, thetas,nDetectors, dSize, cx, cy, ...
  Nx, Ny, dx, dy, window )

  W = Wavelet;
  mkdir('tmp');

  recon = ctIRadon( sinogram, thetas, dSize, cx, cy, Nx, Ny, dx, dy, 'Hanning' );

  metalMask=findMetal(recon,50);
  sinoMask = ctRadon( metalMask, dx, nDetectors, dSize, thetas );
  sinoMask = ( sinoMask == 0 );

  recon=recon-min(recon(:));
  recon=recon./max(recon(:));
  reconSinogram = ctRadon( recon, dx, nDetectors, dSize, thetas );

  tolerance = 1d-3;
  diff = tolerance + 1;
  nIter = 0;
  while( diff > tolerance && nIter < 1000 )
    oldRecon = recon;

    wRecon = W * recon;
    st = softThresh( wRecon, 0.05 );
    stImg = W' * st;

    newSino = ctRadon( stImg, dx, nDetectors, dSize, thetas );

    reconSino = ( reconSinogram .* sinoMask ) + (newSino .* (1-sinoMask));
    recon = ctIRadon( reconSino, thetas, dSize, cx, cy, Nx, Ny, dx, dy, 'Hanning' );

    outSino = reconSino - min(reconSino(:));
    outSino = outSino / max(outSino(:));
    imwrite( outSino, ['tmp/reconSino',num2str(nIter,'%4.4i'),'.jpg'], 'jpeg' );
    imwrite( recon, ['tmp/recon',num2str(nIter,'%4.4i'),'.jpg'], 'jpeg' );

    diff = norm( recon - oldRecon );
    disp(['Difference Value: ', num2str(diff) ]);

    nIter = nIter + 1;
  end

end
