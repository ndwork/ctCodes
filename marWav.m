
function recon = marWav( sinogram, thetas,nDetectors, dSize, cx, cy, ...
  Nx, Ny, dx, dy, window )

  W = wavelet;

%   recon = ctIRadon( sinogram, thetas, dSize, cx, cy, Nx, Ny, dx, dy, 'Hanning' );
load reconimage;
  
  metalMask=findMetal(recon,21000);
  
  sinoMask = ctRadon( metalMask, dx, nDetectors, dSize, thetas );
  sinoMask = ( sinoMask == 0 );
  
  recon=recon-min(recon(:));
  recon=recon./max(recon(:));

  tolerance = 1d-3;
  diff = tolerance + 1;

  while( diff > tolerance )
    oldRecon = recon;

    wRecon = W * recon;

    st = softThresh( wRecon,0.05 );

    stImg = W' * st;
    
    newSino = ctRadon( stImg, dx, nDetectors, dSize, thetas );

    reconSino = ( sinogram .* sinoMask ) + (newSino .* (1-sinoMask));

    recon = ctIRadon( reconSino, thetas, dSize, cx, cy, Nx, Ny, dx, dy, 'Hanning' );

    diff = norm( recon - oldRecon );
    disp(['Difference Value: ', num2str(diff) ]);
    imshow( recon, [] );  drawnow;
  end

end
