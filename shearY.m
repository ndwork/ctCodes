
function sheared = shearY( img, a )

  [Ny Nx] = size(img);

  tx = [0:Nx/2-1 0 -Nx/2+1:-1]';
  ty = [0:Ny/2-1 0 -Ny/2+1:-1]';
  [X,Y] = meshgrid(tx,ty);

  shiftedF = fftshift(img);
  tmp = real( ifft( fft(shiftedF) .* exp(a*2i*pi*X.*Y/Ny) ) );
  sheared = ifftshift( tmp );

end
