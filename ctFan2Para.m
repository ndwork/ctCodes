
function out = ctFan2Para( sinogram, gammas, betas, d, ...
  oThetas, oLines )
  % Written by Nicholas Dwork
  % gammas is a two element array specifying the min and max angle of the
  %   fan beam
  % betas is a two element array specifying the min and max rotation angles
  % d is the distance between the source and the isocenter
  % oBetas is a three element array specifying the min, delta, and max output
  %   rotation angles
  % oLines is the min, delta, and max spacing of projection lines

  sSino = size( sinogram );

  gammasImg = ones(sSino(1),1) * gammas;
  betasImg = betas' * ones(1,sSino(2));

  nThetas = round( ( oThetas(3) - oThetas(1) ) / oThetas(2) );
  nLines  = round( (  oLines(3) -  oLines(1) ) /  oLines(2) );

  thetas1 = linspace( oThetas(1), oThetas(3), nThetas );
  R1 = linspace( oLines(1), oLines(3), nLines );

  thetas = thetas1' * ones(1,nLines);
  Rs = ones(nThetas,1) * R1;

  oGammas = asin( Rs ./ d );
  oBetas = thetas + oGammas - pi/2;
%oGammas = -oGammas;

  minBetas = min(betas(:));
  maxBetas = max(betas(:));
  indxs = find( oBetas > maxBetas );
  if numel( indxs ) > 0 oBetas(indxs)=oBetas(indxs)-2*pi; end;
  indxs = find( oBetas < minBetas );
  if numel( indxs ) > 0 oBetas(indxs)=oBetas(indxs)+2*pi; end;
  oBetas = min( oBetas, maxBetas );

  out = interp2( gammas, betas, sinogram, oGammas, oBetas );
  
  nanIndxs = find( ~isfinite(out) );
  out(nanIndxs) = 0;
end
