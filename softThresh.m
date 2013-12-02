
function xHat = softThresh( y, lambda )
  if nargin<2 lambda=0.1; end;

  xHat = zeros(size(y));
  absY = abs(y);

  gtIndxs = find( absY > lambda );
  if numel(gtIndxs)>0
    xHat(gtIndxs) = ( absY(gtIndxs) - lambda ) ./ absY(gtIndxs) .* y(gtIndxs);
  end

end
