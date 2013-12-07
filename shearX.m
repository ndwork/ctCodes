
function sheared = shearX( f, a )
  sh = shearY( f', a );
  sheared = sh';
end
