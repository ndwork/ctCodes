

function out = isoRot( img, theta )

  absTheta = abs(theta);
  n2Pis = floor( absTheta / (2*pi) );
  equivTheta = absTheta - n2Pis * 2*pi;
  nPiFourths = floor( abs(equivTheta) / (pi/4) );

  if sign( theta ) == 1
    out = img;

    for i=1:nPiFourths
      out = isoRot_theta( out, sign(theta)*pi/4 );
    end

    extra = equivTheta - nPiFourths*pi/4;
    out = isoRot_theta( out, sign(theta)*extra );

  else
    out = img;

    extra = equivTheta - nPiFourths*pi/4;
    out = isoRot_theta( out, -extra );

    for i=1:nPiFourths
      out = isoRot_theta( out, sign(theta)*pi/4 );
    end
  end
    
  %rotation = @(f,t)shearx( sheary( shearx(f,-tan(t/2)) ,sin(t)) ,-tan(t/2));  
end

function out = isoRot_theta( img, theta )
%function out = isoRot( img, theta )
  sh = shearX( img, -tan(theta/2) );
  sh = shearY( sh, sin(theta) );
  out = shearX( sh, -tan(theta/2) );
end
