function out = applyTransposeOfShearletTransform_arrayInput(in,qmf1,qmf2,scale,ndir)

dd1 = in(:,:,1:3); % don't hard-code in 3.
dd2 = in(:,:,4:6);

x1 = half_adj_sampled_DST1(dd1,qmf1,qmf2,scale,ndir,0);
x2 = half_adj_sampled_DST1(dd2,qmf1,qmf2,scale,ndir,1);
out = x1+x2;


end