function out = applyShearletTransform_arrayOutput(I,qmf1,qmf2,scale,ndir)

    [dd1 dd2] = sampled_DST(I,qmf1,qmf2,scale,ndir);
    
    out = cat(3,dd1,dd2);    

end