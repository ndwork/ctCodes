function cost = computeCost_inpaintTightFrame_arrayVersion(x,rho,R,b,applyD)
% D is the analysis operator for the tight frame, which we assume here
% outputs an array.  (Not a cell array.)

% cost = || Dx ||_1 + .5*rho*||R.*x - R.*b||_F^2.

Dx = applyD(x);
cost = norm(Dx(:),1);

term = R.*x - R.*b;

cost = cost + .5*rho*norm(term(:))^2;

end