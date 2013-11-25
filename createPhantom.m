function phntm=createPhantom(delta, metalLoc, metalR)
% delta is resolution of phantom in m
% metalLoc is 2 element array [x,y] indicating center of metal disk in m
% metalR   is radius of metal disk in m

if(length(metalLoc)~=2)
    return
end

[R, C]=size(phantom);

xlen=delta*C;
ylen=delta*R;

x=linspace(-(xlen/2), (xlen/2), R);
y=linspace(-(ylen/2), (ylen/2), C);

xs = ones(numel(y),1) * x;
ys = y' * ones(1,numel(x));

distFromR = sqrt( (xs-metalLoc(1)).^2 + (ys-metalLoc(2)).^2 );
metalMask = distFromR <= metalR;

phntm=struct;
phntm.im=phantom;
phntm.metalMask=metalMask;
phntm.delta=delta;
phntm.x=x;
phntm.y=y;





