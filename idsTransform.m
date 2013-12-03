function T = idsTransform(S,x,y,A,phi)

noiseThresh = 0.02;

T=zeros(length(A),length(phi));

[R, C]=size(S);
if(R~=length(y) || C~=length(x))
    error('x or y is not the right size');
end

% if(~isequal(size(x),size(y)))
%     error('x and y must be the same ');
% end

ytraj=y;
yindxs=1:numel(y);

origS = S;

for i=length(A):-1:1
    disp(['idsTransform: A=' num2str(A(i))]);
    for j=1:length(phi)
        xtraj = A(i) * sin(y + phi(j));
        Svals = interp2(x,y,origS,xtraj,ytraj,'linear',0);
        minSval = min(Svals);
        if minSval <= noiseThresh continue; end;
        xindxs=interp1(x,1:numel(x),xtraj,'nearest');
        indxs1D = (xindxs-1)*R + yindxs;
        S(indxs1D)= S(indxs1D)-Svals;
        T(i,j) = mean(Svals);      
    end
end

end