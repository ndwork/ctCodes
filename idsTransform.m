function T = idsTransform(S,x,y,A,phi)

T=zeros(length(A),length(phi));

[R, C]=size(S);
if(R~=length(y) || C~=length(x))
    error('x or y is not the right size');
end

% if(~isequal(size(x),size(y)))
%     error('x and y must be the same ');
% end

ytraj=y;

for i=1:length(A)
    disp(['idsTransform: A=' num2str(A(i))]);
    parfor j=1:length(phi)
        xtraj = A(i) * sin(y + phi(j));
        Svals = interp2(x,y,S,xtraj,ytraj,'linear');
        T(i,j) = mean(Svals);      
    end
end

end