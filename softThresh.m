function xhat=softThresh(y,lambda)

xhat=zeros(size(y));

for i=1:length(y)
    if(y(i) < -1*lambda)
        xhat(i)=y(i)+lambda;
    elseif(y(i) > lambda)
        xhat(i)=y(i)-lambda;
    elseif(abs(y(i))<lambda)
        xhat(i)=0;
    end
end