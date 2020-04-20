function[y]=lin2circonv(x,N)% x be the linearly convoluted vector
    p=length(x)-N+1;
    y=zeros(1,N); %defines the output buffer
    for i=p:N
        y(i)=x(i);
    end
    for j=1:p-1
        y(j)=x(j)+x(N+j);
    end
    y = y';
end

