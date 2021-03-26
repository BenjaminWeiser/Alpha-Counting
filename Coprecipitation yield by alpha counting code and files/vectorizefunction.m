function [y] = vectorizefunction(f,a,b,c)
%f is the function, 
%a is starting point
%b is space between each point
%c is ending point

n=1;
for i=a:b:c
    y(n,2)=f(i);
    y(n,1)=i;
    n=n+1;
end

