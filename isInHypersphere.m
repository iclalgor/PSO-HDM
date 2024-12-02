function y = isInHypersphere(P, Center, Radius)
% This Matlab function checks the points P whether is in a hypersphere or
% not
% P is mxn vector including m points in n dimensional spaces
[m,n] = size(P);
y=zeros(1,m);
for i=1:m
    if sqrt(sum((Center-P(i,:)).^2)) <= Radius
        y(i) = 1;
    else
        y(i) = 0;
    end
end
