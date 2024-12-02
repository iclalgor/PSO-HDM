function T0 = getBaseTime()
tic;
x= 0.55;
for i=1:1000000
   x = x + x; x=x/2; x=x*x; x=sqrt(x); x=log(x); x=exp(x); x=x/(x+2);
end
T0 = toc;
