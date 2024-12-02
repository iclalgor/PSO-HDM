%% Copyright (c) 2020, Mehdi Ghasri
%% To cooperate in articles, send an email to the following address
% (with Subject=CO Article):
% Email: Eng.mehdighasri@gmail.com


%% Copyright (c) 2019, Mehdi Ghasri
% All rights reserved. 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
% 
% * Redistributions of source code must retain the above copyright notice, this
%   list of conditions and the following disclaimer.
% 
% * Redistributions in binary form must reproduce the above copyright notice,
%   this list of conditions and the following disclaimer in the documentation
%   and/or other materials provided with the distribution
% * Neither the name of IUPUI nor the names of its
%   contributors may be used to endorse or promote products derived from this
%   software without specific prior written permission.
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
% DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
% FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
% DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
% SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
% CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
% OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
% OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
function [lb,ub,D,out] = Func_eng(F)


switch F
    case 'PressureVesselDesign'
        out = @PressureVesselDesign;
        D=4;
        lb=[0.0625 0.0625 10 10];
        ub=[99*0.0625 99*0.0625 200 200];
    case 'StringDesign'
        out = @StringDesign;
        D=3;
        lb=[0.05 0.25 2];
        ub=[2 1.3 15];
    case 'ThreeBarTruss'
        out = @ThreeBarTruss;
        D=2;
        lb=[0 0];
        ub=[1 1];
    case 'GearTrainDesign'
        out = @GearTrainDesign;
        D=4;
        lb=[12 12 12 12];
        ub=[60 60 60 60];
    case 'CantileverBeam'
        out = @CantileverBeam;
        D=5;
        lb=[0.01 0.01 0.01 0.01 0.01];
        ub=[100 100 100 100 100];
    case 'WeldedBeam'
        out = @WeldedBeam;
        D=4;
        lb=[0.1 0.1 0.1 0.1 ];
        ub=[2 10 10 2 ];
     
end
end
function out=PressureVesselDesign(x)

y1=x(:,1);%Ts
y2=x(:,2);%Th
y3=x(:,3);%R
y4=x(:,4);%L
%%% opt
fx=0.6224.*y1.*y3.*y4+...
    1.7781.*y2.*y3.^2+...
    3.1661.*y1.^2.*y4+...
    19.84.*y1.^2.*y3;
%%% const
g(:,1)=-y1+0.0193.*y3;
g(:,2)=-y2+0.0095.*y3;
g(:,3)=-pi.*y3.^2.*y4...
    -(4/3).*pi.*y3.^3 ...
    +1296000;
g(:,4)=y4-240;

%%% Penalty
pp=10^9;
for i=1:size(g,1)
    for j=1:size(g,2)
        if g(i,j)>0
            penalty(i,j)=pp.*g(i,j);
        else
            penalty(i,j)=0;
        end
    end
end

out=fx+sum(penalty,2);


end
function out=WeldedBeam(x)
y1=x(:,1);%W
y2=x(:,2);%L
y3=x(:,3);%d
y4=x(:,4);%h
%%% opt
fx=(y2.*1.1047.*y1.^2)+(0.04811.*y3.*y4.*(14+y2));
%%% const
sigm=504000./(y4.*y3.^2);
q=6000*(14+(y2./2));
D=0.5.*((y2.^2)+(y1+y3).^2).^0.5;
j=2*sqrt(2).*y1.*y2.*((y2.^2./6)+((y1+y3).^2)./2);
delta=65856./(30000.*y4.*y3.^3);
beta=(q.*D)./j;
alfa=6000./(sqrt(2).*y1.*y2);
toa=(alfa.^2+beta.^2+(alfa.*beta.*y2)./D).^0.5;
p=(0.61423*10^6).*((y3.*y4.^3)./6).*(1-(y3.*sqrt(y4.^6.*30/48))./28);

g(:,1)=toa-13600;
g(:,2)=sigm-30000;
g(:,3)=y1-y4;
g(:,4)=(0.1047.*y1.^2.*y2)+(0.04811.*y4.*y3.*(14+y2))-5;
g(:,5)=0.125-y1;
g(:,6)=delta-0.25;
g(:,7)=6000-p;
pp=10^9;
for i=1:size(g,1)
    for j=1:size(g,2)
        if g(i,j)>0
            penalty(i,j)=pp.*g(i,j);
        else
            penalty(i,j)=0;
        end
    end
end

out=fx+sum(penalty,2);

end

function out=StringDesign(x)

y1=x(:,1);%W
y2=x(:,2);%d
y3=x(:,3);%N
%%% opt
fx=(y3+2).*y2.*y1.^2;
%%% const
g(:,1)=1-(y2.^3.*y3)./(71785.*y1.^4);
g(:,2)=(4.*y2.^2-y1.*y2)./...
    (12566.*(y2.*y1.^3-y1.^4))...
    +(1./(5108.*y1.^2))-1;
g(:,3)=1-(140.45.*y1./(y2.^2.*y3));
g(:,4)=(y1+y2)./1.5-1;
%%% Penalty
pp=10^9;
for i=1:size(g,1)
    for j=1:size(g,2)
        if g(i,j)>0
            penalty(i,j)=pp.*g(i,j);
        else
            penalty(i,j)=0;
        end
    end
end

out=fx+sum(penalty,2);

end
function out=ThreeBarTruss(x)

A1=x(:,1);
A2=x(:,2);
%%%opt
fx=(2*sqrt(2).*A1+A2).*100;
%%% const
g(:,1)=2.*(sqrt(2).*A1+A2)./...
    (sqrt(2).*A1.^2+2.*A1.*A2)-2;
g(:,2)=2.*A2./(sqrt(2).*A1.^2+...
    2.*A1.*A2)-2;
g(:,3)=2./(A1+sqrt(2).*A2)-2;
%%% Penalty
pp=10^9;
for i=1:size(g,1)
    for j=1:size(g,2)
        if g(i,j)>0
            penalty(i,j)=pp.*g(i,j);
        else
            penalty(i,j)=0;
        end
    end
end

out=fx+sum(penalty,2);

end
function out=GearTrainDesign(x)

y1=x(:,1);%A
y2=x(:,2);%B
y3=x(:,3);%C
y4=x(:,4);%D
%%% opt
fx=((1/6.931)-((y1.*y2)./(y3.*y4))).^2;
out=fx;
end
function out=CantileverBeam(x)
y1=x(:,1);%1
y2=x(:,2);%2
y3=x(:,3);%3
y4=x(:,4);%4
y5=x(:,5);%5
%%% opt
fx=0.0624.*(y1+y2+y3+y4+y5);
%%%% const
g(:,1)=(61./y1.^3)+(37./y2.^3)+(19./y3.^3)+(7./y4.^3)+(1./y5.^3)-1;
%%% Penalty
pp=10^9;
for i=1:size(g,1)
    for j=1:size(g,2)
        if g(i,j)>0
            penalty(i,j)=pp.*g(i,j);
        else
            penalty(i,j)=0;
        end
    end
end

out=fx+sum(penalty,2);

end
