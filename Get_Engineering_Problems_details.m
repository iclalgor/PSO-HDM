%% Copyright (c) 2020, Mehdi Ghasri
%% To cooperate in articles, send an email to the following address
% (with Subject=CO Article):
% Email: Eng.mehdighasri@gmail.com


function [lb,ub,D,out] = Get_Engineering_Problems_details(F)

switch F

    case 'GearTrainDesign'
        out = @GearTrainDesign;
        D=4;
        lb=[12 12 12 12];
        ub=[60 60 60 60];
    
end
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

