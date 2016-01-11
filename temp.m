function [T] = temp(hx,tend)
%TEMP
close all;

a = 0.1;
v = 1;
L = 3;
Tcool = 5;
Thot = 100;

x = 0:hx:L;
ht = hx/v; %We use the "magic stepsize"
t = 0:ht:tend;

methods = {@upwind,@laxWen};
T = Tcool*ones(length(t),length(x),length(methods));
for i = 1:length(methods)
    meth = methods{i};
    T(:,:,i) = meth(T(:,:,i),hx,ht,a,v,Tcool,Thot);
    figure;
    surf(x,t,T(:,:,i));
    str = func2str(meth);
    title(sprintf('%s - hx = %f',str,hx));

end

end

function T = upwind(T,hx,ht,a,v,Tcool,Thot)
%Solves the problem with the upwind method
[n,m] = size(T);
coeff = v*ht/hx;
for k = 1:n-1
    T(k+1,1) = bcLeft(k*ht,Tcool,Thot);
    T(k+1,2:m) = (1-coeff-a*ht)*T(k,2:m) + coeff*T(k,1:m-1) + a*ht*Tcool;
end
end

function T = laxWen(T,hx,ht,a,v,Tcool,Thot)
%Solves the problem with the Lax-Wendroff method
[n,m] = size(T);
c1 = ht*v*(1+ht*v/hx-a*ht)/(2*hx);
c2 = 1-a*ht-v*v*ht*ht/(hx*hx)+a*a*ht*ht/2;
c3 = ht*v*(-1+ht*v/hx+a*ht)/(2*hx);
c4 = ht*a*(1-ht*a/2)*Tcool;
for k = 1:n-1
    T(k+1,1) = bcLeft(k*ht,Tcool,Thot);
    T(k+1,2:m-1) = c1*T(k,1:m-2) + c2*T(k,2:m-1) + c3*T(k,3:m) + c4;
    T(k+1,m) = 2*T(k+1,m-1) - T(k+1,m-2);
end
end

function T = bcLeft(t,Tcool,Thot)
%Gives the left boundary condition
if t <= 0.5
    T = Tcool + (Thot-Tcool)*sin(pi*t);
elseif t <= 4
    T = Thot;
else
    T = Thot + Tcool*sin(pi*(t-4));
end
end
