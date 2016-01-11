function [] = transport()
%TRANSPORT
close all;
clear all;

T = 1;
a = 1;
N = 100;

hx = 2/N;
x = linspace(0,2,N+1);
sigma = [0.8 1 1.1];
methods = {@upwind,@laxFried,@laxWen};
figure;
for i = 1:3
    meth = methods{i};
    for j = 1:3
        sig =sigma(j);
        ht = hx*sig/a;
        t = 0:ht:2*T;
        m = length(t);
        u = meth(t,sig,N,m);
        subplot(3,3,3*(i-1)+j);
        str = func2str(meth);
        plot(x,u(end,:));title(sprintf('%s - sigma = %f',str,sig));
    
    end
end
end

function u = upwind(t,sigma,N,m)
%Solve the prolem with upwind method
u = zeros(m,N+1);
u(1,1) = bcLeft(t(1));
for k = 1:m-1
    u(k+1,1) = bcLeft(t(k+1));
    u(k+1,2:end) = (1-sigma)*u(k,2:end) + sigma*u(k,1:end-1);
end
end

function u = laxFried(t,sigma,N,m)
%Solve the prolem with Lax-Friendrich
u = zeros(m,N+1);
u(1,1) = bcLeft(t(1));
for k = 1:m-1
    u(k+1,1) = bcLeft(t(k+1));
    u(k+1,2:end-1) = (u(k,3:end)+u(k,1:end-2))/2 - sigma*(u(k,3:end)-u(k,1:end-2))/2;
    u(k+1,end) = 2*u(k+1,end-1)-u(k+1,end-2);
end
end

function u = laxWen(t,sigma,N,m)
%Solve the prolem with Lax-Wendroff
u = zeros(m,N+1);
u(1,1) = bcLeft(t(1));
for k = 1:m-1
    u(k+1,1) = bcLeft(t(k+1));
    u(k+1,2:end-1) = u(k,2:end-1) - sigma*(u(k,3:end)-u(k,1:end-2))/2 + sigma^2*(u(k,3:end)-2*u(k,2:end-1)+u(k,1:end-2))/2;
    u(k+1,end) = 2*u(k+1,end-1)-u(k+1,end-2);
end
end

function [uLeft] = bcLeft(t)
%Gives the left boundary condition
T = 1;
uLeft = (-1).^floor(2.*t./T);
end

