%FokkerPlanck PDE solver and then calculate the average friction
%as a function of speed with a periodic substrate potential
%Written by Zachary Banks Milne, University of Pennsylvania
%Copyright 2018, Zachary Banks Milne

function u=FokkerPlanckPDEMovingBC(v,x,t,ks,kc,Ns,Nc,Xc,STDDIST,DD)
global VV ns nc kB T ks kc nM D StdFokker doNewP PP  
StdFokker=STDDIST;
% PP=pp;
% doNewP=DoNewP;
VV=v;ns=Ns;nc=Nc;kB=1.38e-23;T=300;D=DD;
m = 0;%Do not change. This selects the correct PDE to solve.
sol = pdepe(m,@pdex1pde,@pdex1ic,@pdex1bc,x,t);%The main PDE solver
u = sol(:,:,1);

% --------------------------------------------------------------
function [c,f,s] = pdex1pde(x,t,u,DuDx)
global VV ns nc kB T ks kc nM D
c = ns+nc;
s = -(DuDx*(ns*VV+1e-9*cos((VV*t-x)*2*pi/2e-10)-x*(kc))-(kc-1e-9*2*pi/2e-10*sin((VV*t-x)*2*pi/2e-10))*u);
f = (ns+nc)*D*DuDx;
% --------------------------------------------------------------
function u0 = pdex1ic(x)
global VV ns nc kB T ks kc nM D StdFokker doNewP PP

% if doNewP==0
u0 = normpdf(x,0,StdFokker);
% else
% u0=PP(x);
% end

% --------------------------------------------------------------
function [pl,ql,pr,qr] = pdex1bc(xl,ul,xr,ur,t)
global VV ns nc kB T ks kc nM D

pl = ul;
ql = 0;
pr = ur;
qr = 0;