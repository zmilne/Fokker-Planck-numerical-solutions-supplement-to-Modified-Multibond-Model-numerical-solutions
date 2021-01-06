%Call the FokkerPlanck PDE solver and then calculate the average friction
%as a function of speed for a *periodic* (cosine) substrate potential
%Written by Zachary Banks Milne, University of Pennsylvania
%Copyright 2018, Zachary Banks Milne

%All uniits are standard metric kg, N, m, s, J, K
clc
clear all
close all
figure
%Dynamic and mechanical Parameters
V=1e-7; %the speed 
d=9e-17;% THE DIFFUSION CONSTANT
ks=1.3;%The substrate interaction spring constant
Ns=6e-6;%The substrate interaction damping constant
kc=1;%The cantilever spring constant
Nc=0;%The cantilever damping constant
Xc=2e-10;% the critical stretch length i.e., the maximum length any substrate interaction can stretch
tStep=3.00e-4/V*10e-9;% The time step is scaled by the speed
StdDist=1e-11; %The standard deviation of the initial distriubtion
clear P t
x = linspace(-Xc*4,4*Xc,500);
t = 0:tStep:Xc*10/V;
tLength=length(t);


        Probs.P=FokkerPlanckPDEMovingBCPeriodicPotential(V,x,t,ks,kc,Ns,Nc,Xc,StdDist,d); %This is where P(X) is solved for using FP. This is much more straightforward than having a harmonic potential with moving boundary conditions
        Probs.X=x;%Save the spatial domain in the class 'Probs'
        MaxProb1=max(max(Probs.P));
        for i=1:length(t)
        Probs.Ff(i)=sum(kc*Probs.X(1:end-1).*Probs.P(i,1:end-1).*(diff(Probs.X)));%Solve for the average friction at each time step
        end
 plot(t,Probs.Ff)
 maxFf=max(Probs.Ff)
 PlotAllProbs