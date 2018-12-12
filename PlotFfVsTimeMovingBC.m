%Plots Friction force vs time from FokkerPlanckPDEMovingBC.m
%Written by Zachary Banks Milne, University of Pennsylvania
%Copyright 2018, Zachary Banks Milne

% clc
% % close all
% figure
% % 
% for i=1:qq-2
% PFf(i)=Probs(i).Ff;
% end
% scatter(V*t(1:2:2*length(PFf)),PFf)
%%
PFf(qq)=Probs(qq-1).Ff;
scatter(V*TimeT(2),PFf(qq),50,'b','filled')
ylim([0 4.6e-10])
xlim([-2e-10 2e-9])
title({'Mean F_f','vs. puller position, x'},'fontsize',20)
ylabel('Mean F_f [N]','fontsize',20)
xlabel('puller position, x [m]','fontsize',20)
