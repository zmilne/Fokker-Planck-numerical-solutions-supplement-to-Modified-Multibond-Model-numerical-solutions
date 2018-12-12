%Call the FokkerPlanck PDE solver and then calculate the average friction
%as a function of speed for a *harmonic* (spring) substrate potential
%Written by Zachary Banks Milne, University of Pennsylvania
%Copyright 2018, Zachary Banks Milne

%All uniits are standard metric kg, N, m, s, J, K

clc
clear all
close all
figure

%Dynamic and mechanical Parameters
V=1e-2; %the speed 
d=8e-18;% THE DIFFUSION CONSTANT
ks=1.3;%The substrate interaction spring constant
Ns=6e-6;%The substrate interaction damping constant
kc=1;%The cantilever spring constant
Nc=0;%The cantilever damping constant
Xc=2e-10;% the critical stretch length i.e., the maximum length any substrate interaction can stretch
tStep=3.00e-4/V*10e-9;% The time step is scaled by the speed
StdDist=.5e-11; %The standard deviation of the initial distriubtion

x = linspace(-Xc*1,2.5*Xc,500);%The initial spatial domain.The density of points (third entry) may be changed and should be made denser if P(X) ever becomes negative. 
t = 0:tStep:Xc*105/V;%The total time and time step are determined by the sliding speed V

PNEW=[];%DO NOT CHANGE THIS.
DoNewP=0;%DO NOT CHANGE THIS. It tells the program to initially use a normal distribution for P(X)

doplots=1;% use 1 if you want to see plots or 0 if not
FramesToPlot=1;%How many frames to plot e.g., 1 plots every frame, 2 plots every other frame etc.
dovids=0;%1=make a video, 0=don't make a video

if dovids%initialize video
    vidfile = VideoWriter('testmovie.mp4','MPEG-4');
    open(vidfile);
end


tLength=length(t);
TT=t;
qq=2;%qq is the index which tracks how many iterations of hte while loop below have been performed
XXnotBrokenindices=[];%Tracks the indices of the spatial domain points which contain interactions which have not yet left the boundaries (broken)
XnotBroken=x;%Initially, the spatial domain points which contain interactions which have not yet left the boundaries (broken, is the entire spatial domain)
MaxFf=1;FfNew=2;FfOld=1;

%The main while loop does the following: Starting with the first iteration,
%an initial probability distribution with standard deviation coming from
%StdDist is formed. Each loop solves the FP equation for three time steps.
%The next iteration uses the solution to the probability distribution from
%the third time step of the previous iteration as the first time step for
%its iteration. The spatial domain is redefined at every time step to set
%the boundaries x<=vt-Xc = x>=vt+X_c = 0. The mean friction force at each
%time step, <Ff(t)>=integral_{vt-Xc}^{vt+X_c}kc*P(X,t)dX, is solved at each
%iteration.

while (2*qq+1)<=length(TT)&&length(XnotBroken)>=4&&(FfNew>=5e-13||abs(FfNew-FfOld)>6e-15||qq<=30)
%Runs until the time is greater than the allowed time, until friciton falls
%below 5e-13 N (given at least 31 iterations), or the difference of the
%friction from the previous iteration is less than 6e-15 N (given at least
%31 iterations)
   
    TimeT=TT(2*qq-1:2*qq+1);%Sets the three time steps for this iterations 
    XnotBrokenHold=XnotBroken;
    XnotBroken=x(abs(V*TimeT(2)-x)<=Xc);
    if length(XnotBroken)>=4&&((x(end)-median(XnotBroken))<=3e-10||XnotBroken(end)<=V*TimeT(2))
        x=XnotBroken(1):(x(2)-x(1)):(x(end)+3e-10);
    end
%     length(XnotBroken);
    XXnotBrokenindicesHold=XXnotBrokenindices;
    XXnotBrokenindices=find(x<=max(XnotBroken)&x>=min(XnotBroken));
    if qq<=2%If this is the first two iterations, use a predefined normal distribution
        DoNewP=0;
        Probs(qq-1).P=FokkerPlanckPDEMovingBCHarmonicPotential(V,XnotBroken,TimeT,ks,kc,Ns,Nc,Xc,StdDist,d,DoNewP,PNEW); %This is where P(X) is solved for using FP
        Probs(qq-1).X=XnotBroken; %Save the spatial domain in the class 'Probs'
        MaxProb1=max(max(Probs(qq-1).P));
    else %If this is the 3rd or greater iteration, use the P(X) from the last time step of the previous iteration as the inital distribution for the current time step 
        DoNewP=1;%Tells FokkerPlanckPDEMovingBC that a new P(X) will be used 
        [xData, yData] = prepareCurveData(XnotBrokenHold, Probs(qq-2).P(end,:));%Prepares the new P(X) as initial distribution
        fittedmodel = fit(xData,yData,'linearinterp');%Makes the new P(X) as initial distribution
        Probs(qq-1).P=FokkerPlanckPDEMovingBCHarmonicPotential(V,XnotBroken,TimeT-tStep,ks,kc,Ns,Nc,Xc,StdDist,d,DoNewP,fittedmodel); %This is where P(X) is solved for using FP
        Probs(qq-1).X=XnotBroken; %Save the spatial domain in the class 'Probs'
%         
        if doplots %Make the plots
             if mod(qq,FramesToPlot)==0
                subplot(1,2,1)
%                 plot([1],[1])
                cla('reset')
                scatter(V*TimeT(2),0,'filled')
                box on
                hold on
                plot(fittedmodel,xData,Probs(qq-2).P(end,:))
                title('P(X) vs. X','fontsize',20)
                %This code has issues getting the axes label's to appear in plots (related to the 'cla' command), but they do
                %show up in videos. I haven't been able to figure this one
                %out but it's all good as long as they show up in videos.
        
                xlabel('X: mass position [m]','fontsize',20)
                ylabel('P(X)','fontsize',20)
                ylim([0 10e10])
                xlim([-2e-10 2e-9])
                b=gca;
                legend(b,'off');
            end
        end

        
    end
    FfOld=FfNew;
    %Calculate <Ff> for this iteration
    if size(Probs(qq-1).P,1)<=1
        Probs(qq-1).Ff=sum(kc*Probs(qq-1).X(1:end-1).*Probs(qq-1).P(1,1:end-1).*(diff(Probs(qq-1).X)));
    else
        Probs(qq-1).Ff=.5*sum(kc*Probs(qq-1).X(1:end-1).*Probs(qq-1).P(1,1:end-1).*(diff(Probs(qq-1).X)))+.5*sum(kc*Probs(qq-1).X(1:end-1).*Probs(qq-1).P(2,1:end-1).*(diff(Probs(qq-1).X)));
        
    end
  
    if doplots
        if mod(qq,FramesToPlot)==0
            subplot(1,2,2)
            box on
            hold on
            PlotFfVsTimeMovingBC%plots <Ff> for this iteration
            set(gcf,'Position',[200 200 1000 600])
        end
        if dovids %Make a video
            F(qq) = getframe(gcf);
            writeVideo(vidfile,F(qq));
        end
        
    end

    FfNew=Probs(qq-1).Ff;
    MaxFf=max(Probs(qq-1).Ff);
    qq=qq+1;
end
clc
Ftot=sum([Probs(:).Ff])*2*tStep/TimeT(end);%This calculates the avrage friction for the whole simulation

if doplots==0
    for i=1:qq-2
        PFf(i)=Probs(i).Ff;
    end
end

maxPFf=max(PFf); %Calculates the maximum Ff
C=[V d Ftot maxPFf];

if dovids
    close(vidfile)%com
end