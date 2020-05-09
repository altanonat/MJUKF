% Application of MJUKF to Lorenz System
% Simulation run time
close all
clear all
clc
rng('default');
runtime = 100;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Integration step
global deltat
deltat = 0.01;
% Number of iterations
ni = round(runtime/deltat);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Number of states
global nse
ns  = 3;
nse = 3;
% number of estimated parameters
npar = 3;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initially states and their derivatives are set to zero
state(1:ns)           = zeros(1,ns); 
statederivative(1:ns) = zeros(1,ns);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial states
state(1) = 0.9;
state(2) = 1;
state(3) = 1.1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% time variables
global time 
time     = 0.0;
tplot(1) = time;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters
parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
P = Px;
statecap(1:nse)           = zeros(1,nse);
statecapderivative(1:nse) = zeros(1,nse);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Different initial state estimates
% statecap(1) = 1.5;
% statecap(2) = 1.5;
% statecap(3) = 1.5;
% statecap(1) = 0.5;
% statecap(2) = 0.5;
% statecap(3) = 0.5;
% statecap(1) = 0.1;
% statecap(2) = 0.1;
% statecap(3) = 0.1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
theta(1) = 10;
theta(2) = 28;
theta(3) = 8/3;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
thetacapv = zeros(((2*ns)+1),npar);
% Different initial parameter estimates
% thetacap(1) = 5;
% thetacap(2) = 21;
% thetacap(3) = 1/3;
% thetacap(1) = 10;
% thetacap(2) = 28;
% thetacap(3) = 8/3;
% thetacap(1) = 13;
% thetacap(2) = 33;
% thetacap(3) = 10/3;

% for i=1:npar
%     for j=1:((2*ns)+1)
%         thetacapv(j,i)=(thetacap(i)-.5*(ns+1))+.5*j;
%     end
% end
% thetacap = mean(thetacapv);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
stateplot    = zeros(ns,ni);
statecapplot = zeros(nse,ni);
thetacapplot = zeros(nse,ni);
measurements = zeros(3,ni);
% To measure average elapsed time
REPS = 1;
% REPS = 100;
tElapsed = zeros(REPS,1);
mmse     = zeros(REPS,ni,3);
for j=1:REPS
    tStart = tic;  % TIC, pair 2
    % Initial States
    state(1) = 0.9;
    state(2) = 1;
    state(3) = 1.1;

    % Initial states and parameter estimates
    statecap(1) = 1.5;
    statecap(2) = 1.5;
    statecap(3) = 1.5;
    thetacap(1) = 5;
    thetacap(2) = 21;
    thetacap(3) = 1/3;
    
    for i=1:npar
        for k=1:((2*ns)+1)
            thetacapv(k,i)=(thetacap(i)-.5*(ns+1))+.5*k;
        end
    end
    thetacap = mean(thetacapv);
    
    P = Px;
    
    time         = 0.0;
    tcntr        = 1;
    tplot(tcntr) = time;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for i=1:ni
        stateplot(:,i)    = state;
        thetacapplot(:,i) = thetacap;
        statecapplot(:,i) = statecap;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        statederivative  = mathematicalmodel(state,theta);
        state            = integration(state,statederivative, theta);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % 10% SNR
        noiseSigma          = 0.1 * state(1:3);
        noise               = noiseSigma .* randn(1, length(state(1:3)));
        measurements(1:3,i) = state(1:3)+noise;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % UKF
        [statecap,P,thetacap,thetacapv] = unscentedkalmanfilter(statecap,P,...
            measurements(:,i),thetacapv,thetacap);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        mse(1) = rms(thetacapplot(1,i)-theta(1))^2;
        mse(2) = rms(thetacapplot(2,i)-theta(2))^2;
        mse(3) = rms(thetacapplot(3,i)-theta(3))^2;

        mmse(j,i,:)  = mse;
        tplot(i) = time;
    end
    tElapsed(j) = toc(tStart);  % TOC, pair 2
end
averageTime    = mean(tElapsed);
minaverageTime = min(tElapsed);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
thetaplot(1,:) = 10*ones(1,ni);
thetaplot(2,:) = 28*ones(1,ni);
thetaplot(3,:) = (8/3)*ones(1,ni);

% real \mu for step changes
% thetaplot(1,:)          = theta(1)*ones(size(tplot));
% thetaplot(1,4001:12001) = 7;
% %real \mu for changes
% thetaplot(2,:)          = theta(2)*ones(size(tplot));
% thetaplot(2,4001:12001) = 21;
% %real \mu for changes
% thetaplot(3,:)          = theta(3)*ones(size(tplot));
% thetaplot(3,4001:12001) = 1;

figure(1)
plot(tplot,thetaplot(1,:),'k','LineWidth',1)
hold on
plot(tplot,thetacapplot(1,:),'-.k','LineWidth',0.3)
grid on
ax=gca;
ax.GridLineStyle = '--';
ax.GridColor = [0.6 0.6 0.6];
xlabel('Time (s)','FontSize', 14)
ylabel('$\sigma$','Interpreter','latex','FontSize', 14)
legend({'$\sigma$','$\hat{\sigma}$'},'Interpreter','latex','FontSize', 14)
axis([0 runtime -5 15])
set(gca,'fontsize',14)
set(gca,'XTick',(0:10:runtime))
set(gca,'YTick',(-5:5:15))
figure(2)
plot(tplot,thetaplot(2,:),'k','LineWidth',1)
hold on
plot(tplot,thetacapplot(2,:),'-.k','LineWidth',0.3)
grid on
ax=gca;
ax.GridLineStyle = '--';
ax.GridColor = [0.6 0.6 0.6];
xlabel('Time (s)','FontSize', 14)
ylabel('$\rho$','Interpreter','latex','FontSize', 14)
legend({'$\rho$','$\hat{\rho}$'},'Interpreter','latex','FontSize', 14)
axis([0 time 0 35])
set(gca,'fontsize',14)
set(gca,'XTick',(0:10:runtime))
set(gca,'YTick',(0:7:35))
figure(3)
plot(tplot,thetaplot(3,:),'k','LineWidth',1)
hold on
plot(tplot,thetacapplot(3,:),'-.k','LineWidth',0.3)
grid on
ax=gca;
ax.GridLineStyle = '--';
ax.GridColor = [0.6 0.6 0.6];
xlabel('Time (s)','FontSize', 14)
ylabel('$\beta$','Interpreter','latex','FontSize', 14)
legend({'$\beta$','$\hat{\beta}$'},'Interpreter','latex','FontSize', 14)
axis([0 runtime 0 4])
set(gca,'fontsize',14)
set(gca,'XTick',(0:10:runtime))
set(gca,'YTick',(0:1:4))

figure(4)
plot(tplot,stateplot(1,:),'k','LineWidth',1)
hold on
plot(tplot,statecapplot(1,:),'-.k','LineWidth',0.3)
hold on
plot(tplot,measurements(1,:),'--k','LineWidth',0.3)
grid on
ax=gca;
ax.GridLineStyle = '--';
ax.GridColor = [0.6 0.6 0.6];
xlabel('Time (s)','FontSize', 14)
ylabel('$x_{1},\hat{x}_{1},y_{1}$','Interpreter','latex','FontSize', 14)
legend({'$x_{1}$','$\hat{x}_{1}$','$y_{1}$'},'Interpreter','latex','FontSize', 14)
axis([0 runtime -inf inf])
set(gca,'fontsize',14)
set(gca,'XTick',(0:10:runtime))

figure(5)
plot(tplot,stateplot(2,:),'k','LineWidth',1)
hold on
plot(tplot,statecapplot(2,:),'-.k','LineWidth',0.3)
hold on
plot(tplot,measurements(2,:),'--k','LineWidth',0.3)
grid on
ax=gca;
ax.GridLineStyle = '--';
ax.GridColor = [0.6 0.6 0.6];
xlabel('Time (s)','FontSize', 14)
ylabel('$x_{2},\hat{x}_{2},y_{2}$','Interpreter','latex','FontSize', 14)
legend({'$x_{2}$','$\hat{x}_{2}$','$y_{2}$'},'Interpreter','latex','FontSize', 14)
axis([0 runtime -inf inf])
set(gca,'fontsize',14)
set(gca,'XTick',(0:10:runtime))

figure(6)
plot(tplot,stateplot(3,:),'k','LineWidth',1)
hold on
plot(tplot,statecapplot(3,:),'-.k','LineWidth',0.3)
hold on
plot(tplot,measurements(3,:),'--k','LineWidth',0.3)
grid on
ax=gca;
ax.GridLineStyle = '--';
ax.GridColor = [0.6 0.6 0.6];
xlabel('Time (s)','FontSize', 14)
ylabel('$x_{3},\hat{x}_{3},y_{3}$','Interpreter','latex','FontSize', 14)
legend({'$x_{3}$','$\hat{x}_{3}$','$y_{3}$'},'Interpreter','latex','FontSize', 14)
axis([0 runtime -inf inf])
set(gca,'fontsize',14)
set(gca,'XTick',(0:10:runtime))

figure(7)
plot(tplot,stateplot(1,:),'k','LineWidth',1)
hold on
plot(tplot,statecapplot(1,:),'-.k','LineWidth',0.3)
hold on
plot(tplot,measurements(1,:),'--k','LineWidth',0.3)
grid on
ax=gca;
ax.GridLineStyle = '--';
ax.GridColor = [0.6 0.6 0.6];
xlabel('Time (s)','FontSize', 14)
ylabel('$x_{1},\hat{x}_{1},y_{1}$','Interpreter','latex','FontSize', 14)
legend({'$x_{1}$','$\hat{x}_{1}$','$y_{1}$'},'Interpreter','latex','FontSize', 14)
axis([0 8 -inf inf])
set(gca,'fontsize',14)
set(gca,'XTick',(0:1:8))

figure(8)
plot(tplot,stateplot(2,:),'k','LineWidth',1)
hold on
plot(tplot,statecapplot(2,:),'-.k','LineWidth',0.3)
hold on
plot(tplot,measurements(2,:),'--k','LineWidth',0.3)
grid on
ax=gca;
ax.GridLineStyle = '--';
ax.GridColor = [0.6 0.6 0.6];
xlabel('Time (s)','FontSize', 14)
ylabel('$x_{2},\hat{x}_{2},y_{2}$','Interpreter','latex','FontSize', 14)
legend({'$x_{2}$','$\hat{x}_{2}$','$y_{2}$'},'Interpreter','latex','FontSize', 14)
axis([0 8 -inf inf])
set(gca,'fontsize',14)
set(gca,'XTick',(0:1:8))

figure(9)
plot(tplot,stateplot(3,:),'k','LineWidth',1)
hold on
plot(tplot,statecapplot(3,:),'-.k','LineWidth',0.3)
hold on
plot(tplot,measurements(3,:),'--k','LineWidth',0.3)
grid on
ax=gca;
ax.GridLineStyle = '--';
ax.GridColor = [0.6 0.6 0.6];
xlabel('Time (s)','FontSize', 14)
ylabel('$x_{3},\hat{x}_{3},y_{3}$','Interpreter','latex','FontSize', 14)
legend({'$x_{3}$','$\hat{x}_{3}$','$y_{3}$'},'Interpreter','latex','FontSize', 14)
axis([0 8 -inf inf])
set(gca,'fontsize',14)
set(gca,'XTick',(0:1:8))