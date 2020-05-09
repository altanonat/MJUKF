%% Application of JUKF to Van der Pol Oscillator
%% SNR is 10%
% Simulation run time
close all
clear all
clc
runtime = 100;

% Integration step
global deltat
deltat = 0.1;

% Number of iterations
ni = round(runtime/deltat);

% Number of states
global nse
ns  = 2;
nse = 3;
% number of estimated parameters
global npar
npar = 1;

% Initially states and their derivatives are set to zero
state(1:ns)           = zeros(1,ns); 
statederivative(1:ns) = zeros(1,ns);

% Initial states
state(1) = 1.4;
state(2) = 0;
mu       = 0.2;
% time variables
global time 
time     = 0.0;
tplot(1) = time;

% Parameters
parameters

P = Px;
statecap(1:nse)           = zeros(1,nse);
statecapderivative(1:nse) = zeros(1,nse);

statecap(1) = 0;
statecap(2) = 5;
statecap(3) = 5;

stateplot    = zeros(ns,ni);
statecapplot = zeros(nse,ni);
measurements = zeros(2,ni);
% To measure elapsed time
% REPS = 1000;
REPS = 1;
tElapsed = zeros(REPS,1);
mmse     = zeros(REPS,ni);
for j=1:REPS
    tStart = tic;  % TIC, pair 2
    % Initial states
    state(1) = 1.4;
    state(2) = 0;
    mu       = 0.2;
    
    statecap(1) = 0;
    statecap(2) = 5;
    statecap(3) = 5;
    
    time         = 0.0;
    tcntr        = 1;
    tplot(tcntr) = time;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for i=1:ni
        stateplot(:,i)    = state;
        statecapplot(:,i) = statecap;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        statederivative  = mathematicalmodel(state,mu);
        state            = integration(state,statederivative, mu);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        noiseSigma           = 0.1 * state(1:ns);
        noise                = noiseSigma .* randn(1, length(state(1:ns)));
        measurements(1:ns,i) = state(1:ns)+noise;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [statecap,P] = unscentedkalmanfilter(statecap,P,measurements(:,i));
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        mse = rms(statecapplot(3,i)-mu)^2;
        mmse(j,i) = mse;
        tplot(i)  = time;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tElapsed(j) = toc(tStart);  % TOC, pair 2
end
averageTime    = mean(tElapsed);
minaverageTime = min(tElapsed);
%% Variables to be plotted
figure(1)
plot(tplot,stateplot(1,:),'k','LineWidth',0.5)
hold on
plot(tplot,statecapplot(1,:),'--k','LineWidth',1)
hold on
plot(tplot,measurements(1,:),'LineWidth',0.5,'color',[.8,.8,.8])
grid on
ax=gca;
ax.GridLineStyle = '--';
ax.GridColor = [0.6 0.6 0.6];
xlabel('Time (s)','FontSize', 14)
ylabel('$x_{1}$','Interpreter','latex','FontSize', 14)
legend({'$x_{1}$','$\hat{x}_{1}$','$y_{1}$'},'Interpreter','latex','FontSize', 14)
axis([0 runtime -inf inf])
set(gca,'fontsize',14)
set(gca,'XTick',(0:10:runtime))
set(gca,'YTick',(-5:1:5))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(2)
plot(tplot,stateplot(2,:),'k','LineWidth',0.5)
hold on
plot(tplot,statecapplot(2,:),'--k','LineWidth',1)
hold on
plot(tplot,measurements(2,:),'LineWidth',0.5,'color',[.8,.8,.8])
grid on
ax=gca;
ax.GridLineStyle = '--';
ax.GridColor = [0.6 0.6 0.6];
xlabel('Time (s)','FontSize', 14)
ylabel('$x_{2}$','Interpreter','latex','FontSize', 14)
legend({'$x_{2}$','$\hat{x}_{2}$','$y_{2}$'},'Interpreter','latex','FontSize', 14)
axis([0 runtime -inf inf])
set(gca,'fontsize',14)
set(gca,'XTick',(0:10:runtime))
set(gca,'YTick',(-5:1:5))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(3)
plot(tplot,statecapplot(3,:),'--k','LineWidth',1)
hold on
realval = mu*ones(size(tplot));
plot(tplot,realval,'k','LineWidth',0.5)
grid on
ax=gca;
ax.GridLineStyle = '--';
ax.GridColor = [0.6 0.6 0.6];
xlabel('Time (s)','FontSize', 14)
ylabel('$\hat{\mu}$','Interpreter','latex','FontSize', 14)
legend({'$\hat{\mu}$','$\mu$'},'Interpreter','latex','FontSize', 14)
axis([0 runtime -inf inf])
set(gca,'fontsize',14)
set(gca,'XTick',(0:10:runtime))