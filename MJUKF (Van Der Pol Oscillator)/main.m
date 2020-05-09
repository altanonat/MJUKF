%% Application of MJUKF for Van der Pol Oscillator
close all
clear all
clc
rng('default')
runtime = 100;

% Integration step
global deltat
deltat = 0.1;

% Number of iterations
ni = round(runtime/deltat);

% Number of states
global nse
ns  = 2;
nse = 2;
% number of estimated parameters
npar = 1;

% Initially states and their derivatives are set to zero
state(1:ns)           = zeros(1,ns); 
statederivative(1:ns) = zeros(1,ns);

%% Initial states
state(1) = 1.4;
state(2) = 0;
theta    = 0.2;
%%  time variables
global time 
time     = 0.0;
tplot(1) = time;
%% Parameters
parameters

P = Px;
statecap(1:nse)           = zeros(1,nse);
statecapderivative(1:nse) = zeros(1,nse);
%% Intial State Estimates
statecap(1) = 0;
statecap(2) = 5;
% Different initial state estimates
% statecap(1) = -1;
% statecap(2) = 3;

thetacapv = zeros(((2*ns)+1),npar);
% Initial parameter estimate
thetacap(1) = 5;
% thetacap(1) = -2;

for i=1:npar
    for j=1:((2*ns)+1)
        thetacapv(j,i)=(thetacap(i)-.5*(ns+1))+.5*j;
    end
end
thetacap = mean(thetacapv);

stateplot              = zeros(ns,ni);
statecapplot           = zeros(nse,ni);
thetacapplot           = zeros(npar,ni);
measurements           = zeros(2,ni);

% To measure average elapsed time
% REPS = 100;
REPS = 1;
tElapsed = zeros(REPS,1);
mmse     = zeros(REPS,ni);
tic;  % TIC, pair 1
for j=1:REPS
    tStart = tic;  % TIC, pair 2
    % Initial states
    state(1) = 1.4;
    state(2) = 0;
    theta    = 0.2;
    
    statecap(1) = 0;
    statecap(2) = 5;
    thetacapv   = zeros(((2*ns)+1),npar);
    thetacap(1) = 5;
    P = Px;
    for i=1:npar
        for k=1:((2*ns)+1)
            thetacapv(k,i)=(thetacap(i)-.5*(ns+1))+.5*k;
        end
    end
    thetacap = mean(thetacapv);
    time     = 0.0;
    tplot(1) = time;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for i=1:ni
        stateplot(:,i)    = state;
        thetacapplot(:,i) = thetacap;
        statecapplot(:,i) = statecap;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        statederivative  = mathematicalmodel(state,theta);
        state            = integration(state,statederivative, theta);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        noiseSigma           = 0.1 * state(1:ns);
        noise                = noiseSigma .* randn(1, length(state(1:ns)));
        measurements(1:ns,i) = state(1:ns)+noise;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [statecap,P,thetacap,thetacapv] = unscentedkalmanfilter(statecap,P,...
            measurements(:,i),thetacap,thetacapv);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        mse = rms(thetacapplot(1,i)-theta)^2;

        tplot(i) = time;
        mmse(j,i)   = mse;
    end
    tElapsed(j) = toc(tStart);  % TOC, pair 2  
end
averageTime    = mean(tElapsed);
minaverageTime = min(tElapsed);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%real \mu
realval = theta*ones(size(tplot));

% real \mu for step changes
% realval             = theta*ones(size(tplot));
% realval(1,1:401)    = 0.2;
% realval(1,401:1201) = 2;
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
xlabel('Time (s)','FontSize', 12)
ylabel('$x_{1}$','Interpreter','latex','FontSize', 12)
legend({'$x_{1}$','$\hat{x}_{1}$','$y_{1}$'},'Interpreter','latex','FontSize', 12)
axis([0 runtime -5 5])
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
xlabel('Time (s)','FontSize', 12)
ylabel('$x_{2}$','Interpreter','latex','FontSize', 12)
legend({'$x_{2}$','$\hat{x}_{2}$','$y_{2}$'},'Interpreter','latex','FontSize', 12)
axis([0 runtime -5 5])
set(gca,'fontsize',14)
set(gca,'XTick',(0:10:runtime))
set(gca,'YTick',(-5:1:5))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(3)
plot(tplot,thetacapplot,'--k','LineWidth',1);
hold on
plot(tplot,realval,'k','LineWidth',0.5,'DisplayName','Vector 1');
grid on
ax=gca;
ax.GridLineStyle = '--';
ax.GridColor = [0.6 0.6 0.6];
xlabel('Time (s)','FontSize', 12)
ylabel('$\hat{\mu}$','Interpreter','latex','FontSize', 12)
legend({'$\hat{\mu}$','$\mu$'},'FontSize', 12,'Interpreter','latex')
set(gca,'fontsize',14)
set(gca,'XTick',(0:10:runtime))