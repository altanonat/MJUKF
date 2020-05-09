%% Application of JUKF to Lorenz System
%% SNR is 10%
% Simulation run time
close all
clear all
clc
rng('default')
runtime = 100;

% Integration step
global deltat
deltat = 0.01;

% Number of iterations
ni = round(runtime/deltat);

% Number of states
global nse
ns  = 3;
% Number of estimated states (with parameters)
nse = 6;

% Initially states and their derivatives are set to zero
state(1:ns)           = zeros(1,ns); 
statederivative(1:ns) = zeros(1,ns);

% Initial states
state(1) = 0.9;
state(2) = 1;
state(3) = 1.1;

% time variables
global time 
time         = 0.0;
tcntr        = 1;
tplot(tcntr) = time;

% Parameters
parameters

P = Px;
statecap(1:nse)           = zeros(1,nse);
statecapderivative(1:nse) = zeros(1,nse);

statecap(1) = 1.5;
statecap(2) = 1.5;
statecap(3) = 1.5;
statecap(4) = 5;
statecap(5) = 21;
statecap(6) = 1/3;

global sigma
global rho
global zheta
stateplot              = zeros(ns,ni);
statederivativeplot    = zeros(ns,ni);
statecapplot           = zeros(nse,ni);
statecapderivativeplot = zeros(nse,ni);
measurements           = zeros(3,ni);
% To measure average elapsed time
REPS = 1;
tElapsed = zeros(REPS,1);
mmse     = zeros(REPS,ni,3);
for j=1:REPS
    tStart = tic;  % TIC, pair 2
    % Initial states
    state(1) = 0.9;
    state(2) = 1;
    state(3) = 1.1;
    
    P = Px;
    
    statecap(1) = 1.5;
    statecap(2) = 1.5;
    statecap(3) = 1.5;
    statecap(4) = 5;
    statecap(5) = 21;
    statecap(6) = 1/3;
    % Different initial state estimates
    %     statecap(1) = 0.9;
    %     statecap(2) = 1;
    %     statecap(3) = 1.1;
    %     statecap(4) = 20;
    %     statecap(5) = 10;
    %     statecap(6) = 1/10;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for i=1:ni
        stateplot(:,i)    = state;
        statecapplot(:,i) = statecap;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        statederivative  = mathematicalmodel(state);
        state            = integration(state,statederivative);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % GR Noise added to states measurements
        noiseSigma          = 0.1 * state(1:3);
        noise               = noiseSigma .* randn(1, length(state(1:3)));
        measurements(1:3,i) = state(1:3)+noise;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % UKF
        [statecap,P] = unscentedkalmanfilter(statecap,P,measurements(:,i));
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        mse(1) = rms(statecapplot(4,i)-sigma)^2;
        mse(2) = rms(statecapplot(5,i)-rho)^2;
        mse(3) = rms(statecapplot(6,i)-zheta)^2;

        mmse(j,i,:) = mse;
        tplot(i)    = time;
    end
    tElapsed(j) = toc(tStart);  % TOC, pair 2
end
averageTime    = mean(tElapsed);
minaverageTime = min(tElapsed);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Variables to be plotted
realsigma = 10*ones(size(tplot));
figure(1)
plot(tplot,realsigma,'k','LineWidth',1)
hold on
plot(tplot,statecapplot(4,:),'-.k','LineWidth',0.3)
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
realrho = 28*ones(size(tplot));
figure(2)
plot(tplot,realrho,'k','LineWidth',1)
hold on
plot(tplot,statecapplot(5,:),'-.k','LineWidth',0.3)
grid on
ax=gca;
ax.GridLineStyle = '--';
ax.GridColor = [0.6 0.6 0.6];
xlabel('Time (s)','FontSize', 14)
ylabel('$\rho$','Interpreter','latex','FontSize', 14)
legend({'$\rho$','$\hat{\rho}$'},'Interpreter','latex','FontSize', 14)
axis([0 runtime 0 35])
set(gca,'fontsize',14)
set(gca,'XTick',(0:10:runtime))
set(gca,'YTick',(0:7:35))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
realbeta = 8/3*ones(size(tplot));
figure(3)
plot(tplot,realbeta,'k','LineWidth',2)
hold on
plot(tplot,statecapplot(6,:),'-.k','LineWidth',0.3)
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