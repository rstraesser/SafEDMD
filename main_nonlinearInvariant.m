clear;clc;clearvars;rng(0);close all;

%% System dynamics
lambda = 1; mu = -2;
f = @(x) [mu*x(1);lambda*(x(2)-x(1)^2)];
g = @(x) [0;1];
param.n=2; param.m=1;
ode = @(x,u) f(x)+g(x)*u;

%% Parameters
param.DeltaT = 0.01;
param.d=100;
param.cx=5e-3; 
param.cu=5e-3; 
param.delta = 0.05;
param.Rz=500; 
param.Qz='eye'; % 'eye' | 'optimize'
param.Tsim=10; 
param.xmax=1; param.xmin=-param.xmax;
param.umax=1; param.umin=-param.umax;
param.xplotmax=25; param.plotdist=0.1;
param.x0 = [ 10,-10,  8,- 8,  5,- 5,  5,  5,- 5,- 5;
             18, 18,  0,  0,-12,-12,  5,- 5,  5,- 5];
param.noise_level=0; 

%% Lifting
param.Phi = @(x) [1;x;x(2)-lambda/(lambda-2*mu)*x(1)^2];

%% Data generation
[X0,X1,U] = helperDataCollection(ode,param);

%% Apply SafEDMD
[param,sys,X,Y,compTimeSafEDMD] = helperSafEDMD(X0,X1,param);

%% Controller design
eps.P = 1e-6;
eps.F = 1e-6;
eps.Lambda = 1e-7;
eps.tau = 1e-7;
eps.nu = 1e-7;

[K,Kw,Pinv,sys,compTimeControllerDesign] = helperControllerDesign(sys,eps,param);

%% Computation times
fprintf('Computation time of SafEDMD: %fs\n',compTimeSafEDMD)
fprintf('Computation time of controller design: %fs\n',compTimeControllerDesign)

%% Validation
u = @(x) (eye(param.m)-Kw*kron(eye(param.m),param.hPhi(x))) \ (K*param.hPhi(x));

%% Comparison to EDMDc with LQR
fprintf('Solve EDMDc with LQR... ')
[KLQR,compTimeEDMDcLQR] = helperEDMDcLQR(X0,X1,U,param);
uLQR = @(x) -KLQR*param.hPhi(x);
fprintf('Computation time: %fs\n',compTimeEDMDcLQR)

%% Plotting
fprintf('Plotting...\n')
helperPlotting(ode,param,sys,Pinv,u,uLQR);