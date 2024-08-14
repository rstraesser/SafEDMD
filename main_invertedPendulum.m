clear;clc;clearvars;rng(0);close all;

%% System dynamics
m   = 1;
l   = 1;
b   = 0.5;
g   = 9.81;
f   = @(x) [x(2,:); g/l*sin(x(1,:)) - b/m/l^2*x(2,:)];
g    = @(x) [0; 1/m/l^2];
param.n=2; param.m=1;
ode = @(x,u) f(x)+g(x)*u;

%% Parameters
param.DeltaT = 0.01;
param.d = 6000;
param.cx=3e-4; 
param.cu=3e-4; 
param.delta = 0.05;
param.Rz=2.5;
param.Qz='optimize'; % 'eye' | 'optimize'
param.Tsim=10; 
param.xmax=10; param.xmin=-2;
param.umax=10; param.umin=-param.umax;
param.xplotmax = 10; param.plotdist=0.05;
param.x0 = [ 3.5, 2,-3.5,-2, 3,-3;
            -2  ,-7, 2  , 7,-6, 6];
param.noise_level=1e-2;

%% Lifting
param.Phi = @(x) [1;x;sin(x(1))];

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
[KLQR,compTimeEDMDcLQR] = helperEDMDcLQR(X0,X1,U,param);
uLQR = @(x) -KLQR*param.hPhi(x);

%% Plotting
fprintf('Plotting...\n')
helperPlotting(ode,param,sys,Pinv,u,uLQR);