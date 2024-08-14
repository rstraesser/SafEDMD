function [K,Kw,Pinv,sys,compTime] = helperControllerDesign(sys,eps,param)
fprintf('Design a stabilizing controller...')
% Optimization: Decision variables
lmis = [];
    % (Inverse) Lyapunov matrix
    P = sdpvar(param.N,param.N,'symmetric');
    lmis = [lmis,P >= eps.P*eye(param.N)];
    % State-feedback gain
    L = sdpvar(param.m,param.N,'full');
    % Full-information feedback gain
    Lw = sdpvar(param.m,param.m*param.N,'full');

    % region overapproximating the bilinearity
    if isfield(param,'Sz')
        sys.Pi.Sz = param.Sz; 
    else
        sys.Pi.Sz = zeros(param.N,1);
    end
    if isfield(param,'Rz')
        sys.Pi.Rz = param.Rz; 
    else
        sys.Pi.Rz = 1;
    end
    switch param.Qz
        case 'eye'
            sys.Pi.Qz = -eye(param.N);
        case 'optimize'
            [~,sys.Pi.Qz] = evalc("helperHeuristicalOptimizationDeltaPhi(sys,eps,param)");
        otherwise
            error("Please specify 'param.Qz' as either 'eye' or 'optimize'!")
    end
    % multiplier bilinearity
    if max(max(sys.Pi.Qz - sys.Pi.Qz')) < 1e-10 && max(max(sys.Pi.Rz - sys.Pi.Rz')) < 1e-10
        Piinv = inv([sys.Pi.Qz,sys.Pi.Sz;sys.Pi.Sz',sys.Pi.Rz]);
        if max(max(Piinv - Piinv')) < 1e-10
            tQz = Piinv(1:param.N,1:param.N);
            tSz = Piinv(1:param.N,param.N+1:end);
            tRz = Piinv(param.N+1:end,param.N+1:end);
        else
            error('Matrix Pi^{-1} is not symmetric!')
        end
    else
        error('Matrix Pi is not symmetric!')
    end
    Lambda = sdpvar(param.m,param.m,'symmetric');
    lmis = [lmis,Lambda - eps.Lambda*eye(param.m) >= 0];
    
    if param.cx < 0 || param.cu < 0
        error('The constants for the proportional bound need to be positive!')
    elseif param.cx == 0 && param.cu == 0
        tau = 0;
        F13 = []; F23 = []; F33 = []; F34 = []; F35 = [];
    else
        % multiplier Koopman approximation error
        tau = sdpvar(1);
        lmis = [lmis,tau >= eps.tau];
        
        F13 = zeros(param.N,param.N+param.m);
        F23 = -kron(eye(param.m),tSz')*[zeros(param.N*param.m,param.N),Lw'];
        F33 = tau*0.5*blkdiag(1/param.cx^2*eye(param.N),1/param.cu^2*eye(param.m));
        F34 = [P;L];
        F35 = -[zeros(param.N,param.N*param.m);Lw];
    end
    
    F11 = P - tau*eye(param.N);
    F12 = -sys.tB*kron(Lambda,tSz) - sys.B0*Lw*kron(eye(param.m),tSz);
    F14 = sys.A*P + sys.B0*L;
    F15 = sys.tB*kron(Lambda,eye(param.N)) + sys.B0*Lw;
    F22 = kron(Lambda,tRz) - Lw*kron(eye(param.m),tSz) - (Lw*kron(eye(param.m),tSz))';
    F24 = L;
    F25 = Lw;
    F44 = P;
    F45 = zeros(param.N,param.N*param.m);
    F55 = -kron(Lambda,inv(tQz));

    F = [F11 ,F12 ,F13 ,F14 ,F15 ;
         F12',F22 ,F23 ,F24 ,F25 ;
         F13',F23',F33 ,F34 ,F35 ;
         F14',F24',F34',F44 ,F45 ;
         F15',F25',F35',F45',F55];
    
    lmis = [lmis,F - eps.F*eye(size(F)) >= 0];
    
    % invariance of safe operating region
    nu = sdpvar(1);
    lmis = [lmis,nu >= eps.nu];
    FI11 = nu*tRz - 1;
    FI12 = -nu*tSz';
    FI22 = nu*tQz + P;
    FI = [ FI11 ,FI12 ;
           FI12',FI22];
    lmis = [lmis,FI <= 0];
    
    % Solve optimization
    cost = -trace(P);
    opt = optimize(lmis,cost,sdpsettings('solver','mosek'));

    % Store obtained decision variables
    P = double(P);
    Pinv = P \ eye(param.N);
    Lambda = value(Lambda);

    K = double(L) / P;
    Kw = double(Lw)*kron(inv(Lambda),eye(param.N));
    compTime = opt.solvertime;
    switch opt.problem 
        case -1
            error('Controller design was not successful: %s: ILL_POSED', opt.info)
        case 0 
            fprintf('Controller design completed. %s: PRIMAL_AND_DUAL_FEASIBLE\n', opt.info)
        case 1
            error('Controller design was not successful: %s: DUAL_INFEASIBLE\n', opt.info)
        case 2
            error('Controller design was not successful: %s: PRIMAL_INFEASIBLE\n', opt.info)
        case 4
            error('Controller design was not successful: %s: UNKNOWN\n', opt.info)
        otherwise
            error('Controller design was not successful: %s\n', opt.info)
    end
end