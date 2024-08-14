function QzHeuristics = helperHeuristicalOptimizationDeltaPhi(sys,eps,param)
% Optimization: Decision variables
lmis = [];
    % (Inverse) Lyapunov matrix
    P = sdpvar(param.N,param.N,'symmetric');
    lmis = [lmis,P >= eps.P*eye(param.N)];
    % State-feedback gain
    L = sdpvar(param.m,param.N,'full');
    % Full-information feedback gain
    Lw = sdpvar(param.m,param.m*param.N,'full');

    % multiplier bilinearity
    sys.Pi.Qz = -eye(param.N);
    Piinv = inv([sys.Pi.Qz,sys.Pi.Sz;sys.Pi.Sz',sys.Pi.Rz]);
    tQz = Piinv(1:param.N,1:param.N);
    tSz = Piinv(1:param.N,param.N+1:end);
    tRz = Piinv(param.N+1:end,param.N+1:end);
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
    
    % Solve optimization
    opt = optimize(lmis,[],sdpsettings('solver','mosek'));

    if opt.problem ~= 0
        warning('Pre-computing Qz not successful! Using Qz=-I instead.')
        QzHeuristics = -eye(param.N);
    else
        Pinv = double(P) \ eye(param.N);
        QzHeuristics = -Pinv/norm(Pinv);
    end
end