function [KLQR,compTime] = helperEDMDcLQR(X0,X1,U,param)
compTime = tic;
%% EDMDc with LQR
    X=cell(param.m+1,1);Y=cell(param.m+1,1);
    %% Lift data to the Koopman space
    for k=0:param.m % each input u=0, u=e_1, ..., u=e_m
        X{k+1,1}=cell2mat(arrayfun(@(j) param.hPhi(X0{k+1,1}(:,j)), 1:param.d, 'UniformOutput', false));
        Y{k+1,1}=cell2mat(arrayfun(@(j) param.hPhi(X1{k+1,1}(:,j)), 1:param.d, 'UniformOutput', false));
    end

    %% Apply EDMDc
    X = horzcat(X{:});
    Y = horzcat(Y{:});
    U = horzcat(U{:});
    AB = Y*pinv([X;U]);
    A = AB(:,1:size(X,1));
    B = AB(:,size(X,1)+1:size(X,1)+param.m);
    
    %% Check stabilizability of (A,B0)
    stabilizable = 1;
    if rank(ctrb(A,B)) ~= param.N
        eigA = eig(A);
        eigA = eigA(abs(eigA)>=1);
        for i = 1:length(eigA)
            if rank([eigA(i)*eye(param.N) - A,B]) < param.N
                stabilizable = 0;
            end
        end
    end
    if ~stabilizable
        fprintf(2,'Linear Koopman system is NOT stabilizable.\n')
        KLQR = NaN(param.m,param.N);
    else
        if isfield(param,'LQR_Q')
            Q = param.LQR_Q;
        else
            Q = eye(param.N);
        end
        if isfield(param,'LQR_R')
            R = param.LQR_R;
        else
            R = eye(param.m);
        end
        KLQR = dlqr(A,B,Q,R);
    end
    compTime = toc(compTime);
end