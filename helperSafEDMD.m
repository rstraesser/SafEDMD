function [param,sys,X,Y,compTime] = helperSafEDMD(X0,X1,param)
%% Calculate data-driven approximation of the Koopman operator
compTime = tic;fprintf('Apply SafEDMD to calculate the bilinear surrogate model...')
param.N = size(param.Phi(zeros(param.n,1)),1)-1;
param.hPhi = @(x) [zeros(param.N,1),eye(param.N)]*param.Phi(x);
    X=cell(param.m+1,1);Y=cell(param.m+1,1);
    %% Lift data to the Koopman space
    for k=0:param.m % each input u=0, u=e_1, ..., u=e_m
        if k == 0
            X{k+1,1}=cell2mat(arrayfun(@(j) param.hPhi(X0{k+1,1}(:,j)), 1:param.d, 'UniformOutput', false));
        else
            X{k+1,1}=cell2mat(arrayfun(@(j) param.Phi(X0{k+1,1}(:,j)), 1:param.d, 'UniformOutput', false));
        end
        Y{k+1,1}=cell2mat(arrayfun(@(j) param.hPhi(X1{k+1,1}(:,j)), 1:param.d, 'UniformOutput', false));
    end

    % input u=0
    sys.A = Y{0+1,1}*pinv(X{0+1});
    sys.B0 = NaN(param.N,param.m);
    sys.tB = NaN(param.N,param.m*param.N);
    for i=1:param.m % each input u=e_1, ..., u=e_m
        [B0i_tBi] = Y{i+1,1}*pinv(X{i+1,1});   
        sys.B0(:,i) = B0i_tBi(:,1);
        sys.tB(:,(i-1)*param.N+1:i*param.N) = B0i_tBi(:,2:end) - sys.A;
    end
    sys.A(abs(sys.A)<1e-10) = 0;
    sys.B0(abs(sys.B0)<1e-10) = 0;
    sys.tB(abs(sys.tB)<1e-10) = 0;

    %% Check stabilizability of (A,B0)
    if rank(ctrb(sys.A,sys.B0)) ~= param.N
        eigA = eig(sys.A);
        eigA = eigA(abs(eigA)>=1);
        stabilizable = 1;
        for i = 1:length(eigA)
            if rank([eigA(i)*eye(param.N) - sys.A,sys.B0]) < param.N
                stabilizable = 0;
            end
        end
        if ~stabilizable
            error('Method not applicable since linear system part (A,B0) is not stabilizable!')
        end
    end
    compTime = toc(compTime);fprintf('Done.\n')
end