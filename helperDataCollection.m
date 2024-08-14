function [X0,X1,U] = helperDataCollection(ode,param)
%% Collect data of nonlinear system
randuni = @(dim) (2*rand(dim(1),dim(2))-1);
xwidth = (param.xmax - param.xmin)/2;
xcenter = param.xmax - xwidth;
    X0=cell(param.m+1,1);U=cell(param.m+1,1);X1=cell(param.m+1,1);
    I_m = eye(param.m);
    for k=0:param.m % each input u=0, u=e_1, ..., u=e_m
        X0{k+1,1} = xcenter + xwidth.*randuni([param.n,param.d]);
        if k == 0
            U{k+1,1} = zeros(param.m,param.d);
        else
            U{k+1,1} = repmat(I_m(:,k),[1,param.d]);
        end
        X1{k+1,1}=NaN(param.n,param.d);
        for j=1:param.d % uniformly sampled data
            [~,xnext] =  ode45(@(t,x) ode(x,U{k+1,1}(:,j)),[0,param.DeltaT],X0{k+1,1}(:,j)); 
            X1{k+1,1}(:,j) = xnext(end,:)';
        end
    end
    for k=0:param.m % each input u=0, u=e_1, ..., u=e_m
        W{k+1,1} = param.noise_level*randuni([param.n,param.d]);
        for j=1:param.d % uniformly sampled data
            X1{k+1,1}(:,j) = X1{k+1,1}(:,j) + W{k+1,1}(:,j);
        end
    end
end