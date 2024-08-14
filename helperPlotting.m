function helperPlotting(ode,param,sys,Pinv,u,uLQR)
    if param.n == 2
        figure;grid on;hold all;
        title('Safe operating region for the nonlinear state $x$','interpreter','latex')
        % regions
        param.xplotmin = -param.xplotmax;
        xx = param.xplotmin:param.plotdist:param.xplotmax;
        yy = param.xplotmin:param.plotdist:param.xplotmax;
        [XX,YY] = meshgrid(xx,yy);
        x_test = [horzcat(XX(:))';horzcat(YY(:))'];
        N_test = size(x_test,2);
        xDelta = zeros(param.n,1);
        xSOR = zeros(param.n,1);
        for i = 1:N_test
            if param.hPhi(x_test(:,i))'*Pinv*param.hPhi(x_test(:,i)) <= 1
                xSOR= [xSOR,x_test(:,i)];
            end
            if param.hPhi(x_test(:,i))'*(-sys.Pi.Qz)*param.hPhi(x_test(:,i)) <= sys.Pi.Rz
                xDelta = [xDelta,x_test(:,i)];
            end
        end
        pDelta=plot(xDelta(1,:),xDelta(2,:),'.b');
        boundary_idx = boundary(xDelta(1,:)',xDelta(2,:)',0.95);
        xBoundary = xDelta(:,boundary_idx);
        plot(xBoundary(1,:),xBoundary(2,:),'b')
        %
        pSOR=plot(xSOR(1,:),xSOR(2,:),'.r');
        boundary_idx = boundary(xSOR(1,:)',xSOR(2,:)',0.95);
        xBoundary = xSOR(:,boundary_idx);
        plot(xBoundary(1,:),xBoundary(2,:),'r')
        
        % closed-loop trajectories
        tSteps = ceil(param.Tsim/param.DeltaT);  
        for k=1:size(param.x0,2)
            pX0=plot(param.x0(1,k),param.x0(2,k),'ko');
            xsim = param.x0(:,k);
            for tt = 1:tSteps
                [~,xsimTT] = ode45(@(t,x) ode(x,u(xsim(:,end))),[(tt-1)*param.DeltaT,tt*param.DeltaT],xsim(:,end)); 
                xsim = [xsim,xsimTT(2:end,:)'];
            end
            pSafEDMD=plot(xsim(1,:),xsim(2,:),'k');
            
            % LQR
            xsim = param.x0(:,k);
            for tt = 1:tSteps
                [~,xsimTT] = ode45(@(t,x) ode(x,uLQR(xsim(:,end))),[(tt-1)*param.DeltaT,tt*param.DeltaT],xsim(:,end)); 
                xsim = [xsim,xsimTT(2:end,:)'];
            end
            pLQR=plot(xsim(1,:),xsim(2,:),'k--');
        end
        xlim([param.xplotmin,param.xplotmax])
        ylim([param.xplotmin,param.xplotmax])
        xlabel('$x_1$','interpreter','latex')
        ylabel('$x_2$','interpreter','latex')
        legend([pDelta,pSOR,pX0,pSafEDMD,pLQR],{'$\mathbf{\Delta}_\Phi$','SOR','$x_0$','SafEDMD','LQR'},'location','best','interpreter','latex')
    else
        warning('Plotting only implemented for n=2.')
end