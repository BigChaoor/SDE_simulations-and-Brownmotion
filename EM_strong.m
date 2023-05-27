% EMSTRONG Test strong convergence of Euler-Maruyama
% 
% Solves  dX = lambda*X dt + mu*X dW, X(0) = X0
%       where lambda = 2, mu = 1 and X0 = 1.
%      
% Discretized Brown path over [0,1] has dt = 2^(-9)
% EM uses 5 different timesteps: 16dt, 8dt, 4dt, 2dt, dt
% Examine strong convergence at T = 1: E|X_L-X(T)|

randn('state',100)
lambda = 2; mu = 1; X0 = 1;   % problem parameter
T = 1; N = 2^9; dt = T/N;
M = 1000;                     % number of paths sampled

Xerr = zeros(M,5);            % each cloume of Xerr stores the errors of paths of EM Method on different timesteps

%%%%%%--------------------------------%%%%%
%%%% 以下是自己所写的，对文献中的代码改进，向量运算代替循环可以加快速度 %%%%
dW = sqrt(dt)*randn(M,N);   % true solution and EM 用同一Brown paths
W = cumsum(dW,2);
Xtrue = X0*exp((lambda-0.5*mu^2)+mu*W(:,end)); 
for p = 1:5
    R = 2^(p-1); Dt = R*dt; L = N/R; 
    Xtemp = repmat(X0,[M,1]);  %每个循环需重新赋初值，按不同步长重新计算（不应该放在外循环之外）
    for j = 1:L
        Winner = sum(dW(:,R*(j-1)+1:R*j),2);  % 按第2个维度求和，得到的是列向量
        Xtemp = Xtemp + Dt*lambda*Xtemp + mu*Xtemp.*Winner;  
    end
    Xerr(:,p) = abs(Xtemp - Xtrue);
end  
%%%%--------------------------------%%%%

%%%% 以下被注释屏蔽的是论文中的代码，有三个嵌套循环，运行速度慢于上面的代码。
% for s = 1:M                   % the number of sample paths
%     dW = sqrt(dt)*randn(1,N);
%     W = cumsum(dW);
%     Xtrue = X0*exp((lambda-0.5*mu^2)+mu*W(end));   % the exact solution at t = 1
%     for p = 1:5
%         R = 2^(p-1); Dt = R*dt; L = N/R;      % 取不同步长计算在t = 1处EM的解（要画出Error-Dt的图像，所以需要选取不同的步长计算）
%         Xtemp = X0;
%         for j = 1:L
%             Winner = sum(dW(R*(j-1)+1:R*j));
%             Xtemp = Xtemp + Dt*lambda*Xtemp + mu*Xtemp*Winner;
%         end
%         Xerr(s,p) = abs(Xtemp-Xtrue);       % store the error at t = 1
%     end
% end
%%%% ------------------------------------%%%%%%%

Dtvector = dt*(2.^[0:4]);
subplot(221)
loglog(Dtvector,mean(Xerr),'b*-'); hold on;
loglog(Dtvector,(Dtvector.^(.5)),'r--'); hold off;  %refence slope of 1/2
axis([1e-3 1e-1 1e-4 1]);
xlabel('\Delta t'); ylabel('Sample average of |X(T)-X_L|');
title('EMstrong.m','Fontsize',10);

%%%% Least squares fit of error = C*Dt^q %%%%
A = [ones(5,1),log(Dtvector)']; rhs = log (mean(Xerr)');
sol = A\rhs; q = sol(2)
resid = norm(A*sol - rhs)
