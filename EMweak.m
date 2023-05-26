% EMweak Test weak convergence of Euler-Maruyama
% 
% Solves  dX = lambda*X dt + mu*X dW, X(0) = X0
%       where lambda = 2, mu = 1 and X0 = 1.
% 
% EM uses 5 different timesteps: 2^(p-10), p = 1,2,3,4,5.
% Examine weak convergence at T = 1: |E(X_L)-E(X(T))|
%
% different paths are used for each EM timestep.
randn('state',100);
lambda = 2; mu = 0.1; X0 = 1; T = 1;
M = 50000;     % number of paths sampled

Xem = zeros(5,1);
for p = 1:5
    Dt = 2^(p-10); L = T/Dt;
    Xtemp =X0*ones(M,1);   % or Xtemp = repmat(X0,[M,1]);
    for j = 1:L
%         Winner = sqrt(Dt)*sign(randn(M,1));  %%use for weak E-M %%
        % WEM has weak order of convergence 1, but since it uses no
        % pathwise information, offers no strong convergence. The
        % motivation behind WEM is that random number generators that
        % sample from "two point distribution" can be made effcient than those that sample from N(0,1)
        % however, because we are using the built-in normal random number
        % generator , there is no effciency gain in this case
        Winner = sqrt(Dt)*randn(M,1);      
        % used different paths for each stepsize \Detla t
        % this is reasonable ,because we concerns only the mean of the
        % solution, and the mean of true solution is independent of W(t)
        % so the Brown paths is not required to be consistent to the ... 
        % it's the difference from EMstrong.m
        Xtemp = Xtemp + Dt*lambda*Xtemp + mu*Xtemp.*Winner;
    end
    Xem(p) = mean(Xtemp);
end
Xerr = abs(Xem - exp(lambda));

DtVector =  2.^([1:5]-10);
subplot(222)
loglog(DtVector, Xerr, 'b*-'); hold on;
loglog(DtVector, DtVector, 'r--'); hold off;   % refence slope of 1
axis([1e-3 1e-1 1e-4 1])
xlabel('\Delta t')
ylabel('E(X(T)) - Sample average of X_L')
title('EMweak.m','Fontsize',10)

%%%% Least squares fit of error = C * dt^q %%%%
A = [ones(p,1), log(DtVector)']; rhs = log(Xerr);
sol = A\rhs; q = sol(2);
resid = norm(A*sol - rhs)