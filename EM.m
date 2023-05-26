% EM Euler-Maruyama method on Black-Scholes SDE
%
% Black-Scholes model dX = lamda*X dt+ mu*X dW, X(0) = X0
% where lamda = 2, mu = 1 and X0 = 1
% Discretized Brownian path over [0,1] has dt = 2^(-8)
% Euler-Maruyama uses timestep R*dt.
randn('state',100)
lamda = 2; mu = 1; X0 = 1;   % problem parameters
T = 1; N = 2^8; dt = 1/N;
dW  = sqrt(dt)*randn(1,N);   % Brownian increments
W = cumsum(dW);   % discreatized Brown path

Xtrue = X0*exp((lamda-0.5*mu^2)*(dt:dt:T)+mu*W);
plot(0:dt:T,[X0,Xtrue],'m-'); hold on;

R = 4; Dt = R*dt; L = N/R;   % L is the EM steps of timestep Dt = R*dt
Xem = zeros(1,L);
Xtemp = X0;
% % 另一写法,此处不太适合，因为Xem存储的不包括处置，以下写法适用于将初值存入Xem向量的第一个元素
% Xem(1) = X0;
% for j = 1:L-1
%     Winner = sum(dW(R*(j-1)+1:R*(j)));
%     Xem(j+1) = Xem(j) + Dt*lamda*Xem(j) + mu*Xem(j)*Winner;
% end
for j = 1:L
    Winner = sum(dW(R*(j-1)+1:R*j));
    Xtemp = Xtemp + Dt*lamda*Xtemp + mu*Xtemp*Winner;
    Xem(j) = Xtemp;
end

plot(0:Dt:T,[X0,Xem],'r--*'); hold off;
xlabel('t','FontSize',12)
ylabel('X','FontSize',12,'Rotation',0,'HorizontalAlignment','right')
 
emerr = abs(Xem(end)-Xtrue(end))