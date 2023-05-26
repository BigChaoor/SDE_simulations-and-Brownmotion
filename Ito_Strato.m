% Approximate stochastic integrals模拟随机积分,并求T时刻模拟值的误差
% Ito and stratonovich integrals of W

randn('state',100);
T = 1; N =500; dt = T/N;

dW = sqrt(dt)*randn(1,N);
W = cumsum(dW); 

ito = sum([0,W(1:end-1)].*dW);
strato = sum((0.5*([0,W(1:end-1)]+W) + 0.5*sqrt(dt)*randn(1,N)).*dW);

itoerr = abs(ito - 0.5*(W(end)^2-T))
stratoerr = abs(strato - 0.5*W(end)^2)