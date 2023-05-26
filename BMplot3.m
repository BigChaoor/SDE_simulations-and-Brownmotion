%模拟Brown运动的函数的轨道，并模拟其均值的轨道
randn('state',100);  
T = 1; N = 500; dt = T/N; t = dt:dt:1;

M = 1000;   %M paths simultaneously
dW = sqrt(dt)*randn(M,N);
W = cumsum(dW,2);    %cumulative sums across the 2th dimension(即按矩阵的列累加求和)
U = exp(repmat(t,[M,1])+0.5*W);  %repmat(t,[M,1]) produce an M*N array whose ith row are all copies of t
Umean = mean(U);      %mean(U) 对矩阵的列取均值,对向量（无论行列向量）操作也是取均值
plot([0,t],[1,Umean],'b-'); hold on;    %plot mean over M paths
plot([0,t],[ones(5,1),U(1:5,:)],'r--'); hold off   %plot 5 individual paths
xlabel('t','FontSize',12);
ylabel('U(t)','FontSize',12,'Rotation',0);
legend('mean of 1000 paths','5 individual paths');
averr = norm(Umean-exp(9*t/8),'inf')    %sample error

%%%%笔记
% repmat(U,2) 将U作为元素块复制成2*2的矩阵
% repmat(U,2,3) or repmat(U,[2,3]) 将U作为块复制成2*3矩阵