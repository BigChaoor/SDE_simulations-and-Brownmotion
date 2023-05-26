%模拟Brown运动的轨道
randn('state',100)   %set the state of randn，随机数种子可控制每次运行产生的随机数一样
T = 1; N = 500; dt = T/N;

dW = sqrt(dt)*randn(1,N);  %increments
W = cumsum(dW);      %cumulative sum

plot(0:dt:T,[0,W],'r-');  
xlabel('t','FontSize',16);
ylabel('W(t)','FontSize',16,'Rotation',0);