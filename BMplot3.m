%ģ��Brown�˶��ĺ����Ĺ������ģ�����ֵ�Ĺ��
randn('state',100);  
T = 1; N = 500; dt = T/N; t = dt:dt:1;

M = 1000;   %M paths simultaneously
dW = sqrt(dt)*randn(M,N);
W = cumsum(dW,2);    %cumulative sums across the 2th dimension(������������ۼ����)
U = exp(repmat(t,[M,1])+0.5*W);  %repmat(t,[M,1]) produce an M*N array whose ith row are all copies of t
Umean = mean(U);      %mean(U) �Ծ������ȡ��ֵ,��������������������������Ҳ��ȡ��ֵ
plot([0,t],[1,Umean],'b-'); hold on;    %plot mean over M paths
plot([0,t],[ones(5,1),U(1:5,:)],'r--'); hold off   %plot 5 individual paths
xlabel('t','FontSize',12);
ylabel('U(t)','FontSize',12,'Rotation',0);
legend('mean of 1000 paths','5 individual paths');
averr = norm(Umean-exp(9*t/8),'inf')    %sample error

%%%%�ʼ�
% repmat(U,2) ��U��ΪԪ�ؿ鸴�Ƴ�2*2�ľ���
% repmat(U,2,3) or repmat(U,[2,3]) ��U��Ϊ�鸴�Ƴ�2*3����