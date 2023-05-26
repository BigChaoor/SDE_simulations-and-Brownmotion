% STABLE Mean-Square and asymptotic stability test for EM
% 
% SDE is the linear test problem dX = lambda*X + mu*X dW, X(0) = X0
%     where lambda and mu is constants, X0 = 1

randn('state',100);
T = 20; M = 50000; X0 = 1;
ltype = {'b-', 'r--', 'm-.'};     % linetypes for plot 
% Ҳ�ɻ�����������%% ltype = ['b-', 'r--', 'm-.'] ,�������ʱ�����ltype(k)����
% ѧϰ��̹淶�����ͺ���ɫ��ѡ��

subplot(211)    %%%%% mean-square  %%%%%%
lambda = -3; mu = sqrt(3);     % problem parameters
for k = 1:3
    Dt = 2^(1-k);
    N = T/Dt;
    Xms = zeros(1,N); Xtemp = X0*ones(M,1);
    for j = 1:N
        Winner = sqrt(Dt)*randn(M,1);
        Xtemp = Xtemp + Dt*lambda*Xtemp + mu*Xtemp.*Winner;
        Xms(j) = mean(Xtemp.^2);     % mean-square estimate
    end
    semilogy(0:Dt:T,[X0,Xms],ltype{k},'linewidth',2), hold on   % yȡ������,���Ƶ���semilogx()
end
legend('\Delta t = 1','\Delta t = 1/2','\Delta t = 1/4')  % ѧϰ��ѧ����������룬\Delta t��\surd 3��
title('Mean-square: \lambda = -3, \mu = \surd 3','Fontsize',12)
ylabel('E[X^2]','Fontsize',12), axis([0,T,1e-20,1e+20]),hold off

subplot(212)  %%%% Asymptotic: a single path %%%%
T = 500;
lambda = 0.5; mu = sqrt(6);
for k = 1:3
    Dt = 2^(1-k);
    N = T/Dt;
    Xasy = zeros(1,N); Xtemp = X0;
    for j = 1:N
        Winner = sqrt(Dt)*randn;
        Xtemp = Xtemp + Dt*lambda*Xtemp + mu*Xtemp*Winner;
        Xasy(j) = abs(Xtemp);
    end
    semilogy(0:Dt:T,[X0,Xasy],ltype{k},'Linewidth',2),hold on
end
legend('\Delta t = 1','\Delta t = 1/2','\Delta t = 1/4')
title('Single Path: \lambda = 1/2, \mu = \surd 6','Fontsize',12)
ylabel('|X(t)|','Fontsize',12)
axis([0,T,1e-50,1e+100]),hold off
