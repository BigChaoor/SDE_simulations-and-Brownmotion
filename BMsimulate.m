% simulate Brown Motion
% Ԥ����ģ���������n ,����������ԣ���ͼ�ν���������ʾ�ĵ����С��n
% �����˶���ʱ����s = 0.04
n = 100; s = 0.04;
sqrts = sqrt(s);
x = 0; y = 0;     % �������ڵĳ�ʼλ��
h = plot(x,y,'.'); 
pause(0.5)   % ��ͣ0.5��
axis([-1 1 -1 1])
axis square
grid off
set(h,'EraseMode','xor','Markersize',12);
T = 0;
while T < 10
    x = x + sqrts*randn(n,1);
    y = y + sqrts*randn(n,1);
    set(h,'XData',x,'YData',y)   % ���������ݲ���ͼ
    T =T+s;
    % ���������drawnow����set()֮ǰ�������ʾ��ͼ���ǿհ׵ģ���Ϊִ���˸���
    drawnow   % ����ͼ�δ��ڼ����Ӽ����˴���pause(DeltaT)�滻������ʾ���Ƶ�Ч��
end