% simulate Brown Motion
% 预设所模拟的粒子数n ,但由于随机性，在图形界面中所显示的点可能小于n
% 粒子运动的时间间隔s = 0.04
n = 100; s = 0.04;
sqrts = sqrt(s);
x = 0; y = 0;     % 粒子所在的初始位置
h = plot(x,y,'.'); 
pause(0.5)   % 暂停0.5秒
axis([-1 1 -1 1])
axis square
grid off
set(h,'EraseMode','xor','Markersize',12);
T = 0;
while T < 10
    x = x + sqrts*randn(n,1);
    y = y + sqrts*randn(n,1);
    set(h,'XData',x,'YData',y)   % 设置新数据并画图
    T =T+s;
    % 若将此语句drawnow置于set()之前，最后显示的图像将是空白的，因为执行了更新
    drawnow   % 更新图形窗口及其子级，此处用pause(DeltaT)替换可以显示类似的效果
end