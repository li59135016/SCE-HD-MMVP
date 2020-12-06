function r=Fh(x)%h(t)的积分上限函数
global alpha;
c=2*(1-alpha)^2;a=2*alpha-1;
r=c/(2*(alpha-a)^2)-c./(2*(x-a).^2);



