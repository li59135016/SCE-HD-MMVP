function r=qutile(F)
global alpha;
c=2*(1-alpha)^2;a=2*alpha-1;
r=a+sqrt(c./(2.*((c/(2*(alpha-a)^2))-F)));


