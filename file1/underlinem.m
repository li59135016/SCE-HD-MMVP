function z=underlinem(u)
%��Best62.pdf��17ҳ��step3 ��\underline{m}(u)����
global eigenY;global y;

z=-(1-y)/u+y*sum(eigenY-u)/length(eigenY);

