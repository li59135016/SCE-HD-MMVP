function z=underlinem(u)
%见Best62.pdf中17页，step3 中\underline{m}(u)函数
global eigenY;global y;

z=-(1-y)/u+y*sum(eigenY-u)/length(eigenY);

