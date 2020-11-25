function [PP] = up(g,t,PS,PF)
PP=zeros(size(PS));
if (t<g)+(t>0)*(g==0)   
    PP;
elseif (g==0)*(t==0) 
    PP=eye(size(PP));
else 
part = nsumk(t-g,g); %partitions of survival events around the g fertility events
for i=1:size(part,1)
PP1=eye(length(PS));
for j=1:g
PP1=PP1*(PS^part(i,j))*PF;
end
PP=PP+PP1;
end
PP;
end
end