function [PP] = spup(g,t,PS,PF)
PP=sparse(zeros(size(PS)));
if (t<g)+(t>0)*(g==0)   
    PP;
elseif (g==0)*(t==0) 
    PP=sparse(eye(size(PP)));
else 
part = nsumk(t-g,g); %partitions of survival events around the g fertility events
for i=1:size(part,1)
PP1=sparse(eye(length(PS)));
for j=1:g
PP1=PP1*(PS^part(i,j))*PF;
end
PP=PP+PP1;
end
PP;
end
end