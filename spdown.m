function [PP] = spdown(q,t,S,F,F1,PS,alphamin,alphamax)
if q==0
    PP=S^max(t,0)*((PS')^(-min(t,0)));
else
Ptemp=sparse(zeros(size(S)));
for alpha=alphamin:alphamax
try part=nsumk(t-q+alpha,q);
for i=1:size(part,1)
LPtemp=sparse(eye(size(S)));
    for j=1:q-1
     LPtemp=(S^part(i,j))*F*LPtemp;
     end
Ptemp=Ptemp+LPtemp*(S^part(i,q))*(F1^(alpha==0))*(F^(alpha~=0))*((PS')^max(alpha,0))*(S^(-min(alpha,0)));
end
catch    
end

PP=Ptemp;

end
end