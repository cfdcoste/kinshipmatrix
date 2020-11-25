function [Ssa,Fsa] = stageage(S,F,omega)
s=length(S); %number of stages
ts=[s omega]; % trait structure of stage-age matrix
q=prod(ts); % dimension of stage-age matrix
Ssa=zeros(q,q); Fsa=zeros(q,q);
for i=1:omega-1
Ssa(i*s+1:i*s+s,(i-1)*s+1:(i-1)*s+s)=S;
end
for i=1:omega
Fsa(1:s,(i-1)*s+1:(i-1)*s+s)=F;
end
end