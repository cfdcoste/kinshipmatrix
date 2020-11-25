function [LP] = kin(g,q,S,F,F1,PS,PF,alphamax,tmax)

LP=zeros(size(S));
%for t=max(q+1,g):((1+max(q,g))*omega) % WE NEED Q+1 BECAUSE OF function down where you cant have t=q (cause then you only have fert events and so the MRCA daughter is the ancestor ...)
for t=0:tmax
    
    
    %down(q,t,S,F,PS)*up(g,t,PS,PF)';
LP=LP+ down(q,t,S,F,F1,PS,alphamax)*up(g,t,PS,PF)';
end
LP;
end

