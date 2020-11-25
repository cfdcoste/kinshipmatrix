function [ P,Q] = geneamat( M)

[wM,lamM]=eigs(M,1);
wM=wM/sum(wM);
wM(wM<0)=0;
wM*(wM.^-1)';
P=(1/lamM)*M.*((wM.^-1)*wM');
P(isnan(P))=0; 

[vM,lamM]=eigs(M',1);
vM=vM/(wM'*vM);
%Q=(1/lamM)*M'.*((vM.^-1)*vM');
Q=(1/lamM)*M.*(vM*(vM'.^-1));

end

