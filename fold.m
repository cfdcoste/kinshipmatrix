function [ Mfolded ] = fold(M,w,traitsize, apos) %% folds M(lam,w) on the traits to be found in positions apos of traitsize, i.e. those traits will "disappear"
a=length(apos); % those a traits will disappear
n=length(traitsize); % (n-a) traits will remain

sign=1:n; [Lia,Locb] = ismember(apos,sign);test(min(Lia),1);
sign(Locb)=[];sign=[apos sign] ; % reordering trait order so traits to disappear appear in the first a positions, remaining traits in the (n-a) last positions

order=reshape(1:prod(traitsize),traitsize); %all state-indices of M, in M multidimensional notation

I = eye(n,n);P = I(sign,:); traitsize2=(P*traitsize')';% traitsize vector for new ordering of traits 
order=permute(order,sign); signstate=order(:);I = eye(prod(traitsize),prod(traitsize));PermutStates = sparse(I(signstate,:));       % permutation matrix reordering traits in M to match traitsize2

W=reshape(w,traitsize);% multidimensional version of w
W2=permute(W,sign); % permutation of dimensions in w to match new ordering

wnewstates=sum(reshape(W2,[prod(traitsize2(1:a)) traitsize2(a+1:n)]) ,1); % w for remaining traits after folding over the other traits
wnewstates2=repmat(wnewstates(:),1,prod(traitsize2(1:a)))'; % replication of this new w (a row vector) over folded traits  
weightmat=sparse(repmat( (W2(:)./wnewstates2(:))' , prod(traitsize2),1));weightmat(isnan(weightmat))=0;weightmat=sparse(weightmat);
%matrix of all weights to be applied to M : for each state w(state)/w(all states sharing same remaining traits) 
Mweighted=sparse( (PermutStates*M*PermutStates').*weightmat); %M, reordered and weighted by the appropriate weights
Pr=sparse(kron(eye(prod(traitsize2(a+1:n))),ones(1,prod(traitsize2(1:a))))); %trait reduction "permutation" matrix
Mfolded= Pr *Mweighted* Pr'; 
end

