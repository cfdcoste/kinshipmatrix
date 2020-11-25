clear all; clc

%% Model with more than one offspring per litter: the case of ground squirrels: age III 

S=[ 0 0 0; .29 0 0; 0 .41 .42];
B=[ 0 .5628 .7144;0 0 0; 0 0 0];
Psi=[ 0 2.3 3;0 0 0; 0 0 0];
F=B.*Psi;
F1=max(0,Psi-1);
A=S+F ;
omega=7;
ts=[length(S),omega]; %trait structure (number of stage classes, number of age classes)

[~,F1] = stageage(S,F1,omega);
[S,F] = stageage(S,F,omega); 

A=S+F;

%% Reference Leslie Matrix
[w,lam]=eigs(A,1);
w=w/sum(w);w(w<0)=0;

L=full(fold(A,w,ts, 1));
F1=full(fold(F1,w,ts, 1));
7
S=zeros(size(L)); F=zeros(size(L));
S(2:end,:)=L(2:end,:)
F(1,:)=L(1,:)

A=L;
[w,lam]=eigs(A,1);w=w/sum(w);
 %% genealogical matrices

[P,~]=geneamat(A);
P % backward genealogical markov chain


PF= (F./A).*P;PF(isnan(PF))=0 
PS=P-PF

%%

Ng=(eye(size(A))-PS)^-1; % genealogical fundamental matrix
G=Ng*PF;  %previous-generation genealogical markov chain 
E=((lam*eye(size(S))-S)^-1)*F ; %EulerLotka Matrix
[GE,~]=geneamat(E);

%% computations of kinship


K21age=kin(2,1,S,F,F1,PS,PF,7,7*2) %K(2,1) converged

% ones(1,length(F))*kin(1,0,S,F,F1,PS,PF,-90,90,90*1) % expected number of (1,0) kin by stages of ego
% ones(1,length(F))*kin(2,1,S,F,F1,PS,PF,-90,90,90*2)% expected number of (2,1) kin by stages of ego

%ones(1,length(F))*kin(1,0,S,F,F1,PS,PF,7,7*1)*w % expected number of (1,0) for ego at random
ones(1,length(F))*kin(2,1,S,F,F1,PS,PF,7,7*2)*w% expected number of (2,1) for ego at random

%%
h=bar3(K21age);
shading interp
for i = 1:length(h)
     zdata = get(h(i),'Zdata');
     set(h(i),'Cdata',zdata)
     set(h,'EdgeColor','k')
end
xlabel('age of ego') 
ylabel('age of aunt') 
zlabel('expected number of aunts') 
%legend({'age 1','age 2','age 3','stable-state abundances'},'Location','southeast')
