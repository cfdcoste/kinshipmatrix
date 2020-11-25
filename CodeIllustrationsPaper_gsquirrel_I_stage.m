clear all; clc

%% Model with more than one offspring per litter: the case of ground squirrels I 

S=[ 0 0 0; .29 0 0; 0 .41 .42]

B=[ 0 .5628 .7144;0 0 0; 0 0 0]
Psi=[ 0 2.3 3;0 0 0; 0 0 0]
F=B.*Psi

F1=max(0,Psi-1)

A=S+F 

%% assmptotic lambda and abundances

[w,lam]=eigs(A,1);
w=w/sum(w);

 %% genealogical matrices

[P,~]=geneamat(A);
P % backward genealogical markov chain


PF= (F./A).*P;PF(isnan(PF))=0 
PS=P-PF


%%
up(1,3,PS,PF) %U(1,3)

%%

Ng=(eye(size(A))-PS)^-1 % genealogical fundamental matrix
G=Ng*PF  %previous-generation genealogical markov chain 
G^10
E=((lam*eye(size(S))-S)^-1)*F ; %EulerLotka Matrix
[GE,~]=geneamat(E)

%% computations of kinship

down(2,6,S,F,F1,PS,7) % D(2,50) converged

down(2,6,S,F,F1,PS,100) % D(2,50) converged

sum(sum(abs(down(2,6,S,F,F1,PS,100) -down(2,6,S,F,F1,PS,7))))


%%

kin(1,0,S,F,F1,PS,PF,7,7*1) %K(1,0) with alphamin=-7 and alphamax=7 and tmax=7

sum(sum(abs(kin(1,0,S,F,F1,PS,PF,7,7*1) -kin(1,0,S,F,F1,PS,PF,100,100*1))))

%%


kin(2,1,S,F,F1,PS,PF,7,7*2) %K(2,1) with alphamin=-7 and alphamax=7 and tmax=7x2

sum(sum(abs(kin(2,1,S,F,F1,PS,PF,7,7*2) -kin(2,1,S,F,F1,PS,PF,100,100*2) )))


%%
ones(1,length(F))*kin(1,0,S,F,F1,PS,PF,7,7*1) % expected number of (1,0) for ego at random
ones(1,length(F))*kin(2,1,S,F,F1,PS,PF,7,7*2)% expected number of (2,1) for ego at random


%%
ones(1,length(F))*kin(1,0,S,F,F1,PS,PF,7,7*1)*w % expected number of (1,0) for ego at random
ones(1,length(F))*kin(2,1,S,F,F1,PS,PF,7,7*2)*w% expected number of (2,1) for ego at random

