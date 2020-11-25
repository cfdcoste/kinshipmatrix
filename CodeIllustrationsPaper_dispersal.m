clear all; clc

%% Model with several classes of offspring: a theoretical model with dispersion

S=[ 0 0 0 0; .4*(1-.5) .5 .3*.4 0;0 0 0 0;  .4*.5 0 .3*(1-.4) 0.4];

B=[ 0 .8 0 0; 0 0 0 0; 0 0 0 .7; 0 0 0 0];
Psi=[ 0 2.3 0 0; 0 0 0 0; 0 0 0 2.3; 0 0 0 0];
F=B.*Psi;

F1=max(0,Psi-1);

%sum(S) %survival probabilities
%S./repmat(sum(S),4,1);% dispersal rates

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

Ng=(eye(size(A))-PS)^-1; % genealogical fundamental matrix
G=Ng*PF  %previous-generation genealogical markov chain 
G^10
E=((lam*eye(size(S))-S)^-1)*F  %EulerLotka Matrix
[GE,~]=geneamat(E)

%% computations of kinship

kin(1,0,S,F,F1,PS,PF,10,10*1); 

kin(2,1,S,F,F1,PS,PF,10,10*2); 

ones(1,length(F))*kin(1,0,S,F,F1,PS,PF,10,10*1)*w % expected number of (1,0) for ego at random
ones(1,length(F))*kin(2,1,S,F,F1,PS,PF,10,7*10)*w% expected number of (2,1) for ego at random


%%  patch A matrix
A=A(1:2,1:2);
S=S(1:2,1:2);
F=F(1:2,1:2);
F1=F1(1:2,1:2);

[w,lam]=eigs(A,1);
w=w/sum(w);

[P,Q]=geneamat(A);
PF= (F./A).*P;PF(isnan(PF))=0;
PS=P-PF;
kin(1,0,S,F,F1,PS,PF,10,10*1); 

kin(2,1,S,F,F1,PS,PF,10,10*2); 

ones(1,length(F))*kin(1,0,S,F,F1,PS,PF,10,10*1)*w % expected number of (1,0) for ego at random
ones(1,length(F))*kin(2,1,S,F,F1,PS,PF,10,10*2)*w% expected number of (2,1) for ego at random
