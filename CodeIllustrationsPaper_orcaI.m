clear all; clc

%% Stage-based model: the case of Orcas I
% lefkovitch matrix
S=[ 0 0 0 0; 0.9775 0.9111 0 0;0 0.0736 0.9534 0; 0 0 0.0452 0.9804]; 
F=[0 0.0043 0.1132 0  ;0 0 0 0; 0 0 0 0; 0 0 0 0]; 
F1=zeros(size(F));

%%

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


down(2,50,S,F,F1,PS,90) % 
down(2,50,S,F,F1,PS,900) 

sum(sum(abs(down(2,50,S,F,F1,PS,90) -down(2,50,S,F,F1,PS,900))))


%%
kin(1,0,S,F,F1,PS,PF,90,90*1) %K(1,0) with alphamin=-90 and alphamax=90 and tmax=90
kin(1,0,S,F,F1,PS,PF,900,900*1) %K(1,0) converged

sum(sum(abs(kin(1,0,S,F,F1,PS,PF,90,90*1) -kin(1,0,S,F,F1,PS,PF,900,900*1))))

%%

kin(2,0,S,F,F1,PS,PF,90,90*2) %K(1,0) with alphamin=-90 and alphamax=90 and tmax=90
kin(2,0,S,F,F1,PS,PF,900,900*2) %K(1,0) converged

sum(sum(abs(kin(2,0,S,F,F1,PS,PF,90,90*2) -kin(2,0,S,F,F1,PS,PF,900,900*2))))


%%
kin(2,1,S,F,F1,PS,PF,90,90*2) %K(1,0) with alphamin=-90 and alphamax=90 and tmax=90
kin(2,1,S,F,F1,PS,PF,900,900*2) %K(1,0) converged

sum(sum(abs(kin(2,1,S,F,F1,PS,PF,90,90*2) -kin(2,1,S,F,F1,PS,PF,900,900*2))))

%%




ones(1,length(F))*kin(1,0,S,F,F1,PS,PF,90,90*1) % expected number of (1,0) kin by stages of ego
ones(1,length(F))*kin(2,0,S,F,F1,PS,PF,90,90*2)% expected number of (2,0) kin by stages of ego
ones(1,length(F))*kin(2,1,S,F,F1,PS,PF,90,90*2)% expected number of (2,1) kin by stages of ego


%%
ones(1,length(F))*kin(1,0,S,F,F1,PS,PF,90,90*1)*w % expected number of (1,0) for ego at random
ones(1,length(F))*kin(2,0,S,F,F1,PS,PF,90,90*2)*w % expected number of (2,0) for ego at random
ones(1,length(F))*kin(2,1,S,F,F1,PS,PF,90,90*2)*w% expected number of (2,1) for ego at random


%%  North Atlantic Right Whale instead

S=[ 0 0 0 0; 0.92 0.86 0 0;0 0.08 0.80 0.83; 0 .02 0.19 0.0]; % for lefkovitvh
F=[0 0 0.0 0.335  ;0 0 0 0; 0 0 0 0; 0 0 0 0]; % for lefkovitvh
F1=zeros(size(F));

%%
A=S+F 

%% assmptotic lambda and abundances

[w,lam]=eigs(A,1);
w=w/sum(w);

 %% genealogical matrices

[P,~]=geneamat(A);
PF= (F./A).*P;PF(isnan(PF))=0 ;
PS=P-PF;

kin(2,1,S,F,F1,PS,PF,90,90*2) %K(1,0) with alphamin=-90 and alphamax=90 and tmax=90
kin(2,1,S,F,F1,PS,PF,900,900*2) %K(1,0) converged
sum(sum(abs(kin(2,1,S,F,F1,PS,PF,90,90*2) -kin(2,1,S,F,F1,PS,PF,900,900*2))))

ones(1,length(F))*kin(2,1,S,F,F1,PS,PF,90,90*2)% expected number of (2,1) kin by stages of ego

ones(1,length(F))*kin(2,1,S,F,F1,PS,PF,90,90*2)*w% expected number of (2,1) for ego at random


