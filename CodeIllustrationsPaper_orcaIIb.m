clear all; clc

%% Stage-based model via stage x age model: the case of Orcas II
% lefkovitch matrix
S=[ 0 0 0 0; 0.9775 0.9111 0 0;0 0.0736 0.9534 0; 0 0 0.0452 0.9804]; 
F=[0 0.0043 0.1132 0  ;0 0 0 0; 0 0 0 0; 0 0 0 0]; 
F1=zeros(size(F));

%% stage x age matrix
omega=90;
ts=[length(S),omega]; %trait structure (number of stage classes, number of age classes)
[~,F1] = stageage(S,F1,omega);
[S,F] = stageage(S,F,omega); 
%%
A=S+F;

%% assmptotic lambda and abundances
[w,lam]=eigs(A,1);
w=w/sum(w); w(w<0)=0;

%% genealogical matrices
[P,~]=geneamat(A);
P; % backward genealogical markov chain
P(isinf(P))=0; 

PF= (F./A).*P;PF(isnan(PF))=0; 
PS=P-PF;

%%

Ng=(sparse(eye(size(A)))-PS)^-1; % genealogical fundamental matrix
G=Ng*PF;  %previous-generation genealogical markov chain 
E=((lam*sparse(eye(size(S)))-S)^-1)*F; %EulerLotka Matrix
[GE,~]=geneamat(E);


%% computations of kinship


ones(1,length(F))*kin(1,0,S,F,F1,PS,PF,90,90*1)*w % expected number of (1,0) for ego at random
ones(1,length(F))*kin(2,0,S,F,F1,PS,PF,90,90*2)*w % expected number of (1,0) for ego at random
ones(1,length(F))*kin(2,1,S,F,F1,PS,PF,90,90*2)*w % expected number of (1,0) for ego at random

%% plot number of expected alive granddaughters in 50yrs time per stage and age
figure
plot(full(reshape(sum(spdown(2,50,S,F,F1,PS,-90,90)),ts))' )% D(2,50) with alphamin=-90 and alphamax=90 summed over all kin stages... number of grd-chil per stage of ego
xlabel('age') 
ylabel('expected number of granddaughters at age 50') 
legend({'offspring','juvenile','reproductive adult','post-reproductive adult'},'Location','northwest')

%% plot probabilty of mother alive per age and stage 
vec=reshape((ones(1,length(F))*spkin(1,0,S,F,F1,PS,PF,-90,90,90*1)),ts);
vec(vec()==0)=nan;
figure
plot(vec')
xlabel('age') 
ylabel('probabilty of mother alive') 
legend({'offspring','juvenile','reproductive adult','post-reproductive adult'},'Location','southwest')

%%

