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
S=sparse(S);
F=sparse(F);
F1=sparse(F1);
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

%% S^omgea=0 and PS^omgea=0 
%sum(sum(S^(omega)))
%sum(sum(PS^(omega)))

%% computations of kinship

spup(1,3,PS,PF); %U(1,3)
spdown(1,2,S,F,F1,PS,-100,+100); % D(1,2)

 
spdown(2,100,S,F,F1,PS,-100,-1) % Y(2,100)
spdown(2,50,S,F,F1,PS,0,0) % YO(2,50)

% spdown(2,50,S,F,F1,PS,-90,90) % D(2,50) with alphamin=-90 and
% ...converged
% spdown(2,50,S,F,F1,PS,-900,900) % D(2,50) converged
sum(sum(spdown(2,50,S,F,F1,PS,-90,90) -spdown(2,50,S,F,F1,PS,-100,100))) %indeed

% spkin(1,0,S,F,F1,PS,PF,-90,90,90*1) %K(1,0) with alphamin=-90 and alphamax=90 and tmax=90 ...converged
% spkin(1,0,S,F,F1,PS,PF,-100,100,2*100) %K(1,0) converged
spkin(1,0,S,F,F1,PS,PF,-90,90,90*1) -spkin(1,0,S,F,F1,PS,PF,-100,100,2*100) %K(1,0) indeed


%spkin(2,0,S,F,F1,PS,PF,-90,90,90*2) %K(2,0) with alphamin=-90 and alphamax=90 and tmax=180

%spkin(2,1,S,F,F1,PS,PF,-90,90,90*2) %K(2,1) with alphamin=-90 and alphamax=90 and tmax=180

%ones(1,length(F))*spkin(1,0,S,F,F1,PS,PF,-90,90,90*1) % expected number of (1,0) kin by stages of ego
reshape((ones(1,length(F))*spkin(1,0,S,F,F1,PS,PF,-90,90,90*1))',ts);
%ones(1,length(F))*spkin(2,0,S,F,F1,PS,PF,-90,90,90*2)% expected number of (2,0) kin by stages of ego
%ones(1,length(F))*spkin(2,1,S,F,F1,PS,PF,-90,90,90*2)% expected number of (2,1) kin by stages of ego

ones(1,length(F))*spkin(1,0,S,F,F1,PS,PF,-90,90,90*1)*w % expected number of (1,0) for ego at random
ones(1,length(F))*spkin(2,0,S,F,F1,PS,PF,-90,90,90*2)*w % expected number of (2,0) for ego at random
ones(1,length(F))*spkin(2,1,S,F,F1,PS,PF,-90,90,90*2)*w% expected number of (2,1) for ego at random

%% plot number of expected alive granddaughters in 50yrs time per stage and age
figure
plot(full(reshape(sum(spdown(2,50,S,F,F1,PS,-90,90)),ts))' )% D(2,50) with alphamin=-90 and alphamax=90 summed over all kin stages... number of grd-chil per stage of ego
xlabel('age') 
ylabel('expected number of granddaughters at age 50') 
legend({'offspring','juvenile','reproductive adult','post-reproductive adult'},'Location','northwest')

%% plot probabilty of mother alive per age and stage 
vec=reshape((ones(1,length(F))*spkin(1,0,S,F,F1,PS,PF,-90,90,90*1)),ts)
vec(vec()==0)=nan
figure
plot(vec')
xlabel('age') 
ylabel('probabilty of mother alive') 
legend({'offspring','juvenile','reproductive adult','post-reproductive adult'},'Location','southwest')

%%

