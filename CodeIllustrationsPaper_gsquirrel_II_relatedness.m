clear all; clc

%% Model with more than one offspring per litter: the case of ground squirrels II relatedness  

S=[ 0 0 0; .29 0 0; 0 .41 .42]
B=[ 0 .5628 .7144;0 0 0; 0 0 0]
Psi=[ 0 2.3 3;0 0 0; 0 0 0]
F=B.*Psi;
F1=max(0,Psi-1);
A=S+F 
[w,lam]=eigs(A,1);
w=w/sum(w);
[P,~]=geneamat(A);
P % backward genealogical markov chain
PF= (F./A).*P;PF(isnan(PF))=0; 
PS=P-PF;

%% relatedness computations

%load('Kstock.mat','Kstock')
dmax=6; % 

%return
Kstock=zeros([size(F) dmax+1 dmax+1]);


for d=0:dmax %0
d
for g=0:d
g
q=d-g;
Kstock(:,:,g+1,q+1)=kin(g,q,S,F,F1,PS,PF,7,7*max(g,q));
%save('Kstock.mat','Kstock');
end
end

%%
Kd=zeros([size(F) dmax+1]);

for d=0:dmax
for j=0:d
    Kd(:,:,d+1)=Kd(:,:,d+1)+Kstock(:,:,j+1,d-j+1);
end
end

Kd(:,:,7) % is K(6) 


%%  individual toal coefficient of relationship

TR=zeros(size(Kd));
cum=TR(:,:,1);
for d=1:dmax
cum= (2^(-d))*Kd(:,:,d+1) + cum ; 
TR(:,:,d+1)=cum;
end

TR(:,:,7)% is TR(6)

%%
tr=[];meantr=[];
for d=1:dmax
tr(d,:)=ones(1,length(S))*TR(:,:,d);
meantr(d)=ones(1,length(S))*TR(:,:,d)*w;
end
tr(:,end+1)=meantr';


tr=[];meantr=[];
for d=1:dmax
tr(d+1,:)=ones(1,length(S))*TR(:,:,d+1);
meantr(d+1)=ones(1,length(S))*TR(:,:,d+1)*w;
end
tr(:,end+1)=meantr';

%% figure
bar(tr(2:end,:))
xlabel('dmax') 
ylabel('cumulative (over d up to dmax) expected number of kin') 
legend({'age 1','age 2','age 3','stable-state abundances'},'Location','southeast')


