function [x] = nsumrecuk(n,k)
if (k==0)*(n==0)
   x=[0]; 
elseif (k==1)*(n>=0)     
   x=[n];
else
    temp=[];
    for a=0:n
    temp1=nsumrecuk(a,k-1);   
    temp1=[temp1 ones(size(nsumrecuk(a,k-1),1),1)*(n-a)];
    temp=[temp;temp1];
    end
x=sortrows(temp);
end