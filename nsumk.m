function [x] = nsumk(n,k)
if (k==0)*(n==0)
   x=[0]; 
elseif (k==1)*(n>=0)    
   x=[n];
else    
m = nchoosek(k+n-1,k-1);
div = [zeros(m,1),nchoosek((1:(k+n-1))',k-1),ones(m,1)*(k+n)];
x = diff(div,1,2)-1;
end

