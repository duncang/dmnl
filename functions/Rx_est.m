%autocorrelation function

function [Rxall]=Rx_est(X,M)
N=length(X);
Rx=zeros(1,M+1);    	              
for m=1:M+1,                          
  for n=1:N-m+1,
    Rx(m)=Rx(m)+X(n)*X(n+m-1);        
  end;
  Rx(m)=Rx(m)/(N-m+1); 		      
end;
for i=1:M,     
    Rxall(i)=Rx(M+2-i);
end
Rxall(M+1:2*M+1)=Rx(1:M+1); 
