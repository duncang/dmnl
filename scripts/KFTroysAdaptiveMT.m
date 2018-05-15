%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   This script file is the KF data generation and solution 
%   for Question B. 
%
%   Author: Jason Ford, Oct. 2007.
%
%   Released to students.
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%



% Some fixed parameters
T=10000;                     % 1000 seconds
dt=1;                       % Sample period
A=[1,0,dt,0;0,1,0,dt;0,0,1,0;0,0,0,1];
Q=diag([0,0,0.01,0.01]);



R=diag([10,10]);


%Q=diag([0,0,0,0]);

%R=diag([1,1]);

%Create some data storage vectors (state + measurment)
xv=zeros(4,T);                          %Note, this is a 4xT vector used to store the trajectory
yv=zeros(2,T);                          % and a 2xT vector to store measurements.

% Initial state value
x0=[10000,10000,1,0]';

% Uncomment the below line, if you want to use.
 randn('state',0)            %This forces matlab to repeat the random number sequence.


%% MAIN LOOP to generate the state and measurement processes.

%  Initial state and measurement value (k=1 values)
%           Need to do this before entering the main loop.
xv(:,1)=x0;                                                   % the initial state value;    
yv(1,1)=x0(1)   +sqrt(R(1,1))*randn(1,1);    % and initial measurement value;
yv(2,1)=x0(2)   +sqrt(R(2,2))*randn(1,1);

xold=xv(:,1);           % need to initilise before entering the main loop

for k=2:T                       %note: we start loop at k=2.
    %generate the noises.
    vk  = sqrt(Q)*randn(4,1);
    wk1 = sqrt(R(1,1))*randn(1,1);
    wk2 = sqrt(R(2,2))*randn(1,1);

    %generate new state and measurements.
    xnew=A*xold                         +vk;                    %  x_k+1=A*x_k    +   v_k
%     ynew(1,1)=sqrt(xnew(1)^2+xnew(2)^2) +wk1;                   %  y_k^1= range   + w_k^1
%     ynew(2,1)=atan(xnew(2)/xnew(1))     +wk2;                   %  y_k^2= bearing + w_k^2
     ynew(1,1)=xnew(1) +wk1;                   %  y_k^1= range   + w_k^1
     ynew(2,1)=xnew(2) +wk2;                   %  y_k^2= bearing + w_k^2


    %store and set for next loop;
    xv(:,k)=xnew;
    yv(:,k)=ynew;
    xold=xnew;
end










% Alternative compact version of the code. Is equivalent to the longer "for loop" code above.
% Does data storage at the same time as calculation.  A bit harder to read.
% for k=2:T
%        xv(:,k)=A*xv(:,k-1)                 +sqrt(Q)*randn(4,1);        %  x_k+1=A*x_k    +   v_k
%        yv(1,k)=sqrt(xv(1,k)^2+xv(2,k)^2)   +sqrt(R(1,1))*randn(1,1);   %  y_k^1= range   + w_k^1
%        yv(2,k)=atan(xv(2,k)/xv(1,k))       +sqrt(R(2,2))*randn(1,1);   %  y_k^2= bearing + w_k^2
% end

%% PLOTTING

% Generated some plots.

kv=1:T;       % a time based for plots.

%some plot commands (to check that measurements match what they should).
figure(1)
plot(kv,yv(1,:),kv,xv(1,:));
xlabel('Time (s)')
ylabel('x_k^1 (m)')
legend('Measurement','Truth')

figure(2)
plot(kv,yv(2,:),kv,xv(2,:));
xlabel('Time (s)')
ylabel('x_k^2 (m)')
legend('Measurement','Truth')

pause;
close all;

%% KF






%apriori  data

%assume initial estimates  (ie. this is the true position, plus an error
xh0=x0+[100,100,0.1,-0.1]';
p0=1000*eye(4,4);



% some storage
xhv=zeros(4,T);
xhv(:,1)=xh0;

% a few temp values used.
Pko(:,:,1)=p0;


    H=[1,0,0,0; ...
        0,1,0,0 ];
    
       
        
        
%my wrong guess of Q and R (assume I don't know them)..

Q=diag([0,0,0,0]);
R = diag([100,100 ]);
        
  
        
        
     % Qk_hat(:,:,1) = zeros(4,4);
    
    startadaptive = 200;
            
    Qk_hat(:,:,startadaptive) = Q;    
     
    Rk_hat(:,:,startadaptive) = R;
    
     rk_hat(:,startadaptive) = zeros(1,2); 
    
     qk_hat(:, startadaptive) = zeros(1,4);
  
     
     PHIk = A;          

for k=2:T
    
    
    
   
   
    
    %else
        
     % xt=A*xhv(:,k-1) + qk(k-1);      
      
     
     
     if k > startadaptive
         
         
          %state propagation       
    %if k<=100
    xt(:,k)=A*xhv(:,k-1);% +qk(:,k-1);    % the one-step ahead prediction.
    Pkminus(:,:,k) =A*(Pko(:,:,k-1)) *A' + Qk_hat(:,:,k-1);   
     else
         
         
          %state propagation       
    %if k<=100
    xt(:,k)=A*xhv(:,k-1) ;   % the one-step ahead prediction.
         
          Pkminus(:,:,k) =A*(Pko(:,:,k-1)) *A' + Q;     
          
     end
 
    
    
    %compute observation noise
    
     Nr = 10;
         
     lr = Nr;
     
    %innovation     
    rk(:,k) = yv(:,k)-H*xt(:,k);
    gammak(:,:,k) = H*Pkminus(:,:,k)*H'; 
    
    
    
             
     
    %compute observation noise
   if k > startadaptive     
         
   rk_hat(:,k) = rk_hat(:,k-1) + (1/lr)*(rk(:,k) - rk(:,k-lr));        
        
   %could use k+1 here because R computed before Q in the MT algorithm
   Rk_hat(:,:,k) = Rk_hat(:,:,k-1) + (1/(lr-1))*(       (rk(:,k)- rk_hat(:,k))*(rk(:,k)- rk_hat(:,k))' - (rk(:,k-lr)- rk_hat(:,k))*(rk(:,k-lr)- rk_hat(:,k))' + (1/lr)*(rk(:,k)- rk(:,k-lr))*(rk(:,k)- rk_hat(:,k-lr))'+ ((lr-1)/lr)*(gammak(:,:,k-lr) - gammak(:,:,k))    );
         
   
   
   
    %Rk_hat(:,:,k) = R;
         
    
            %set the diagonal elements to absolute value of their estimates
            Rk_hatdiag  = abs(diag(Rk_hat(:,:,k)));            
            Rk_hat(1,1,k) = Rk_hatdiag(1);            
            Rk_hat(2,2,k) = Rk_hatdiag(2);            
            
   end
   
   
   %compute kalman gain
    if k > startadaptive
   
     Kk=Pkminus(:,:,k)*H'*inv(gammak(:,:,k) +Rk_hat(:,:,k)); 
   
    else
        
          Kk=Pkminus(:,:,k)*H'*inv(gammak(:,:,k) +R); 
    end
     %state estimation 
     
      % xhv(:,k)=xt(:,k)+ Kk*(yv(:,k)-H*xt);  % measurement update.
    
     %  xhv(:,k)=xt(:,k)+ Kk*( rk(:,k)- rk_hat(:,k));  % measurement update.
    
     xhv(:,k)=xt(:,k)+ Kk*( rk(:,k));  % measurement update.
    
   
    Pko(:,:,k) = Pkminus(:,:,k) - Kk*H*Pkminus(:,:,k);
    
    
    
    %compute state noise
    
    N = 10;
    lq = N;
    
    
     qk(:,k) =  xhv(:,k)  - PHIk* xhv(:,k-1);
            
     delk(:,:,k) = PHIk*Pko(:,:,k-1)*PHIk' - Pko(:,:,k);
     
     
     
  if k > startadaptive  
    
        qk_hat(:,k) = qk_hat(:,k-1) + (1/lq)*(qk(:,k) - qk(:,k-lq));   %this is the mean over the data
            
            
            Qk_hat(:,:,k) = Qk_hat(:,:,k-1) + (1/(lq-1))*(     (qk(:,k) - qk_hat(:,k))*(qk(:,k) - qk_hat(:,k))' - (qk(:,k-lq) - qk_hat(:,k))*(qk(:,k-lq) - qk_hat(:,k))' + (1/lq)*(qk(:,k) - qk(:,k-lq))*(qk(:,k) - qk(:,k-lq))' + ((lq-1)/lq)*(delk(:,:,k-lq) - delk(:,:,k)) );
            
            
           % Qk_hat(:,:,k) = Q;
            
            
            %set the diagonal elements to absolute value of their estimates 
            Qk_hatdiag  = abs(diag(Qk_hat(:,:,k)));            
            Qk_hat(1,1,k) = Qk_hatdiag(1);            
            Qk_hat(2,2,k) = Qk_hatdiag(2);            
            Qk_hat(3,3,k) = Qk_hatdiag(3);            
            Qk_hat(4,4,k) = Qk_hatdiag(4);
            
      
            
            
            
  end
       % end
        
        
        
        
        
        %test
        
        
        if  k > startadaptive  
        
         for i = k-N+1:k
       % Jtemp(i) = rk(:,i)'*inv(Hsave(:,:,i)*Pksave(:,:,i)*Hsave(:,:,i)' + Rsave(:,:,i))*rk(:,i); 
        
        
         Jtemp(i) = rk(:,i)'*inv(H*Pkminus(:,:,i)*H' + Rk_hat(:,:,i))*rk(:,i); 
         
       %  Jtemp(i) = rk(:,i)'*inv(H*Pkminus(:,:,i)*H' + R)*rk(:,i); 
        
        
         end                          
    
        
        J =  (1/N)*(sum(Jtemp(k-N+1:k)));
        
        Jsave(k) = J;
        
        end
     
     
    
    

    %Pko=Pk;  % for next loop
    
   % Pkosave(:,:,k) = Pko;
    
      
     
     
     % Pkminus=A*(Pko) *A' + Q;   
   end
    
 
         
    
     
 
     
     
 
   
    
    %get P from inverse of information matrix 
    for i = k-N+1:k
   P0temp(:,:,i) = PHIk'*H'*inv(Rk_hat(:,:,i))*H*PHIk; 
    
    end
    
    
    %for i = k-N:k
        
             
        
    sumPtemp = zeros(4,4); %initialise
   % sumPtemp(:,:) = sum(P0temp(:,:,k-N:k));
    for i = k-N+1:k
        
        sumPtemp = sumPtemp + P0temp(:,:,i) ;
        
     %sumPtemp = P0temp(:,:,k) + P0temp(:,:,k-1) + P0temp(:,:,k-2) + P0temp(:,:,k-3) + P0temp(:,:,k-4);
    
     
    end
    P0est(:,:,k) = pinv((1/N)*sumPtemp);   %used pinv here becaue it was going nan
   % end
    
    



%this is the mean and std AFTER the adaptive filter is applied
meanRMS = mean((xhv(1,startadaptive:T) -  xv(1,startadaptive:T)).^2)
stdRMS = std((xhv(1,startadaptive:T) -  xv(1,startadaptive:T)).^2)


figure;


for kk = 1:T
    
    varPko1(kk) = abs(Pko(1,1,kk));
    stdPko1(kk) = sqrt(varPko1(kk));
    
    
    varPko2(kk) = abs(Pko(2,2,kk));
    stdPko2(kk) = sqrt(varPko2(kk));
    
end



figure;
plot((2*stdPko1),'r');
hold
plot((xhv(1,2:T) -  xv(1,2:T))); 



figure;
plot((2*stdPko2),'r');
hold
plot((xhv(2,2:T) -  xv(2,2:T))); 




 
 figure;
plot(Jsave,'g');

tilefigs
pause;
close all


