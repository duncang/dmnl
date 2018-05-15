

feature accel on

runload = 1;

if runload == 1
    
clear all
close all


dT = 1;    %sample time

Sigma_wt = 8;

Sigma_vt = 5


endepoch = 10;

%endepoch = 100;

wt(1) = Sigma_wt*randn(1);
vt(1) = Sigma_vt*randn(1);


endepochHighRate = endepoch*100;

wt(1) = Sigma_wt*randn(1);
vt(1) = Sigma_vt*randn(1);


for i = 2:endepochHighRate
    
    %generate measurement noise, white noise
    
[vt(i),variance,randU] = GaussMarkov_Process(vt(i-1),100000,Sigma_vt,1/dT);          


[wt(i),variance,randU] = GaussMarkov_Process(wt(i-1),100000,Sigma_wt,1/dT);
           

 
end

vt = vt*(1/dT);
wt = wt*(1/dT);


wn = pi/24; %natural frequency
damp = 0.2; %damping ratio



%true process model

xt(1:3,1) = [0, 0.5, 0]'; %initial P, V, Accel
z(1:3,1) = [1, 0, 0]';
%generate truth
for i = 2:endepochHighRate

        
xtdot(1:3,i) = [ 0 1 0 ; 0 0 1; 0 -wn^2 -2*damp*wn;]*xt(1:3,i-1) + [ 0 0 vt(i)]';
   

xt(1:3,i) = xt(1:3,i-1) + xtdot(1:3,i)*dT;

%measurements

z(i) = [1 0 0]*xt(1:3,i) + wt(i);
        
    
end


                          
end %if runload == 1
  

  
  %number of particles
  M = 2000;
  
                          
 %Initialise filter
    
 %state
   %x1_0(1:2,1) = [1, 0];
  
  %x1_0(1:2,1) = [xt(1,1) , xt(2,1)];
   
   for m = 1:M
                
      
     x1_0(1:2,m) = [xt(1,1), xt(2,1)]';
     
     x1p(1,1:2,m) = [xt(1,1), xt(2,1)]';
 
   end
  
  
  %Use SISR algorithm
  
 
  %generate particles
   W_1pf(1,1) = 20^2;
      W_1pf(2,2) = 20^2;
      
      
      R_1pf(1,1) = 20^2;

  
  %driving noise for process    
  xm(1,1,1:M) = sqrt(W_1pf(1,1) )*randn(M,1);    %position  (i, variable, particle number)
  xm(1,2,1:M) = sqrt(W_1pf(2,2))*randn(M,1);     %velocity
  
  
  %driving noise for measurements
  zm(1,1,1:M) = sqrt(R_1pf(1,1) )*randn(M,1);    %position  (i, variable, particle number)
 % zm(1,2,1:M) = sqrt(R_1pf(2,2))*randn(M,1);     %velocity
  
  
  %not sure if this is right
  %generate particles at i = 1
%   for m = 1:M
%       
%   x1p(1,1:2,m) = [1, dT; 0, 1]*x1_0(1:2,m) + [(dT^2)/2; dT].*xm(1,1:2,m)' ;
%   
%   
%   end

for i = 2:endepochHighRate
  
  xm(i,1,1:M) = sqrt(W_1pf(1,1) )*randn(M,1);    %position  (i, variable, particle number)
  xm(i,2,1:M) = sqrt(W_1pf(2,2))*randn(M,1);     %velocity
 
  
  %driving noise for measurements
  zm(i,1,1:M) = sqrt(R_1pf(1,1) )*randn(M,1);    %position  (i, variable, particle number)  

 
end

 %wkstar(1,1,1:M) = rand(M,1);  %sample from uniform distribution
 
 
 
 %initialise weights, all assigned with equal weighting ie uniform
 %distribution? 
 for m = 1:M
 
 wkstar(1,m) = 1/M; 
 
 end

  
  for i = 2:endepochHighRate
        
   
      
  M_save(i) = M;  %save copy of number of particles    
      
      
      %for model 1      
      
      %predict x1, apriori state estimate
               
  for m = 1:M
      
  %x1p(i,1:2,m) = [1, dT; 0, 1]*x1p(i-1,1:2,m)' + [(dT^2)/2; dT].*xm(1,1:2,m)' ;   %apriori state estimate
  
  %x1p(i,1:2,m) = [1, dT; 0, 1]*x1p(i-1,1:2,m)' + [(dT^2)/2; dT].*xm(i,1:2,m)' ;   %apriori state estimate using different xm each time
  
  x1p(i,1:2,m) = [1, dT; 0, 1]*x1p(i-1,1:2,m)' + xm(i,1:2,m)' ;   %apriori state estimate using different xm each time
  
  %get copy of state prediction
  x1p_minus_save(i,:,m) = x1p(i,1:2,m) ;
  end
           

      

      for m = 1:M
      
     % zp(i,m) = [1 0]*x1p(i,1:2,m)' + zm(1,1,m);     %run apriori state estimate through measurements 
      
        zp(i,m) = [1 0]*x1p(i,1:2,m)' + zm(i,1,m);   % try with zm difference each epoch
      end      
      
      %normalized importance ratio
      
      %calculate importance weights
      for m = 1:M      
          
      %qptemp(m) = 1/(sqrt(R_1pf)*2*pi)*exp( -0.5*(R_1pf^-1)*(z(i)-zp(i,m))^2);       %this formula is simply the gauss pdf equation with variance R_1 and z measurements zp, 
     
      qptemp(m) = normpdf(z(i),zp(i,m),sqrt(R_1pf));
      
        %qptemp3(m) = normpdf(X,MU,SIGMA)
        
        
      %qptemp2(m) = (z(i)-zp(i,m)) + sqrt(R_1pf)*randn(1);
      
      end
      
      
     % qp(i,:) = qptemp./(sum(qptemp));   %the sum of this is 1. 
            
      
      for m = 1:M
          times(m) = qptemp(m)*x1p(i,1,m);
          times2(m) = qptemp(m)*x1p(i,2,m);
      end
      
      
      %compute weights
     for m = 1:M
         
        wkstar(i,m) = wkstar(i-1,m)*qptemp(m);
        
     end
     
     %normalise weights
     wk(i,:) = wkstar(i,:)/sum(wkstar(i,:));
    
     
     %compute estimate
     
     for m = 1:M
         
         summer(m) = wk(i,m)*x1p(i,1,m);

         summer2(m) = wk(i,m)*x1p(i,2,m);
       
     end
     
     %I THINK A PROBLEM IS THAT I HAVENT CONSIDERED THE WHOLE VECTOR< IE
     %IVE DONE PV SEPARATE< SO THE V ACCURACY IS BAD
     
     
     %x1totalest(i) = sum(wk(i,:)*x1p(i,1,:));
          
     x1totalest(i) = sum(summer); 
    % x1totalvar(i) = var(summer);
     x1totalest2(i) = sum(summer2);   
     
     %calculate variance, from wikipedia
     
     for m = 1:M
         
%       vartemp(m) = ((wk(i,m) - x1totalest(i))^2)*x1p(i,1,m); 
%          
%       vartemp2(m) = ((wk(i,m) - x1totalest2(i))^2)*x1p(i,2,m);
      
                 
%       vartemp(m) = ((summer(m) - x1totalest(i))^2)*x1p(i,1,m); 
%          
%       vartemp2(m) = ((summer2(m) - x1totalest2(i))^2)*x1p(i,2,m); 
      
      
      
     vartemp(m) = wk(i,m)*(x1p(i,1,m)-x1totalest(i))*(x1p(i,1,m)-x1totalest(i));
         
      vartemp2(m) = wk(i,m)*(x1p(i,2,m)-x1totalest2(i))*(x1p(i,2,m)-x1totalest2(i));
      
     end     
     
     x1totalvar(i) = sum(vartemp);
      x1totalvar2(i) = sum(vartemp2);
         
     
     
     %velocity
    
     %x1totalvar2(i) = var(summer2);    
      
     
     Neff(i) = 1/(sum(wk(i,:).^2));    
     
     
     if Neff(i) < 150
         
         %resample from gauss dist
       %  hi = 44;
         
      
       
       %initalise CDF 
      cdf1(1) = 0;
       wkstar(i,1) = 1/M;
       
      for m = 2:M          
          
          %construct cdf
          
          cdf1(m) = cdf1(m-1)+wk(i,m); 
          
      end
      
      %sample from uniform dist on interval [0, 1/M]
      a = 0; b = 1/M;
      u1(1) = a + (b-a) * rand(1);      
   
      m = 2; %start at bottom of cdf
      
      for j = 1:M
          
          um(j) = u1(1) + (1/M)*(j-1);
          
          while um(j) > cdf1(m)  %plot cdf, is a straight line
              m = m+1;
              
              if m >= M
                  m = M;
                  break;
              end
              
              
          end
         
          %assign sample
          x1p(i,1,j) =  x1p(i,1,m) ;
          x1p(i,2,j) =  x1p(i,2,m) ;
          
          %assign weight
          wkstar(i,j) = 1/M;
          
      end
          
          
          
      
             %reset weights
          %wkstar(i,m) = 1/M;
      
      
       
%        
%        for m = 1:M      
%           
% 
%       
% %       x1p(i,1,m) = x1totalest(i) + sqrt(x1totalvar(i))*randn(1);   
% %       x1p(i,2,m) =   x1totalest2(i)   + sqrt(x1totalvar2(i))*randn(1);    
%       
%        %this isn't exactly right my variances here aren't from the cdf of
%        %the weights but just my initial values
%          
%       x1p(i,1,m) = x1totalest(i) + sqrt(W_1pf(1,1) )*randn(1);   
%       x1p(i,2,m) =   x1totalest2(i)   + sqrt(W_1pf(2,2) )*randn(1);    
%       
%             
%             
%       %reset weights
%       wkstar(i,m) = 1/M;  %sum of all these will be 1
%             
%       
%       
%       
%       end
%          
         
         
         
     end
     
         
     
     
      
%      testweight1(i,:) = times./(sum(qptemp.*times));
%      
%      testweight2(i,:) = times2./(sum(qptemp.*times2));
    
     
      %this must be for direct method PF not the SIS version
      %update states, and resample
%       up = rand(M+1,1);    %sample from uniform disttibution
%       tp = -log(up);   %this is a property of the uniform distribution, I THINK??
%       
%       Tp = cumsum(tp);
%       Qp = cumsum(qp(i,:))';
%       
%       k = 1;
      
%       for m = 1:M      
%      
%       %SIR 
%       %this resamples to discard low weights, so the size of the particles
%       %will be less
%       
% %       if m == 255
% %           hi = 2
% %       end
% 
% 
% 
% 
%   
% 
% 
% 
% 
%         %calculate degenerency
%       
%       %calculate N effective approximation, 
%       
%       
%       
%      % if Neff(i) < 50
% 
%      
% 
% %       
% %       if Qp(m,1)*Tp(M,1)> Tp(k,1)     
% %       
% %           x1ptemp(i,1:2,k) = x1p(i,1:2,m);
% %           
% %           x1p(i,1:2,k) = x1p(i,1:2,m);  %save the particles with high weight
% %           k = k+1;
% %       end  
% %       
% 
%      % end
%       
%       
%       
%       end      
      
      %without resampling ie SIS algorithm
      
%       for m = 1:M
%           
%       x1p(i,1,m) = qp(i,m)*x1p(i,1,m);   %update state estimate
%       x1p(i,2,m) = qp(i,m)*x1p(i,2,m);
%           
% %       x1p(i,1,m) = qp(i,m)*x1p(i,1,m) + x1p(i,1,m);   %update state estimate
% %       x1p(i,2,m) = qp(i,m)*x1p(i,2,m) + x1p(i,2,m); ;
%       
% %       x1p(i,1,m) = testweight1(i,m);   %update state estimate
% %       x1p(i,2,m) = testweight2(i,m);
%           
%       
%      % x1ptest(i,1,m) = qp(i,m)*x1p(i,1,m);
%           
%      
%          
% %       x1p(i,1,m) = (qptemp(m)*x1p(i,1,m))./(sum(qptemp(m)*x1p(i,1,m)));
% %       x1p(i,2,m) = (qptemp(m)*x1p(i,2,m))./(sum(qptemp(m)*x1p(i,2,m)));
%        
%       %the new number of samples, M    
%       %M = k; 
%       end   
      %Neff(i) = 1/(sum(qp(i,:).^2));    
      
      
      
  end
          
      
      
  
for i = 1:endepochHighRate
  x1total(i) = mean(x1p(i,1,:));
end



figure;
hold;
plot(x1totalest(1:endepochHighRate));
plot(xt(1,1:endepochHighRate),'r');






figure;
hold;
plot(x1totalest2(1:endepochHighRate));
plot(xt(2,1:endepochHighRate),'r');

tilefigs;
pause;

for jj = 1:50
    
figure;   
hist(x1p_minus_save(jj,1,:));
figure;
hist(x1p(jj,1,:));
tilefigs;
pause;
close all

end
      
      
      
%       
%       
%       
%       
%       
% 
% x1in(1:2,1) = x_hat_out(1:2);
% 
% 
% 
%     
% x1in(1:2,2) = [1, dT; 0, 1]*x1in(1:2,1);
%       
% 
% 
% %for model 2   
% 
%       %predict x21, apriori state estimate
% x2in(1:3,2) = [1, dT,(dT^2)/2;   0, 1, dT;  0 , 0 , 1;]*x2in(1:3,1);           
%        
% 
% PHI_C = [ 1 dT 0 0 0;
%           0 1 0 0 0;
%              0 0 1 dT (dT^2)/2;
%              0 0 0 1 dT;
%                  0 0 0 0 1;];
%              
%                  
%  PHI_1 =  [1 dT; 0 1]; 
%                  
%  PHI_2 =  [1, dT,(dT^2)/2;   0, 1, dT;  0 , 0 , 1;];
%                  
%   
% %  
% %  W_1(1,1) = 2*Sigma_vt^2;
% %   W_1(2,2) = 2*Sigma_vt^2;
% % 
% %    W_2(1,1) = 2*Sigma_vt^2;
% %     W_2(2,2) = 2*Sigma_vt^2;
% %      W_2(3,3) = 2*Sigma_vt^2;
% 
%   W_1(1,1) = 100*0.3*1^2;
%   W_1(2,2) = 100*0.3*1^2;
% 
%    W_2(1,1) = 100*0.3*1^2;
%     W_2(2,2) = 100*0.3*1^2;
%      W_2(3,3) =100*2*Sigma_vt^2;
% 
%  
%      G_1 = zeros(2);
%      G_2 = zeros(3);
%      
%     G_1(1,1) = (dT^2)/2;     
%     G_1(2,2) = dT;
%     
%      G_2(1,1) = (dT^3)/6;     
%     G_2(2,2) = (dT^2)/2;
%      G_2(3,3) = dT;
%                  
%      
%      
%      
%  
%  
%  Q_1 = PHI_1*G_1*W_1*G_1'*PHI_1'*dT;
%  Q_2 = PHI_2*G_2*W_2*G_2'*PHI_2'*dT;
%  
%  
% %   
% %  Q_1 = PHI_1*W_1*PHI_1'*dT;
% %  Q_2 = PHI_2*W_2*PHI_2'*dT;
% %  
%  
%  CorrVal = 1;  
%   
%  
%   Qk1k2 = [0 0 0;
%           0 0 0;];
%       
%   
%  Qk1k2(1,1) = CorrVal*sqrt(abs(Q_1(1,1)))*sqrt(abs(Q_2(1,1)));
%  Qk1k2(2,2) = CorrVal*sqrt(abs(Q_1(2,2)))*sqrt(abs(Q_2(2,2)));
%           
%  
%  
%  Qk1k2 = [0 0 0;
%           0 0 0;];
%      
%   
%  
%            
%   Q_C = [Q_1, Qk1k2;
%           Qk1k2', Q_2;];  
%         
%   
%           
%  
% %Q_C = eye(5,5);    
% % Q_C(1,1) = 0.5;
% % Q_C(2,2) = 0.05;
% % Q_C(3,3) = 0.5;
% % Q_C(4,4) = 0.05;
% % Q_C(5,5) = 0.005;
% 
% 
% % Q_C(1,1) = 0.5;
% % Q_C(2,2) = 0.05;
% % Q_C(3,3) = 0.3;
% % Q_C(4,4) = 0.025;
% % Q_C(5,5) = 0.005;
% 
% 
%                           
% P_minus_out_CTemp = PHI_C * P_minus_in_CTemp * PHI_C' + Q_C;
% 
% 
% % end  %end k
% 
% P_minus_C = P_minus_out_CTemp;
%              
% x_hat_minus = [x1in(1:2,2);x2in(1:3,2);];
%  
% 
% 
% 
% 
% H_C = [1 0 0 0 0;
%        0 0 1 0 0;];
%    
%        
%        R = eye(1);
%        
%        
%        
% %        R_C = [R, 0, R;
% %               0, 0, 0;
% %               R, 0, R;]; 
%             
%                 
%        R_C = [R, R;
%               R, R;]; 
%       
% %              
% %        R_C = [R, 0;
% %               0, R;]; 
% % %           
% 
%      NumberStates  = 2+3;
% 
% %apriori state covariance 
%               
%         
%        z_C = [z(i), z(i)]';      
%                                                  
%                 
% 
% 
%     V_C = H_C * P_minus_C * H_C' + R_C;  %this is the innovations  covariance
%         
%     
%     K_C = P_minus_C * H_C' * pinv(V_C);    
%     
%     
%     v = z_C - H_C * x_hat_minus;           
%     
%     x_hat_out = x_hat_minus + K_C * (v);
%   
%     
% 
%     %x_hat_out_C = K_C *(z_C);    
% 
%    P_out_C = (eye(NumberStates) - K_C*H_C)*P_minus_C*(eye(NumberStates)-K_C*H_C)' + K_C*R_C*K_C';       
%     
%     
%    x_hat_out_save(:,i) = x_hat_out; 
%    P_out_C_save(:,:,i) = P_out_C;
%       
%       
%   end
%   
%       
%   %do some plots  
%   %plot positions  
% 
% posstd1 = sqrt(P_out_C_save(1,1,:));
% posstd2 = sqrt(P_out_C_save(3,3,:));
% 
% velstd1 = sqrt(P_out_C_save(2,2,:));
% velstd2 = sqrt(P_out_C_save(4,4,:));
% 
% 
% accstd2 = sqrt(P_out_C_save(5,5,:));
% 
% figure();
% plot(x_hat_out_save(1,:),'b-','LineWidth',2); title 'Position Error(m)';xlabel('Time');
% hold;
% plot(x_hat_out_save(3,:),'r-','LineWidth',2); 
% 
% % plot(xt(1,:) + 2*posstd1(:),'b--');
% % plot(xt(1,:) - 2*posstd1(:),'b--');
% % 
% % plot(xt(1,:) + 2*posstd2(:),'r--');
% % plot(xt(1,:) - 2*posstd2(:),'r--');
% 
% %plot truth
% plot(xt(1,:),'g');
% 
% 
% 
% poserr1 = xt(1,:) - x_hat_out_save(1,:);
% poserr2 = xt(1,:) - x_hat_out_save(3,:);
% 
% velerr1 = xt(2,:) - x_hat_out_save(2,:);
% velerr2 = xt(2,:) - x_hat_out_save(4,:);
% 
% accerr2 = xt(3,:) - x_hat_out_save(5,:);
% 
% 
% figure();
% hold;
% plot(abs(poserr1),'b-','LineWidth',2); title 'Position Error(m)';xlabel('Time');
% plot(abs(poserr2),'r-','LineWidth',2); 
% 
% plot(abs(velerr1),'b-','LineWidth',1); 
% plot(abs(velerr2),'r-','LineWidth',1); 
% plot(abs(accerr2),'k-','LineWidth',1); 
% 
% plot(0 + 2*posstd1(:),'b--','LineWidth',2);
% %plot(0 - 2*posstd1(:),'b--');
% 
% plot(0 + 2*posstd2(:),'r--','LineWidth',2);
% %plot(0 - 2*posstd2(:),'r--');
% 
% 
% plot(0 + 2*velstd1(:),'b--','LineWidth',1);
% %plot(0 - 2*posstd1(:),'b--');
% 
% plot(0 + 2*velstd2(:),'r--','LineWidth',1);
%      
% 
% plot(0 + 2*accstd2(:),'k--','LineWidth',1);
% tilefigs;
% pause;
% %close all;
% 
% 
% 
% % 
% % 
% % 
% % %plot  P and V like Julier did
% % 
% % 
% % posstd1 = sqrt(P_out_C_save(1,1,:));
% % posstd2 = sqrt(P_out_C_save(3,3,:));
% % 
% % figure();
% % plot(x_hat_out_save(1,:),'b-','LineWidth',2); title 'Position Error(m)';xlabel('Time');
% % hold;
% % plot(x_hat_out_save(3,:),'r-','LineWidth',2); 
% % 
% 
% 
% 
% %calculate correlations
% 
% 
% %calculate correlation coefficients
% 
% for i = 1:endepoch
% 
% corrcoeffPos(i) =  P_out_C_save(1,3,i)./(  sqrt(P_out_C_save(1,1,i)).*sqrt(P_out_C_save(3,3,i)));
% corrcoeffVel(i) =  P_out_C_save(2,4,i)./(  sqrt(P_out_C_save(2,2,i)).*sqrt(P_out_C_save(4,4,i)));
% 
% end
% 
% 
% 
% 

