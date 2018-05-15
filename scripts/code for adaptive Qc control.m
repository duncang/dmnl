
    
    %Adaptively Control Q_C
    
    
     
     %Now adaptively control Q_C for next time around
     
     %note, that Q_C is calculated at 100Hz, but this adaptive bit is done
     %at 1 Hz.
     
     
       Q_Cuse = 0; %if this is 1, then use the adaptive value of Q_C  
         if ControlQ_C == 1 && i>=10
         
         %do it over the last 5 epochs. NOTE that if control Q_C is
         %enabled i dont start it until 10 epochs have passed.
         
         numepochs = 5;
         for k = 1:numepochs
             
             deltaxsquared(:,:,k) = x_hat_save(:,i-k+1)*x_hat_save(:,i-k+1)';
             
         end
         
         for k = 1:34
             for j = 1:34
                 Cdelta_xk(k,j) = mean(deltaxsquared(k,j,:));
             end
         end
         
       
         
         Q_Ctemp = Cdelta_xk +  P_save(:,:,i) -  PHI_C*P_save(:,:,i-1)*PHI_C';
    
         
         Q_Cuse = 1;
         
                   
        QK1K2actual =  Q_C(1:17,18:34);
         QK1K2actualprevious =  Q_Ctemp(1:17,18:34);
         
         Q_1 = Q_Ctemp(1:17,1:17);
           Q_2 = Q_Ctemp(18:34,18:34);
         
         %calculate correlation 
         
                          
         
                Q_C = Q_Ctemp;
                
                Q_C(1:17, 18:34) = zeros(17,17);
                
                Q_C(18:34, 1:17) = zeros(17,17);
                
%                 [Q_1, zeros(17,17);
%                    zeros(17,17), Q_2;];
                   
         
%          Q_C  = Q_1
         
         
             for jup = 1:17
                 
                % QK1K2(jup,jup) = CorrVal*sqrt(abs(Q_INS(jup,jup)))*sqrt(abs(Q_M(jup,jup)));
                 
                 %CorrVal(jup) = QK1K2actualprevious(jup,jup)/sqrt(abs(Q_1(jup,jup)))*sqrt(abs(Q_2(jup,jup)));
                 
                              
                 
             end
               %mean of CorrVal is about 0.2
         
               
         end
         