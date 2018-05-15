% 
% 


%Read in pseudoranges
% 

Speedoflight = 2.99792458e8; 


Filename = 'GPSl0770.txt'



[GPSTime_Week, GPSTime_Sec,NumberRinexObsTypes,ApproxPos, DATA_STRUCT] = ReadRinexGRS(Filename)




%calculate satellite position





%True user position




%least squares best fit


Signal = DATA_STRUCT.C1(8,:)
 [CurveFit, BestFit] = LeastSquaresBestFit(Signal, 2, 1);
 
%  
 SP3FILE = 'RAPID';
 %read SP3 file
 
 
 %for i = 1:100
 
  [NumberSVs, VehicleIDs, EpochNumber, GPSTime, SV_X_Data, SV_Y_Data, SV_Z_Data, SV_T_Data] = readsp3(SP3FILE);
  
  
 %end
 
 
 [GPS_WEEK, GPSsecs] = ftime([08 03 17 00 30  0.0000000] )
 
  
 

 
 
 %TRUE 
 
 userpos = [-5046773.357, 2568446.084, -2925289.017];
 
 
 
 
 PosTruthSub = [userpos, 0];
 %calculate receiver clock bias
 
 PRNvec = [10, 3, 25, 11, 28, 23, 27, 19, 13];
 
 
 
  for k = 1:9
 data = SV_X_Data(1:50,PRNvec(k));
x = 1:50;
y = data;
xx = 1:0.0011:50;
SVX(k,:) = spline(x,y,xx);


data = SV_Y_Data(1:50,PRNvec(k));
x = 1:50;
y = data;
xx = 1:0.0011:50;
SVY(k,:)  = spline(x,y,xx);


data = SV_Z_Data(1:50,PRNvec(k));
x = 1:50;
y = data;
xx = 1:0.0011:50;
SVZ(k,:)  = spline(x,y,xx);



data = SV_T_Data(1:50,PRNvec(k));
x = 1:50;
y = data;
xx = 1:0.0011:50;
SVT(k,:)  = spline(x,y,xx);
 
  end
  
 GPSConstants;
 
 for i = 1:1000
     
          
     NSub(i) = 9;
     
     for k = 1:9
     PRMeasured_SimulatedLSQSub(k) = DATA_STRUCT.C1(PRNvec(k),i) ;
     

     
     SatPosLSQSub(k,1:4) = [SVX(k,i), SVY(k,i),SVZ(k,i),SVT(k,i)*Speedoflight];
     
     end
 
     
    [GPSPosLSQSub VarSolutionVec_Sub  NumIterations_Sub  ResVec_Sub  M_Sub LSQ_FailSub limit_Sub DOP_Sub(1:5)] = GARD_LSQ(PosTruthSub,NSub(i),PRMeasured_SimulatedLSQSub,SatPosLSQSub);

  
        XPos(i) = GPSPosLSQSub(1);
            YPos(i) = GPSPosLSQSub(2);
                ZPos(i) = GPSPosLSQSub(3);
                    dtPos(i) = GPSPosLSQSub(4);
        
 end
 
 
 
 
 
 
 for i = 1:1000
     
     PRtrue(i) = sqrt(   (SVX(3,i) - userpos(1))^2 + (SVY(3,i) - userpos(2))^2 + (SVZ(3,i) - userpos(3))^2);
     
 end
 
 
 
 PR_OEMV_25 = DATA_STRUCT.C1(25,1:1000) -dtPos(1:1000)  + Speedoflight*SVT(3,1:1000);
 
 
 
 
 
 
 
 
 
 