function [SVPos, PRMeasured_Simulated,PRMeasured_SimulatedCSAC, GeometricRange, SVVelocity] = GARDSim_GenerateMeasurementsSimulated(i,GPStimeAMicro,UserPos,N,N_Prev, Sats, SV_X_Data, SV_Y_Data, SV_Z_Data, SV_T_Data, SV_Xvel_Data, SV_Yvel_Data, SV_Zvel_Data, SV_Tvel_Data, ALPHA, BETA, dTposCSAC);


global GPS_PI OMEGAedot mu Earthradius Speedoflight c F L1_f L2_f gamma L1_Wavelength;




%Version 1.00
%Troy Bruggemann 11 August 2005

%This function Generates Matrix of satellite positions and Predicted and
%Observed Pseudorange Vectors from given input as described below:
       

%INPUTS 
% UserPos - [XPos,YPos,ZPos,dTPos] Estimated User Navigation State Vector(m)
% N - number of observations to use in solution (number of satellites)
% PRMeasured - Vector of Measured Pseudoranges [1..N]
% SVPos - Matrix of Satellite poisitions [1..N][Xs,Ys,Zs,dTs] (m)
%========================================================================
% OUTPUTS 
%SolutionVec = [X,Y,Z,dt];
% X - estimated position ECEF (m)
% Y - estimated position ECEF (m)
% Z - estimated position ECEF (m)
% dt - estimated receiver clock bias (m)
% 
%VarSolutionVec = [var_x,var_y,var_z,var_dt];
% var_x - estimated X position variance (m)
% var_y - estimated X position variance (m)
% var_z - estimated X position variance (m)
% var_dt - estimated X position variance (m)
% NumIterations - Number of iterations used for the solution
% ResidualVector - [1..N] Vector of pseudorange residuals (m)
% M - Least squares design matrix [1..M][1..4]



% Speedoflight = 2.99792458e8;


%--------------------------------------------------------------------------
        %Generate Matrix of satellite positions and Predicted and Observed Pseudorange Vectors.
        %---------------------------------------------------------------------
        for PRNCounter = 1:32   %for all satellites in constellation

            for k = 1:N   %run N times
                if PRNCounter == Sats(k);  %Find the satellite we want to use

                    SV(i,k,1:4) = [SV_X_Data(i,PRNCounter), SV_Y_Data(i,PRNCounter), SV_Z_Data(i,PRNCounter), SV_T_Data(i,PRNCounter)*Speedoflight ] ;
                    %Predicted Pseudoranges for N satellites
                    
                    SatPos(1) =   SV(i,k,1);
                    SatPos(2) =   SV(i,k,2);
                    SatPos(3) =   SV(i,k,3);
                    SatPos(4) =   SV(i,k,4);

                    
                    SVVel(i:i+5,k,1:4) = [SV_Xvel_Data(i:i+5,PRNCounter), SV_Yvel_Data(i:i+5,PRNCounter), SV_Zvel_Data(i:i+5,PRNCounter), SV_Tvel_Data(i:i+5,PRNCounter)*Speedoflight ] ;
                    
                      
                    SatVel(1) =   SVVel(i+2,k,1);
                    SatVel(2) =   SVVel(i+2,k,2);
                    SatVel(3) =   SVVel(i+2,k,3);
                    SatVel(4) =   SVVel(i+2,k,4); 
                     
                                      
                    %UserPos(4) is modelled clock bias..only varies
                    %between epoch not between each measurement
                    %----------------------------------------------------
                    %Form Simulated (Simulated) observed PR's
                    
                    
                    
                      %To add ionospheric delay to Simulated PRs
                    IonoDelay = ionomodel(GPStimeAMicro, UserPos(1:3), [SV(i,k,1) SV(i,k,2) SV(i,k,3)] , ALPHA, BETA);

                    
                    time_step = 1;
                    
                   
                    
                    [PRMeasured_Simulated(k),PRMeasured_SimulatedCSAC(k),GeometricRange(k)] = GARDSim_PRPredictSimulated(SatPos,UserPos,IonoDelay, dTposCSAC);

                                 

                    
                    
                    %Add various errors onto the Simulated PR
                                                 
                    %PRMeasured_Simulated(k) = PRMeasured_Simulated(k) + IonoDelay;

                    %have to add this to the previous value, for ICP, and
                    %then add the wavelength*Namb because it is time
                    %independant, unless cycle slips
                    % CPMeasPredict(k,i) = CPMeasPredict(k) + NAmb

                    %Introduce cycle slips?

                    
                      %Put error on PRmeasurement - Test

                    %SV to put error on is SV 6
                    
%                    
                     %if i == 1 | i == 2 |i == 3 | i == 4 |i == 5 | i == 6 | i == 7
%                      
%                      if i >1 & i < 300
%                        if PRNCounter == 6 
%                           % PRMeasured_Simulated(k) = PRMeasured_Simulated(k) + 8000;
%                         
%                            %add ramp error
%                            PRMeasured_Simulated(k) = PRMeasured_Simulated(k) + 5*i;
%                             
%                        end
%                       end
%                  


%                      if i == 1 | i == 2
%                        if PRNCounter == 6 
%                            PRMeasured_Simulated(k) = PRMeasured_Simulated(k) + 400000;
%                            
%                        end
%                       end
%                     
%                     
                    %Error put on for 835 Assignment - error is already in
                    %835 rinex files

%                     if i == 2
%                       if PRNCounter == 1 
%                           PRMeasured_Simulated(k) = PRMeasured_Simulated(k) + 400;
%                           PRMeasured_Observed(k) = PRMeasured_Observed(k) + 400;
% 
%                           % PRRangeRate(k) = PRMeasured_Observed(k) + ;
%                       end
%                      end                                     

                  
                    %generate matrix For least square solution
                    SVPos(k,1) =  SV(i,k,1);
                    SVPos(k,2) =  SV(i,k,2);
                    SVPos(k,3) =  SV(i,k,3);
                    SVPos(k,4) =  SV(i,k,4);
                    
                    SVVelocity(k,1) =  SVVel(i,k,1);
                    SVVelocity(k,2) =  SVVel(i,k,2);
                    SVVelocity(k,3) =  SVVel(i,k,3);
                    SVVelocity(k,4) =  SVVel(i,k,4);


                end
            end

        end
        
        
        
        
        
        
        
        
     