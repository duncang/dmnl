function [SVPos, PRMeasured_Observed, PRMeasured_Simulated, PRRateMeasured_Simulated,ICPMeasured_Observed, ICPMeasured_Simulated, CPMeasured_Simulated, NAmb, GeometricRange, SVVelocity, PRRangeRate] = GARDSim_GenerateMeasurementsObserved(i,GPStimeAMicro,UserPos,N, Sats, C1_PRNAMicro, L1_PRNAMicro, SV_X_Data, SV_Y_Data, SV_Z_Data, SV_T_Data, SV_Xvel_Data, SV_Yvel_Data, SV_Zvel_Data, SV_Tvel_Data, ALPHA, BETA,PRMeas_previous, ICPMeas_previous,CPMeas_previous,GeometricRange_Previous,EstimatedClockBias);


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
                    
                    [PRMeasured_Simulated(k),PRRateMeasured_Simulated(k),ICPMeasured_Simulated(k),CPMeasured_Simulated(k),NAmb(k),GeometricRange(k)] = GARDSim_PRPredict(SatPos,UserPos,IonoDelay,time_step,PRMeas_previous(k), ICPMeas_previous(k), CPMeas_previous(k),GeometricRange_Previous(k),EstimatedClockBias);

                    
                  
                    %function [PRMeasPred, PRRMeasPredict, CPMeasPredict,NAmb, GeometricRange] = GARDSim_PRPredict(SatPos,UserPos, IonoDelay, time_step, PRMeas_previous, CPMeas_previous,GeometricRange_Previous);

                    
                    %Add various errors onto the Simulated PR
                                                 
                    %PRMeasured_Simulated(k) = PRMeasured_Simulated(k) + IonoDelay;

                    %have to add this to the previous value, for ICP, and
                    %then add the wavelength*Namb because it is time
                    %independant, unless cycle slips
                    % CPMeasPredict(k,i) = CPMeasPredict(k) + NAmb

                    %Introduce cycle slips?

                    %Observed from rinex file
                    PRMeasured_Observed(k) = C1_PRNAMicro(PRNCounter,i);
                                       
                    ICPMeasured_Observed(k) = L1_PRNAMicro(PRNCounter,i);
                    
                    PRRangeRate(k) = GPSVelocityVector(i+1,UserPos,N, Sats, C1_PRNAMicro(PRNCounter,i+1-1:i+1+1), L1_PRNAMicro(PRNCounter,i+1-1:i+1+1), SV_X_Data, SV_Y_Data, SV_Z_Data, SV_T_Data);
% 
%                     %                     SV(i,k,1:4) = [SV_X_Data(i,PRNCounter), SV_Y_Data(i,PRNCounter), SV_Z_Data(i,PRNCounter), SV_T_Data(i,PRNCounter)*Speedoflight ] ;
%                    
%                     
                    PRRangeRate(k) = PRRangeRate(k)*L1_Wavelength; %in metres per second
%                     
                    %-------------------------------------------------
                    %Smooth pseudoranges with L1
                    %Carrier phase smoothing , Hatch implementation
                    %                     if i == 1
                    %                         PRMeasured_Observed(k) = C1_PRNAMicro(PRNCounter,i)
                    %                         PreviousSmoothed(k) = PRMeasured_Observed(k);
                    %                         R_Ex = PRMeasured_Observed(k);
                    %
                    %                     else
                    %
                    %                     R_Ex = PreviousSmoothedGlobal(i-1,k) + (L1_PRNAMicro(PRNCounter,i) - L1_PRNAMicro(PRNCounter,i-1));
                    %
                    %                     %smoothed
                    %                     PRMeasured_Observed(k) = 0.5*(C1_PRNAMicro(PRNCounter,i) + R_Ex);
                    %
                    %
                    %                     PreviousSmoothed(k) = PRMeasured_Observed(k);
                    %                     end
                    %---------------------------------------------------

                      %Put error on PRmeasurement - Test

                    %SV to put error on is SV 6
                    
%                    
%                      if i == 1 | i == 2 |i == 3 | i == 4 |i == 5 | i == 6 | i == 7
%                        if PRNCounter == 6 
%                            PRMeasured_Simulated(k) = PRMeasured_Simulated(k) + 4000;
%                            PRMeasured_Observed(k) = PRMeasured_Observed(k) + 4000;
%                             
%                        end
%                       end
                 


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
        
        
        
        
        
        
        
        
     