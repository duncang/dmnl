function [FinalSolution_ObservedFinal, RAIM_ALERT_FDE_ObsFinal, GoodSatsSats_FDEFinal, Sats_FDE_Obs, DOP] = GARDSim_FDE(SimulatedOrObserved_Flag,i,GPStimeAMicro,UserPos,N, Sats,C1_PRNAMicro,L1_PRNAMicro, SV_X_Data, SV_Y_Data, SV_Z_Data, SV_T_Data,SV_Xvel_Data, SV_Yvel_Data, SV_Zvel_Data, SV_Tvel_Data,PFalseALERT,SigmaS,ALERT_Limit, ALPHA,BETA,PRMeas_previous, ICPMeas_previous,CPMeas_previous,GeometricRange_Previous,EstimatedClockBias)
%$Id: GARDSim_FDE.m 1863 2008-07-14 07:02:29Z greerd $
%Troy Bruggemann 11 August 2005

%This function performs fault identification and exclusion 
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



%Do RAIM FDE solution for Observed PRs

%FDE taken from Parkinson Volume 2

%if failure...go through all the available subsets...
% if 6 satellites, can only detect one error at a time.

% if error occurs on more than one subset then have to assume that more
% than one satellite failure has occurred.


%Do FDE solution
%if failure:


        %Do FDE - ONLY DETECTS ONE SATELLITE FAILURE AT A TIME AT THE MOMENT

        %Form N-1 Subsets
        IsAValue = 0;

        N_FDE = N-1;  %Detect one error

        for i_FDE = 1:N  %number of times have to calculate subset solution.

            %make a new vector of satellites to use in the solution,ie a new Sats(k)
            %circular shift and drop last element

            B = circshift(Sats,[0 i_FDE]);
            % SatsFDECopy(i,i_FDE,1:N_FDE) = B(1:N_FDE);  %Get global copy of satellites used for FDE

            Sats_FDE_Obs(i_FDE,:) = B(1:N_FDE);

             %Generate new Matrix of satellite positions and Predicted Pseudoranges

            [SVPos_FDE,PRMeasured_FDE_Observed,PRMeasured_FDE_Simulated,PRRate_FDE_Simulated,ICPMeasured_ObservedFDE,ICPMeasured_SimulatedFDE,CPMeasuredSimulated_FDE,NAmb_FDE, GeometricRange_FDE SVVelocity_FDE PRRangeRate_FDE] = GARDSim_GenerateMeasurements(i,GPStimeAMicro,UserPos,N_FDE, Sats_FDE_Obs(i_FDE,:),C1_PRNAMicro,L1_PRNAMicro, SV_X_Data, SV_Y_Data, SV_Z_Data, SV_T_Data,SV_Xvel_Data, SV_Yvel_Data, SV_Zvel_Data, SV_Tvel_Data, ALPHA, BETA,PRMeas_previous, ICPMeas_previous,CPMeas_previous,GeometricRange_Previous,EstimatedClockBias);
      
            
%             
%          [SVPos,PRMeasured_Observed,PRMeasured_Simulated,PRRate_Simulated,ICPMeasured_Observed,ICPMeasured_Simulated,CPMeasured_Simulated,NAmb, GeometricRange, SVVel, PRRate_Observed] = GARDSim_GenerateMeasurements(i,GPStimeAMicro,UserPos,N, Sats,C1_PRNAMicro,L1_PRNAMicro, SV_X_Data, SV_Y_Data, SV_Z_Data, SV_T_Data, SV_Xvel_Data, SV_Yvel_Data, SV_Zvel_Data, SV_Tvel_Data, ALPHA, BETA,aa,qq,bb,cc,EstimatedClockBias(i));
%         
            
            
            %if should do FDE For observed measurements:
            if SimulatedOrObserved_Flag == 0

                PRMeasured_FDE = PRMeasured_FDE_Observed;
                %otherwise do FDE for simulated measurements:
            elseif SimulatedOrObserved_Flag == 1

                PRMeasured_FDE = PRMeasured_FDE_Simulated;

            end


            [SolutionVec_FDE(i_FDE,1:4) VarSolutionVec_FDE NumIterations_FDE ResVec_FDE M_FDE LSQ_Fail_FDE(i_FDE) limit_FDE DOP] = GARD_LSQ(UserPos,N_FDE,PRMeasured_FDE,SVPos_FDE);
          
            
            if LSQ_Fail_FDE(i_FDE) == 0

                   
%                 size(UserPos)
%                 size(SVPos_FDE)
%                 size(N)
%                 
                
                [HPL_H0(i_FDE) VPL_H0(i_FDE)] = GARD_HPLVPLGRAS(N_FDE,UserPos, SVPos_FDE);
                
                %[BadGeometry_FDE_Obs(i_FDE) RAIM_ALERT_FDE_Obs(i_FDE) SLOPE_Max_FDE_Obs(i_FDE) r_FDE_Obs(i_FDE) rd_FDE(i_FDE) ARP_FDE(i_FDE)] = GARD_RAIMGRAS(N_FDE,PFalseALERT,SigmaS,ALERT_Limit,ResVec_FDE,M_FDE,HPL_H0,VPL_H0);

                 %Using parkinsons method (for 835 assignment)
                
                [BadGeometry_FDE_Obs(i_FDE), RAIM_ALERT_FDE_Obs(i_FDE),  SLOPE_Max_FDE_Obs(i_FDE), r_FDE_Obs(i_FDE), rd_FDE(i_FDE),  ARP_FDE(i_FDE)] = GARD_RAIM(N_FDE,PFalseALERT,SigmaS,ALERT_Limit,ResVec_FDE,M_FDE);
                
                
                
                
                              
                SolutionFDE(i_FDE,1:4) = SolutionVec_FDE(i_FDE,1:4);

                %check if a good solution exists:

                if RAIM_ALERT_FDE_Obs(i_FDE) ~= 100
                              
                              i_FDEGood = i_FDE;  
                              IsAValue = 1;                              
                             
                end
                 
                
            else %If least squares does not converge, still need to fill in the values at this i_FDE value otherwise vector sizes will be wrong
            BadGeometry_FDE_Obs(i_FDE) = 0;
            RAIM_ALERT_FDE_Obs(i_FDE)= 0;
            SLOPE_Max_FDE_Obs(i_FDE)= 0;
            r_FDE_Obs(i_FDE)= 0;
            rd_FDE(i_FDE) = 0;
            ARP_FDE(i_FDE) = 0;
            HPL_H0(i_FDE) = 0;
            VPL_H0(i_FDE) = 0;
            
            SolutionFDE(i_FDE,1:4) = 0;
                

            end %if LSQ_Fail == 0             


        end % for i_FDE = 1:N       
                
        
        %output results
        
        if IsAValue == 1
            
        
         GoodSatsSats_FDEFinal = Sats_FDE_Obs(i_FDEGood,:);
        
         FinalSolution_ObservedFinal(1) = SolutionFDE(i_FDEGood,1);
         FinalSolution_ObservedFinal(2) = SolutionFDE(i_FDEGood,2);
         FinalSolution_ObservedFinal(3) = SolutionFDE(i_FDEGood,3);
         FinalSolution_ObservedFinal(4) = SolutionFDE(i_FDEGood,4);
         
         RAIM_ALERT_FDE_ObsFinal = 0;
                                  
         
      
        elseif IsAValue == 0
            
         GoodSatsSats_FDEFinal = zeros(1,N_FDE);
         
         RAIM_ALERT_FDE_ObsFinal =100;
        
         FinalSolution_ObservedFinal(1) = 0;
         FinalSolution_ObservedFinal(2) = 0;
         FinalSolution_ObservedFinal(3) = 0;
         FinalSolution_ObservedFinal(4) = 0;
                    
       
        end        
        
        %To improve function:
        %Identify number of 'good' solutions obtained
        %calculate whether identification of more than 1 failed satellite
        %can be achieved.        
        
        
%Verification
%a)
%This function has been tested with one error put on one satellite's psuedoranges, the
%satellite was correctly identified and a correct navigation solution still
%outputted.
%b)
%This was tested with all least squares not converging for all subsets, and
%the function correctly output 0,0,0,0 as the solution (no solution). 
%
%c) This function was tested for the case where the least squares did not converge, with the one subset that did not contain the
%failed satellite. The function correctly output 0,0,0,0 as the result

%d) this function doesn't work yet if more than one satellite has failed






   
