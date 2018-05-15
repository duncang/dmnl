function [PRMeasPredict,PR_RateMeasPredict, GeometricRange] = GARDSim_PRPredictSimulated(SatPos, SatVel, UserPos, UserVel, GPSPr_noiseError, GPSPr_Rate_noiseError, PRMeasPredict_Previous);
% function [PRMeasPredict,PR_RateMeasPredict, GeometricRange] =
% GARDSim_PRPredictSimulated(SatPos, SatVel, UserPos, UserVel, 
%                           GPSPr_noiseError, GPSPr_Rate_noiseError, 
%                           PRMeasPredict_Previous);
% Troy Bruggemann 30 June 2005
% $Id: GARDSim_PRPredictSimulated.m 1880 2008-07-15 05:21:10Z n2523710 $
%
%This function generates simulated measured raw Pseudoranges and Carrier Phase and Doppler based on GPS Satellite position in
%ECEF , and given user position in ECEF
%Incorporates the true geometric range to the satellite + additional biases
% %does for Single Epoch Only and does only One Satellite
%
%
%INPUT
%SatPos -  satellite positions in ECEF [X Y Z dT] (m)
%UserPos - User position in ECEF [X Y Z dt] (m)
%
%OUTPUT 
%PRMeasPred - Predicted Pseudorange measurement (m)
%PRRMeasPredict - Predicted Pseudorange Rate measurement (m/s)
%CPMeasPredict - Predicted integrated carrier phase measurement (m)
%NAmb - Integer Ambiguity 
%GeometricRange - true range between satellite and user (m)
% time step is the sample time period in seconds
% PRMeas_previous, CPMeas_previous are the previous pseudorange and carrier
% phase meassurements
% 
% Speedoflight = 2.99792458e8; %m/s
% 

%SensorNoiseParameters;

global GPS_PI OMEGAedot mu Earthradius Speedoflight c F L1_f L2_f gamma L1_Wavelength;
               
                                            
                     SatVector_ECEF = SatPos(1:3) - UserPos(1:3);
                     GeometricRange = norm(SatVector_ECEF);            %geometric range from Satellite to user
                                              
                                                   
                     
                     %add a earth rotation error as well to the PR..                                      
                     
                     delta_pr_omegaedot = -(OMEGAedot / Speedoflight) * (SatPos(1) * UserPos(2) - SatPos(2) * UserPos(1));
                                                  
                     %delta_pr_omegaedot = 0;                     
                     %Because you add the delta_pr_omegaedot to the PR's to correct for the earth rotation, you have to subtract it here if you want to add the earth 
                     %rotation error.                   
                     
                                                  
                     
                     
                     % sat clock here is actually in seconds - might need
                     % to multiply by Speedoflight to convert to metres...                     
                     
                     PRMeasPredict = GeometricRange + UserPos(4) - SatPos(4)  - delta_pr_omegaedot + GPSPr_noiseError;
                  
                                       
                     
         %PR rate calculations                         
        
         r_VecCalcVel = (SatVel(1) - UserVel(1))*(SatPos(1)-UserPos(1)) + (SatVel(2) - UserVel(2))*(SatPos(2)-UserPos(2)) + (SatVel(3) - UserVel(3))*(SatPos(3)-UserPos(3));
     
         Relative_Velocity = r_VecCalcVel/GeometricRange;
         
       

        PR_RateMeasPredict = Relative_Velocity + UserVel(4) - SatVel(4) + GPSPr_Rate_noiseError;
         
                
         
         
         
         
         
        
        %PR_RateMeasPredict = PRMeasPredict - PRMeasPredict_Previous  ; %form PR rates by the difference between current and previous PR's. 
       

       
