function [PRMeasPred, PRRMeasPredict, ICPMeasPredict, CPMeasPredict,NAmb, GeometricRange] = GARDSim_PRPredict(SatPos,UserPos, IonoDelay, time_step, PRMeas_previous, ICPMeas_previous,CPMeas_previous,GeometricRange_Previous,EstimatedClockBias);



global GPS_PI OMEGAedot mu Earthradius Speedoflight c F L1_f L2_f gamma L1_Wavelength;
% persistent PRMeasPred
% can use mlock 
% %store persistent variable in this function ~~??
% if isempty PRMeasPred 

%$Id: GARDSim_PRPredict.m 1880 2008-07-15 05:21:10Z n2523710 $

%$Id: GARDSim_PRPredict.m 1880 2008-07-15 05:21:10Z n2523710 $

%Troy Bruggemann 30 June 2005
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



global GPS_PI OMEGAedot mu Earthradius Speedoflight c F L1_f L2_f gamma L1_Wavelength;



      

                    %Increase these values if you want the position solution to be noisier (ie larger jumps between solutions)
                    % a = 1; b = 5;  %range of values for random number %note. b has to be the negative value one. 
                     
                    a = 1; b = 5;  %note, changing these does not give any bias on the position solution, just the random noise. 
                    MeasurementPrecision = a + (b-a) * rand(1);
                    %MeasurementPrecision = 3;  %metres. If you want this to be a fixed value. The amount of noise on the pseudorange measurement due to the measurement precision.
                                            
                     SatVector_ECEF = SatPos(1:3) - UserPos(1:3);
                     Slant_range = norm(SatVector_ECEF);            %geometric range from Satellite to user
                                              
                     GeometricRange = Slant_range;                                   
                     
                     %add a earth rotation error as well to the PR..                                      
                     
                     delta_pr_omegaedot = -(OMEGAedot / Speedoflight) * (SatPos(1) * UserPos(2) - SatPos(2) * UserPos(1));
                                                  
                                          
                     %Because you add the delta_pr_omegaedot to the PR's to correct for the earth rotation, you have to subtract it here if you want to add the earth 
                     %rotation error. 
                     
                     NoisePR = (GeometricRange - delta_pr_omegaedot - SatPos(4) + IonoDelay)/MeasurementPrecision;
                     
                     NoisePR = floor(NoisePR);
                     
                     Slant_range_with_noise = MeasurementPrecision*NoisePR + UserPos(4) ;                                        
                     
                     PRMeasPredict = Slant_range_with_noise;                          
                     
                     
                     %------------------Calculate PR Rate -----------------------
                     % Added DG 6 Oct 05 - calcualte pseudorange rate
                     PRRMeasPredict = (PRMeasPredict - PRMeas_previous) / time_step; 
                     
                     %------------------Calculate Carrier Phase-----------------------
                     
                     NAmb = 0;  %Integer Ambiguity, is time independant. 
                     
                   
                     %generate a cp measurement by finding the relative
                     %change in distance between current epoch and last
                     %epoch. 
                     
                     
                     %increase these values if noisier carrier phase
                     %measurements are required. It's based on that a
                     %receiver can measure carrier phase to 1% precision.
                     
                     %following values are fractions of an L1 cycle
                     %aCP = 0.045; bCP = 0.055;  %range of values for random number %note. b has to be the negative value one. 
                     
                     %aCP = 0.5; bCP = 0.4; 
                     %MeasurementPrecisionCPTemp = aCP + (bCP-aCP) * rand(1);
                     
                     MeasurementPrecisionCPTemp = 1.0;  %This seems to give the noise on the carrier phase to be looking like ht ashtech micro Z. 
        

                     MeasurementPrecisionCP = MeasurementPrecisionCPTemp*L1_Wavelength; %can measure carrier phase to about 1% precision of a wavelength

                     
                     
                     %Express carrier phase as a range
                     %CPMeasCurrent = (GeometricRange - delta_pr_omegaedot + UserPos(4) - SatPos(4) - IonoDelay) + NAmb;
                     
                     
                     %The UserPos(4) corrupts the measurement with the
                     %clock bias, and the EstimatedClockBias is the
                     %estimate of this clock bias. 
                     
%                      plop = SatPos(4)
%                      jump = UserPos(4)
%                      goo = EstimatedClockBias
%                      
%                      billy = UserPos(4)-EstimatedClockBias
                     
                     
                     %CPMeasCurrent = (GeometricRange + UserPos(4)-EstimatedClockBias - SatPos(4)) + NAmb;
                     
                     
                     %here it is assuming that the receiver clock bias can
                     %be estimated accuractly and therefore is not seen on
                     %the carrier phase measurement. If it is included in
                     %the carrier phase measurement then the carrier phase
                     %is very noisy. 
                     
                     
                     CPMeasCurrent = (GeometricRange - SatPos(4) - IonoDelay) + NAmb; %the satellite clock error is very small
                     
                     
                     
                     %CPMeasCurrent = GeometricRange;
                     %not sure if UserPos should be added separate or
                     %included in the above
                     
                     
                     NoiseCP = (CPMeasCurrent)/MeasurementPrecisionCP;
                     
                     NoiseCP = floor(NoiseCP);
                    
                     CPMeasCurrentWithNoise = MeasurementPrecisionCP*NoiseCP;   
                     
                     
                     CPMeasCurrentWithNoise = CPMeasCurrentWithNoise;% + UserPos(4); %add receiver clock bias
                     
                     
                     %now subtract from the previous to give the relative
                     %change in carrier phase between previous epoch and
                     %current epoch.
                     Rel_Dist_with_noise = CPMeasCurrentWithNoise - CPMeas_previous;
                                       
                                          
                     CPMeasPredict = CPMeasCurrentWithNoise;       %this is the carrier phase expressed as a range                       
                                          
                    
                     %Add Rel_Dist_with_noise to the previous ICP measurement to give
                     %integrated carrier phase
                     ICPMeasPredict = ICPMeas_previous + Rel_Dist_with_noise; %integrated carrier phase in L1 CYCLES
                     
                     
                    
                     
      
%For function output.                     
PRMeasPred = PRMeasPredict;




