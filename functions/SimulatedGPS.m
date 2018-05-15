function [GPSPosLSQ, SatPos,SatVel,SV_AboveElevationMask,PRMeasured_Simulated,PRRateMeasured_Simulated,N, Sat_PRN_Vec] = SimulatedGPS(PosTruth, VelTruth, GPSWeek, GPSSec,ElevationMask,NavData,PHI,THETA,PSI);

%$Id: SimulatedGPS.m 1883 2008-07-15 05:53:55Z n2523710 $
%Troy Bruggemann  11 August 2005
GPSConstants; %All generic constants put in this script, as global variables;
%turn JIT accelerator on (speeds up code processing)


feature accel on %only works in matlab 7


%VelTruth is in ECEF x, y, z



%this is the GPS processor which calculates simulated pseudoranges and
%satellite positions for the  GPSandModelEKF.m


%------------------------------------------------------------------------
%Atmospheric Models
%-----------------------------------------------------------------------

%IONOSPHERIC PARAMETERS for IONO Model
ALPHA = [0.1118e-07  -0.7451e-08  -0.5960e-07   .01192e-06];%          ION ALPHA
BETA = [0.1167e+06  -0.2294e+06  -0.1311e+06   .01049e+07]; %         ION BETA


    %--------------------------
    %Optimised GPS Constellation - DO229C WAAS MOPS
    %---------------------------
    
        n = 1;
        
         Sat_PRN_Vec = zeros(1,32);  %initialise the vector with zeros. Columns of this vector is the PRN number of the satellite
         
        for SV = 1:32            
            
            [SV_X_Data(SV) SV_Y_Data(SV) SV_Z_Data(SV) SV_T_Data(SV) ValidData_Satellite(SV)] = GPSOrbitPropagatorOptimal(GPSWeek, GPSSec, SV, NavData);

            [SV_Xvel_Data(SV) SV_Yvel_Data(SV) SV_Zvel_Data(SV) SV_Tvel_Data(SV) SV_Xacc_Data(SV) SV_Yacc_Data(SV) SV_Zacc_Data(SV) SV_Tacc_Data(SV) ValidDataSatVels(SV)] = GPSOrbitPropagatorOptimalVelocities(GPSWeek, GPSSec, SV, NavData);
            
            
            
%             if ValidData_Satellite(SV) == 1 && ValidDataSatVels(SV) == 1
%             
%               
%             
%             
%             end
        
            
            
            
                %calculate satellite elevations 
                   
                 [Azimuth(SV), Elevation(SV)] = AzEl(PosTruth(1:3), [SV_X_Data(SV), SV_Y_Data(SV), SV_Z_Data(SV)]);   
                 Azimuth(SV) = rad2deg(Azimuth(SV));
                 Elevation(SV) = rad2deg(Elevation(SV));
                 
                 
                 [AzimuthBody(SV), ElevationBody(SV)] = AzElAntenna(PosTruth(1:3), [SV_X_Data(SV), SV_Y_Data(SV), SV_Z_Data(SV)],PHI,THETA,PSI);
                 
                 AzimuthBody(SV) = rad2deg(AzimuthBody(SV));
                 ElevationBody(SV)  = rad2deg(ElevationBody(SV));
                 
               
                 
                 
                 if (Elevation(SV) > ElevationMask)  && (ValidData_Satellite(SV) == 1 && ValidDataSatVels(SV) == 1)      
                     
                    
                    
                   SV_AboveElevationMask(n) = SV ;  %PRN numbers of the satellites above the elevation mask             
                    
                    
                    SatPos(SV,1) = SV_X_Data(SV);
                    SatPos(SV,2) = SV_Y_Data(SV);
                   SatPos(SV,3) = SV_Z_Data(SV);
                    SatPos(SV,4) = c*SV_T_Data(SV); %convert to metres
                     
                     
                    SatVel(SV,1) = SV_Xvel_Data(SV);
                    SatVel(SV,2) = SV_Yvel_Data(SV);
                   SatVel(SV,3) = SV_Zvel_Data(SV);
                   SatVel(SV,4) = c*SV_Tvel_Data(SV); %convert to metres
                     
                     
                    Sat_PRN_Vec(SV) = 1; 
                    
                    
                    
                    
                    %set ionodelay to 0 for now
              IonoDelay = 0;  
              [PRMeasured_Simulated(SV),PRRateMeasured_Simulated(SV),GeometricRange(SV)] = GARDSim_PRPredictSimulated([SV_X_Data(SV) SV_Y_Data(SV) SV_Z_Data(SV) SV_T_Data(SV)] ,[SV_Xvel_Data(SV) SV_Yvel_Data(SV) SV_Zvel_Data(SV) SV_Tvel_Data(SV)] ,PosTruth,VelTruth,IonoDelay);  
            
                    
                                   
                    
                    
                   
                    
                    n = n+1;        
                    
                    
                    
                        
                    
     
                 end                
                 
           end
        
        
         %number of satellites above elevation mask
         N = length(SV_AboveElevationMask);
    
        %N = 9; %always keep it at 9 satellites for now
        
        
%         %generate PR's for each satellite
%         
%         for SV = 1:32  
%         
%             
%             if Sat_PRN_Vec(SV) == 1
%                 
%                          IonoDelay = 0;  
%                    [PRMeasured_SimulatedTemp(SV),GeometricRangeTemp(SV)] = GARDSim_PRPredictSimulated(SatPos(SV,:),PosTruth,IonoDelay);    
%                 
%             end
%             
%             
%         end
%         
        
        
        
   
    
%pick out the satellites we want   for the Least squares



k = 1;
for SV = 1:32
        
    if Sat_PRN_Vec(SV) == 1

        PRMeasured_SimulatedLSQ(k) =  PRMeasured_Simulated(SV);

        
         SatPosLSQ(k,1) = SV_X_Data(SV);
                    SatPosLSQ(k,2) = SV_Y_Data(SV);
                   SatPosLSQ(k,3) = SV_Z_Data(SV);
                    SatPosLSQ(k,4) = c*SV_T_Data(SV); %convert to metres                                                             

        k = k+1;
    end
end     
    
      
     
         [GPSPosLSQ VarSolutionVec_Observed NumIterations_Observed_Vel ResVec_Observed M_Observed LSQ_Fail_Observed limit_Observed DOP_Observed(1:5)] = GARD_LSQ(PosTruth,N,PRMeasured_SimulatedLSQ,SatPosLSQ);   

        
        %UserPos = ([PosTruth(1,i),PosTruth(1,i),PosTruth(1,i),c*dTposTCXO(i)]); %this will be a copy of the values set by the user in the flight planning section
        
                
  



