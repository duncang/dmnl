

function gps = GARD_ReadNovatelLogData(inputfile)

%  gps = GARD_ReadNovatelLogData(inputfile)
%
%  inputfile is a comma delimited log file from ESM

%inputfile = 'log_20090722_gps.out';

fd = fopen(inputfile);

RangeRecord = 0;
EphemRecord = 0;
BestXYZRecord = 0;
IONUTCRecord = 0;

while ~feof(fd)
    
   newline = fgetl(fd);
   
   % check for version 
   if strcmp(newline(1:3),'37,')
        disp(sprintf('Logged VERSION string: %s', newline));
   else
       
       linedata = str2num(newline);

       switch linedata(1)
           case 7,
               % GPS EPHEM - size = 36
               if length(linedata) ~= 36
                   break;
               else
                   prn = linedata(5);
                   if (prn ~= 0)
                       disp(sprintf('Got EPHEM for PRN %d',prn)); 
                       EphemRecord = EphemRecord+1;
                       gps.GPSEphem(EphemRecord).prn = prn;
                       gps.GPSEphem(EphemRecord).rtTimestamp = linedata(2);
                       gps.GPSEphem(EphemRecord).GPSWeek = linedata(3);
                       gps.GPSEphem(EphemRecord).GPSSec = linedata(4);
                       gps.GPSEphem(EphemRecord).dTOW = linedata(6);
                       gps.GPSEphem(EphemRecord).ulHealth = linedata(7);
                       gps.GPSEphem(EphemRecord).ulIODE1 = linedata(8);
                       gps.GPSEphem(EphemRecord).ulIODE2 = linedata(9);
                       gps.GPSEphem(EphemRecord).ulGPSWeek = linedata(10);
                       gps.GPSEphem(EphemRecord).ulZWeek = linedata(11);
                       gps.GPSEphem(EphemRecord).dTOE = linedata(12); % /* Time of Ephemeris  */
                       gps.GPSEphem(EphemRecord).dA = linedata(13); % 		/* Semi-major axis, metres */
                       gps.GPSEphem(EphemRecord).dDeltaN = linedata(14); % 		/* Mean Motion Difference, radians/sec */
                       gps.GPSEphem(EphemRecord).dM0 = linedata(15); % 		/* mean anomoly, radians */
                       gps.GPSEphem(EphemRecord).dEccentricity = linedata(16); % 	/* eccentricity */
                       gps.GPSEphem(EphemRecord).dOmega = linedata(17); % 		/* argument of perigee */
                       gps.GPSEphem(EphemRecord).dcuc = linedata(18); % 
                       gps.GPSEphem(EphemRecord).dcus = linedata(19); % 
                       gps.GPSEphem(EphemRecord).dcrc = linedata(20); % 
                       gps.GPSEphem(EphemRecord).dcrs = linedata(21); % 
                       gps.GPSEphem(EphemRecord).dcic = linedata(22); % 
                       gps.GPSEphem(EphemRecord).dcis = linedata(23); % 
                       gps.GPSEphem(EphemRecord).dInclination0 = linedata(24); % /* inclination  */
                       gps.GPSEphem(EphemRecord).dInclination_dot = linedata(25); % 	/* inclination rate */
                       gps.GPSEphem(EphemRecord).dOmega0 = linedata(26); % 		/* right ascention */	
                       gps.GPSEphem(EphemRecord).dOmega_dot = linedata(27); % 
                       gps.GPSEphem(EphemRecord).ulIODC = linedata(28); % 
                       gps.GPSEphem(EphemRecord).dTOC = linedata(29); % 		/* SV clock correction, seconds */
                       gps.GPSEphem(EphemRecord).dTGD = linedata(30); % 		/* estimated group delay */
                       gps.GPSEphem(EphemRecord).dA_f0 = linedata(31); %	/* clock aging parameter, seconds */
                       gps.GPSEphem(EphemRecord).dA_f1 = linedata(32); % 		/* clock aging parameter, seconds/second */
                       gps.GPSEphem(EphemRecord).dA_f2 = linedata(33); % 		/* clock aging parameter, seconds/second/second */
                       gps.GPSEphem(EphemRecord).ulAntiSpoofing = linedata(34); % 	/* anti-spoofing flag 0=FALSE, 1=TRUE */
                       gps.GPSEphem(EphemRecord).dN = linedata(35); % 		/* corrected mean motion, radians/sec */
                       gps.GPSEphem(EphemRecord).dURA = linedata(36); % 		/* User-Range Accuracy	 */    

                   end


               end
           case 43,
               % range - length = 244
               if length(linedata == linedata(5) * 10 + 5)
                  RangeRecord = RangeRecord+1;
                  gps.RangeData(RangeRecord).rtTimeStamp = linedata(2);
                  gps.RangeData(RangeRecord).GPSWeek = linedata(3);
                  gps.RangeData(RangeRecord).GPSSec = linedata(4);
                  gps.RangeData(RangeRecord).lNumberObservations = linedata(5);
                  
                  for i = 1:gps.RangeData(RangeRecord).lNumberObservations
                     gps.RangeData(RangeRecord).Obs(i).usPRN = linedata(5+((i-1)*10)+1);
                     gps.RangeData(RangeRecord).Obs(i).usGlonassFrequency = linedata(5+((i-1)*10)+2);
                     gps.RangeData(RangeRecord).Obs(i).dPseudorange = linedata(5+((i-1)*10)+3);
                     gps.RangeData(RangeRecord).Obs(i).fPseudorangeSigma = linedata(5+((i-1)*10)+4);
                     gps.RangeData(RangeRecord).Obs(i).dCarrierPhase = linedata(5+((i-1)*10)+5);
                     gps.RangeData(RangeRecord).Obs(i).fCarrierPhaseSigma = linedata(5+((i-1)*10)+6);
                     gps.RangeData(RangeRecord).Obs(i).fDoppler = linedata(5+((i-1)*10)+7);
                     gps.RangeData(RangeRecord).Obs(i).fCNo = linedata(5+((i-1)*10)+8);
                     gps.RangeData(RangeRecord).Obs(i).fLockTime = linedata(5+((i-1)*10)+9);
                     gps.RangeData(RangeRecord).Obs(i).ulTrackingStatus = linedata(5+((i-1)*10)+10);
                  end
                  
               end
               
           case 241,
               if length(linedata == 31)
                    BestXYZRecord = BestXYZRecord + 1;
                   gps.BestXYZData(BestXYZRecord).rtTimeStamp = linedata(2);
                   gps.BestXYZData(BestXYZRecord).GPSWeek = linedata(3);
                   gps.BestXYZData(BestXYZRecord).GPSSec = linedata(4);
                   gps.BestXYZData(BestXYZRecord).lPosSolutionStatus = linedata(5);
                   gps.BestXYZData(BestXYZRecord).lPosiitionType = linedata(6);
                   gps.BestXYZData(BestXYZRecord).dPosX = linedata(7);
                   gps.BestXYZData(BestXYZRecord).dPosY = linedata(8);
                   gps.BestXYZData(BestXYZRecord).dPosZ = linedata(9);
                   gps.BestXYZData(BestXYZRecord).fPosXSigma = linedata(10);
                   gps.BestXYZData(BestXYZRecord).fPosYSigma = linedata(11);
                   gps.BestXYZData(BestXYZRecord).fPosZSigma = linedata(12);
                   gps.BestXYZData(BestXYZRecord).lVelSolutionStatus =linedata(13);
                   gps.BestXYZData(BestXYZRecord).lVelocityType =linedata(14);
                   gps.BestXYZData(BestXYZRecord).dVelX = linedata(15);
                   gps.BestXYZData(BestXYZRecord).dVelY = linedata(16);
                   gps.BestXYZData(BestXYZRecord).dVelZ = linedata(17);
                   gps.BestXYZData(BestXYZRecord).fVelXSigma = linedata(18);
                   gps.BestXYZData(BestXYZRecord).fVelYSigma = linedata(19);
                   gps.BestXYZData(BestXYZRecord).fVelZSigma =linedata(20);
                   gps.BestXYZData(BestXYZRecord).fVelocityLatency = linedata(25);
                   gps.BestXYZData(BestXYZRecord).fDifferentialAge =linedata(26);
                   gps.BestXYZData(BestXYZRecord).fSolutionAge =linedata(27);
                   gps.BestXYZData(BestXYZRecord).ucNumberObservationsTracked =linedata(28); 
                   gps.BestXYZData(BestXYZRecord).ucNumberL1ObservationsUsed =linedata(29);
                   gps.BestXYZData(BestXYZRecord).ucNumberL1ObservationsAboveRTKMaskAngle =linedata(30);
                   gps.BestXYZData(BestXYZRecord).ucNumberL2ObservationsAboveRTKMaskAngle =linedata(31);
                   
               end

           case 8,
               if (length(linedata == 21))
                   IONUTCRecord = IONUTCRecord+1;
                   gps.IONUTCData(IONUTCRecord).rtTimeStamp = linedata(2);
                   gps.IONUTCData(IONUTCRecord).GPSWeek = linedata(3);
                   gps.IONUTCData(IONUTCRecord).GPSSec = linedata(4);
                   gps.IONUTCData(IONUTCRecord).a0 = linedata(5);
                   gps.IONUTCData(IONUTCRecord).a1 = linedata(6);
                   gps.IONUTCData(IONUTCRecord).a2 = linedata(7);
                   gps.IONUTCData(IONUTCRecord).a3 =linedata(8);
                   gps.IONUTCData(IONUTCRecord).b0 = linedata(9);
                   gps.IONUTCData(IONUTCRecord).b1 = linedata(10);
                   gps.IONUTCData(IONUTCRecord).b2 =linedata(11);
                   gps.IONUTCData(IONUTCRecord).b3 = linedata(12);
                   gps.IONUTCData(IONUTCRecord).ulUTCWeekNumber =linedata(13); 
                   gps.IONUTCData(IONUTCRecord).ulUTCReferenceTime = linedata(14);
                   gps.IONUTCData(IONUTCRecord).A0 = linedata(15);
                   gps.IONUTCData(IONUTCRecord).A1 = linedata(16);
                   gps.IONUTCData(IONUTCRecord).ulUTCFutureWeekNumber = linedata(17);
                   gps.IONUTCData(IONUTCRecord).ulUTCFutureDayNumber = linedata(18);
                   gps.IONUTCData(IONUTCRecord).lUTCLeapSeconds = linedata(19);
                   gps.IONUTCData(IONUTCRecord).lUTCFutureLeapSeconds = linedata(20);
                   gps.IONUTCData(IONUTCRecord).ulUTCDeltaT = linedata(21);
               end
             

           otherwise,
               disp(sprintf('Unknown message %d', linedata(1)));

       end % switch
   
   end % strcmp
   
end % while

disp(sprintf('Got %d Range Obs',RangeRecord));
disp(sprintf('Got %d BestXYZ Obs', BestXYZRecord));

fclose(fd);





