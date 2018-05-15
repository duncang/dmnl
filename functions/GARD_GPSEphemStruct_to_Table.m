function ephem = GARD_GPSEphemStruct_to_Table(GPSEphem)

% ephem = GARD_GPSEphemStruct_to_Table(GPSEphem)
% converts a GPSEphem struct from the Novatel to an ephemerides table
% which is used directly by GPSOrbitPropagator
%
% $Id$


% GPSEphem is a 1xN struct array with fields 
%     rtTimestamp
%     GPSWeek
%     GPSSec
%     dTOW
%     ulHealth
%     ulIODE1
%     ulIODE2
%     ulGPSWeek
%     ulZWeek
%     dTOE
%     dA
%     dDeltaN
%     dM0
%     dEccentricity
%     dOmega
%     dcuc
%     dcus
%     dcrc
%     dcrs
%     dcic
%     dcis
%     dInclination0
%     dInclination_dot
%     dOmega0
%     dOmega_dot
%     ulIODC
%     dTOC
%     dTGD
%     dA_f0
%     dA_f1
%     dA_f2
%     ulAntiSpoofing
%     dN
%     dURA


%   Ephemerides table is a 32 x M array where M is the number of
%   ephemerides stored (may be duplicates)
%     1: PRN_nav_vec
%     2: toc_nav_sec 
%     3: svclkbias_nav 
%     4: svclkdrift_nav 
%     5: svclkdriftrate_nav 
%     6: IODE_nav
%     7: crs_nav 
%     8: deltan_nav 
%     9: M0_nav 
%     10: cuc_nav 
%     11: e_nav 
%     12: cus_nav 
%     13: sqrta_nav 
%     14: toe_nav 
%     15: cic_nav 
%     16: OMEGA0_nav
%     17: cis_nav 
%     18: i0_nav 
%     19: crc_nav 
%     20: omega_nav 
%     21: OMEGA_DOT_nav 
%     22: idot_nav 
%     23: CodeL2_nav 
%     24: GPSwknum_nav 
%     25: L2Pdataflag_nav 
%     26: SVaccuracy_nav 
%     27: SVhealth_nav 
%     28: tgd_nav 
%     29: IODC_nav
%     30: txtimeofmessage_nav 
%     31: fitinterval_nav
%     32:  toc_nav_week


EphemIndex = 0;
for EphemIndex = 1:length(GPSEphem)  %% max PRN used
    
        ephem(EphemIndex,1) = GPSEphem(EphemIndex).prn;
        ephem(EphemIndex,2) = GPSEphem(EphemIndex).dTOC;
        ephem(EphemIndex,3) = GPSEphem(EphemIndex).dA_f0;
        ephem(EphemIndex,4) = GPSEphem(EphemIndex).dA_f1;
        ephem(EphemIndex,5) = GPSEphem(EphemIndex).dA_f2;
        ephem(EphemIndex,6) = GPSEphem(EphemIndex).ulIODE1;
        ephem(EphemIndex,7) = GPSEphem(EphemIndex).dcrs;
        ephem(EphemIndex,8) = GPSEphem(EphemIndex).dDeltaN;
        ephem(EphemIndex,9) = GPSEphem(EphemIndex).dM0;
        ephem(EphemIndex,10) = GPSEphem(EphemIndex).dcuc;
        ephem(EphemIndex,11) = GPSEphem(EphemIndex).dEccentricity;
        ephem(EphemIndex,12) = GPSEphem(EphemIndex).dcus;
        ephem(EphemIndex,13) = GPSEphem(EphemIndex).dA^0.5;
        ephem(EphemIndex,14) = GPSEphem(EphemIndex).dTOE;
        ephem(EphemIndex,15) = GPSEphem(EphemIndex).dcic;
        ephem(EphemIndex,16) = GPSEphem(EphemIndex).dOmega0;
        ephem(EphemIndex,17) = GPSEphem(EphemIndex).dcis;
        ephem(EphemIndex,18) = GPSEphem(EphemIndex).dInclination0;
        ephem(EphemIndex,19) = GPSEphem(EphemIndex).dcrc;
        ephem(EphemIndex,20) = GPSEphem(EphemIndex).dOmega;
        ephem(EphemIndex,21) = GPSEphem(EphemIndex).dOmega_dot;
        ephem(EphemIndex,22) = GPSEphem(EphemIndex).dInclination_dot;
        ephem(EphemIndex,23) = 0;
        ephem(EphemIndex,24) = GPSEphem(EphemIndex).ulGPSWeek;
        ephem(EphemIndex,25) = 0;
        ephem(EphemIndex,26) = GPSEphem(EphemIndex).dURA;
        ephem(EphemIndex,27) = GPSEphem(EphemIndex).ulHealth;
        ephem(EphemIndex,28) = GPSEphem(EphemIndex).dTGD;
        ephem(EphemIndex,29) = GPSEphem(EphemIndex).ulIODC;
        ephem(EphemIndex,30) = GPSEphem(EphemIndex).dTOW;
        ephem(EphemIndex,31) = 0;
        ephem(EphemIndex,32) = GPSEphem(EphemIndex).ulZWeek;
    
end