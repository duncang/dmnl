function [SV_X_Data, SV_Y_Data, SV_Z_Data, SV_T_Data, ValidData, URA] = GPSOrbitPropagator(SVWeek, SVTime, SV, NavigationData, EphValidTime)
% [SV_X_Data, SV_Y_Data, SV_Z_Data, SV_T_Data, ValidData] =
% GPSOrbitPropagator(SVWeek, SVTime, SV, NavigationData, EphValidTime)
%
%  This function calculates the predicted X, Y, Z and T for the specified
%  satelltie based on the broadcast ephemeris data
%
%  SVTime is the epoch at which the orbit position is required
%  SV is the desired satellite to use.
%  NavigationData is the navigation data read by 'freadnav'.
%  EphValidTime is the time that the ephemeris data is valid in seconds.
%  if not supplied 7500 seconds is used (2 hours 5 mins).
%
%  Written by Duncan Greer and Troy Bruggemann (c) QUT GRAS Project 2005
%  Last Modified:  $Id: GPSOrbitPropagator.m 4294 2011-01-30 05:50:56Z greerd $
%
% constants used by this function
% pi = 3.14159126535898; % GPS value for PI
% OMEGAedot = 7.2921151467e-5; % Earth rotation rate in radians per second
% mu = 3.986005e14; % WGS-84 valeu of the earths universal gravitational parameter in m^3/s^2
% c = 2.99792458e8;  % speed of light m/s
% F = -4.442807633e-10; % a random number from the ICD page 88
% 
% navigation data returned by freadnav:
    %     1: PRN_nav_vec(i) 
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
    %     32: toc_nav_week 

    
GPSConstants;
global GPS_PI OMEGAedot mu Earthradius Speedoflight c F L1_f L2_f gamma L1_Wavelength;    
  
    
    % find the indexes of the satellite in the navigation data which has
    % ephemeris within +/- 2 hours of the data epoch we wantt o calculate for
    
    %Find where the occurances of the satellite are in the file
    SVIndexArray = find(NavigationData(:,1) == SV);
    % check if SV does not exist in this nav file
    
    if(isempty(SVIndexArray))
        %disp(sprintf('No SV data for SV %d',SV));
        
        SV_X_Data = 0;
        SV_Y_Data = 0;
        SV_Z_Data = 0;
        SV_T_Data = 0;
        ValidData = 0;

        return;
    end


    %Check the times

    %SVIndex is an array with all the indices of where those satellites are

    %check for size of SVIndex

    sizeSVIndexArray = size(SVIndexArray);


    for p = 1:sizeSVIndexArray(1)

        SVIndexTemp = SVIndexArray(p);   %get the index number

        %toe for the satellite is:
        toe = NavigationData(SVIndexTemp,14);
        toe_week = NavigationData(SVIndexTemp,24);

        %check that ephemeris falls within +/- 2 hours 5 minutes

        Week = SVWeek - toe_week;
        Secs = SVTime - toe;
        %number of seconds away

        TimeAway(p) = abs(Week*604800 + Secs);
     
    end
    
    %find the one which is closest to the epoch we want to calculate the
    %satellite position for
    
    closest = min(TimeAway);

    if(~exist('EphValidTime'))
        EphValidTime = 7500;
    end
    
    if (closest < EphValidTime) %7500 is 2 hours and 5 mins

        Index = find(TimeAway == closest);

        SVIndex = SVIndexArray(Index(1));  %use Index(1) because its the first because Index might be an array if theres more than one closest time.

        %check the SV health flag, 0 means its okay
        if (NavigationData(SVIndex,27) == 0)

            DataValid = 1;  %there is data in the right time epoch and is a healthy satellite

        else
            disp(sprintf('SV%d: Unhealthy',SV));
            DataValid = 0;
        end

    else
        %disp(sprintf('SV%d: Closest ephem is %d sec old',SV,closest));
        DataValid = 0;
    end
  
    
    
    
    if (DataValid == 1)

        %     % get the satellite approximate time of transmission -
        %     ts = SVRXTime - PR(SV,2) / c;
        %ts = SVTime - PR/c;
        ts = SVTime;
        % make initial estimate of Ek - this is revised later with the
        % corrected tk

        toe = NavigationData(SVIndex,14);
        A = NavigationData(SVIndex,13)^2;
        n0 = sqrt(mu/A^3);
        tk = ts - toe;

        % check for end of week roll-over
        if tk > 302400
            tk = tk - 604800;
        elseif tk < -302400
            tk = tk + 604800;
        end

        

        n = n0 + NavigationData(SVIndex,8);
        Mk = NavigationData(SVIndex,9) + n*tk;
        e = NavigationData(SVIndex,11);
        Ek = KepplerSolver(Mk,e);


        % get the clock model coefficients
        af0 = NavigationData(SVIndex,3);
        af1 = NavigationData(SVIndex,4);
        af2 = NavigationData(SVIndex,5);

        % calculate the relativistic correction
        delta_tr = F * e * sqrt(A) * sin(Ek);

        toc = NavigationData(SVIndex,2);

        % calculate the satellite time offset
        delta_ts = af0 + af1*(ts - toc) + af2*(ts - toc)^2 + delta_tr;

        % correct delta_ts for iono group delay
        TGD = NavigationData(SVIndex,28);
        delta_ts = delta_ts - TGD;

        
        %-----added by TB 21.9.05 -----------
        ts = SVTime - delta_ts;
       
        % make initial estimate of Ek - this is revised later with the
        % corrected tk

        toe = NavigationData(SVIndex,14);
        A = NavigationData(SVIndex,13)^2;
        n0 = sqrt(mu/A^3);
        tk = ts - toe;

        % check for end of week roll-over
        if tk > 302400
            tk = tk - 604800;
        elseif tk < -302400
            tk = tk + 604800;
        end


        n = n0 + NavigationData(SVIndex,8);
        Mk = NavigationData(SVIndex,9) + n*tk;
        e = NavigationData(SVIndex,11);
        Ek = KepplerSolver(Mk,e);


        % get the clock model coefficients
        af0 = NavigationData(SVIndex,3);
        af1 = NavigationData(SVIndex,4);
        af2 = NavigationData(SVIndex,5);

        % calculate the relativistic correction
        delta_tr = F * e * sqrt(A) * sin(Ek);

        toc = NavigationData(SVIndex,2);

        % calculate the satellite time offset
        delta_ts = af0 + af1*(ts - toc) + af2*(ts - toc)^2 + delta_tr;

        % correct delta_ts for iono group delay
        TGD = NavigationData(SVIndex,28);
        delta_ts = delta_ts - TGD;          
                     
        URA = sqrt(NavigationData(SVIndex,26));
        %------------------------------------
        
        

        % calculate the GPS system time
        GPSTime = ts - delta_ts;                   

        % calculate the orbital coeficients or extract from navigation
        % message

        % revise estimate of Ek based on updated tk
        tk = GPSTime - toe;
        
        %-----added by TB 21.9.05 -----------
        % check for end of week roll-over
        if tk > 302400
            tk = tk - 604800;
        elseif tk < -302400
            tk = tk + 604800;
        end
        %-------------------------------------
        
        
        n = n0 + NavigationData(SVIndex,8);
        Mk = NavigationData(SVIndex,9) + n*tk;
        e = NavigationData(SVIndex,11);
        Ek = KepplerSolver(Mk,e);

        vk = atan2((sqrt(1 - e^2)*sin(Ek)/(1-e*cos(Ek))),((cos(Ek) - e)/(1 - e*cos(Ek))));
        PHIk = vk + NavigationData(SVIndex,20);

        % orbital pertubations
        Cus = NavigationData(SVIndex,12);
        Cuc = NavigationData(SVIndex,10);
        Crs = NavigationData(SVIndex,7);
        Crc = NavigationData(SVIndex,19);
        Cis = NavigationData(SVIndex,17);
        Cic = NavigationData(SVIndex,15);


        delta_uk = Cus * sin(2*PHIk) + Cuc * cos(2*PHIk);
        delta_rk = Crs * sin(2*PHIk) + Crc * cos(2*PHIk);
        delta_ik = Cis * sin(2*PHIk) + Cic * cos(2*PHIk);

        IDOT = NavigationData(SVIndex,22);
        i0 = NavigationData(SVIndex,18);

        uk = PHIk + delta_uk;
        rk = A*(1 - e*cos(Ek)) + delta_rk;
        ik = i0 + delta_ik + IDOT * tk;

        xk_dash = rk * cos(uk);
        yk_dash = rk * sin(uk);

        OMEGA0 = NavigationData(SVIndex,16);
        OMEGAdot = NavigationData(SVIndex,21);

        OMEGAk = OMEGA0 + (OMEGAdot - OMEGAedot) * tk - OMEGAedot*toe;

        % final satellite positions in ecef
        xk = xk_dash * cos(OMEGAk) - yk_dash * cos(ik)*sin(OMEGAk);
        yk = xk_dash * sin(OMEGAk) + yk_dash * cos(ik)*cos(OMEGAk);
        zk = yk_dash * sin(ik);

        %SVPosition(SV,:) = [xk  yk zk];

        SV_X_Data = xk;
        SV_Y_Data = yk;
        SV_Z_Data = zk;
        SV_T_Data = delta_ts;
        ValidData = 1;
       
    else  % DataValid = 0, no ephemeris was found for the satellite within the right time
        SV_X_Data = 0;
        SV_Y_Data = 0;
        SV_Z_Data = 0;
        SV_T_Data = 0;
        URA = 0;
        ValidData = 0;
       % return;

    end
    
       
    %verified ephemeris and satellite selection from navigation data :
    %6/9/05 TB
    %successfully (1)ignores satellite ephemeris which has bad health flag
    %(2) doesn't use ephemeris more than 2 hours 5 min old
    %(3) ignores satellites which aren't in the ephemeris
    

    
    