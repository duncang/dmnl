function [SV_X_Data, SV_Y_Data, SV_Z_Data, SV_T_Data, SV_Xvel_Data, SV_Yvel_Data, SV_Zvel_Data, SV_Tvel_Data, SV_Xacc_Data, SV_Yacc_Data, SV_Zacc_Data, SV_Tacc_Data, ValidData] = GPSOrbitPropagatorOptimalVelocitiesFixed(SVWeek,SVTime, SV, NavigationData)
%
%
%  This function calculates the predicted X, Y, Z and T velocities and accelerations for the specified
%  satellite based on the broadcast ephemeris data
%
%  SVTime is the epoch at which the orbit position is required
%  SV is the desired satellite to use.
%  NavigationData is the navigation data read by 'freadnav'.
%
%  Written by Troy Bruggemann (c) QUT GRAS Project 2005
%  Last Modified:  1 September 2005
%
% constants used by this function
% pi = 3.14159265358979; % GPS value for PI
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
%
% $Id: GPSOrbitPropagatorOptimalVelocities.m 1883 2008-07-15 05:53:55Z n2523710 $
%
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
        
        
        SV_Xvel_Data = 0;
        SV_Yvel_Data = 0;
        SV_Zvel_Data = 0;
        SV_Tvel_Data = 0;

        SV_Xacc_Data = 0;
        SV_Yacc_Data = 0;
        SV_Zacc_Data = 0;
        SV_Tacc_Data = 0;
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

    if (closest < 75000000000000000000000000000000000000000000000000000000) %don't want this check for the optimal constellation, applies to broadcast ephemeris only%7500 is 2 hours and 5 mins


        Index = find(TimeAway == closest);

        SVIndex = SVIndexArray(Index(1));  %use Index(1) because its the first because Index might be an array if theres more than one closest time.

        %check the SV health flag, 0 means its okay
        if (NavigationData(SVIndex,27) == 0)

            DataValid = 1;  %there is data in the right time epoch and is a healthy satellite

        else
            DataValid = 0;
        end

    else
        DataValid = 0;
    end


    if (DataValid == 1)
        %     % get the satellite approximate time of transmission -
        %     ts = SVRXTime - PR(SV,2) / c;
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


        %Mean anomaly (rads) and its derivative (rads/s) at tk.
        
        
            Mk = NavigationData(SVIndex,9) + n*tk;
        
        Mkdot = sqrt(mu/(A*A*A)) + NavigationData(SVIndex,8); %derivative of Mk in rad/s
        
       % Mk = NavigationData(SVIndex,9) + Mkdot*tk;
        
        
        
        % Mk = NavigationData(SVIndex,9) + Mkdot*tk;  %Is possible to use value
        % of Mkdot here, improved accuracy?


        e = NavigationData(SVIndex,11);
        Ek = KepplerSolver(Mk,e);
        Ekdot = Mkdot/(1.0-e*cos(Ek));

        Ekdotdot = -(Mkdot*e*sin(Ek)*Ekdot)/(1-e*cos(Ek))^2;


        % get the clock model coefficients
        af0 = NavigationData(SVIndex,3);
        af1 = NavigationData(SVIndex,4);
        af2 = NavigationData(SVIndex,5);

        % calculate the relativistic correction, will be zero because e is
        % zero for optimal constellation
        delta_tr = F * e * sqrt(A) * sin(Ek);

        delta_tr_dot = F * e * sqrt(A) * cos(Ek)*Ekdot; %I don't know if this is correct or not but think it is TB 1/9/05
        delta_tr_dotdot = F * e * sqrt(A) * cos(Ek)*Ekdotdot - Ekdot*F*e*sqrt(A)*sin(Ek)*Ekdot;                 %I don't know if this is correct or not but seems to be, in the same order of magnitude if you just
        %subtract the velocities TB 1/9/05..about 3e-16 seconds which is negligable

        toc = NavigationData(SVIndex,2);

        % calculate the satellite time offset
        delta_ts = af0 + af1*(ts - toc) + af2*(ts - toc)^2 + delta_tr;

        % calculate the satellite time drift
        delta_tsdot = af1 + 2*af2*(ts - toc) + delta_tr_dot;

        % calculate the satellite time drift rate (acceleration)
        delta_tsdotdot =  2*af2 + delta_tr_dotdot;

        % correct delta_ts for iono group delay
        TGD = NavigationData(SVIndex,28);
        delta_ts = delta_ts - TGD;
        
        
             
           
        
        
        
        

        % store the delta_ts to apply to PRs
       % SVBias(SV) = c * delta_ts;

        % calculate the GPS system time
        GPSTime = ts - delta_ts;

        % calculate the orbital coeficients or extract from navigation
        % message

        % revise estimate of Ek based on updated tk
        tk = GPSTime - toe;
        
        
        
                
%         %-----added by TB 21.9.05 -----------
%         % check for end of week roll-over - Doesn't apply for optimal
%         constellation - broadcast ephemeris only
%         if tk > 302400
%             tk = tk - 604800;
%         elseif tk < -302400
%             tk = tk + 604800;
%         end
%         %-------------------------------------
        
        
        
        %=================================================
        %added NOVEMBER 2007 - this 3 lines is the bits that were in
        %gpsorbitpropoptimal but werent in this file
        Week = SVWeek - toe_week;
         Secs = SVTime - toe;
        
         tk = Week*604800 + Secs;  %total number of seconds between reference epoch and epoch we want to calculate for
         %=================================================
        
        
             
        
        
        
        
        
        
        n = n0 + NavigationData(SVIndex,8);
        Mk = NavigationData(SVIndex,9) + n*tk;
        Mkdot = sqrt(mu/(A*A*A)) + NavigationData(SVIndex,8); %derivative of Mk in rad/s
		%Mk = NavigationData(SVIndex,9) + Mkdot*tk; 

        e = NavigationData(SVIndex,11);
        Ek = KepplerSolver(Mk,e);
        Ekdot = Mkdot/(1.0-e*cos(Ek));

        %vk = atan2((sqrt(1 - e^2)*sin(Ek)/(1-e*cos(Ek))),((cos(Ek) - e)/(1 - e*cos(Ek))));

        vk = atan2((sqrt(1 - e^2)*sin(Ek)),(cos(Ek) - e)); %note this is not the same as pk (Phik below)

     
        PHIk = vk + NavigationData(SVIndex,20);  
        
        
        PHIkdot = sqrt(1 - e^2)*Ekdot/(1.0-e*cos(Ek));
        
        

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
        ukdot = PHIkdot*(1.0 + 2.0*(Cus*cos(2*PHIk) - Cuc*sin(2*PHIk)));


        rk = A*(1 - e*cos(Ek)) + delta_rk;
        rkdot = A*e*sin(Ek)*Ekdot + 2.0*PHIkdot*(Crs*cos(2*PHIk) - Crc*sin(2*PHIk));


        ik = i0 + delta_ik + IDOT * tk;

        ikdot = IDOT + 2.0*PHIkdot*(Cis*cos(2*PHIk) - Cic*sin(2*PHIk));


        % Compute the satellite's position vector in its orbital plane
        %   and its derivative.
        xk_dash = rk * cos(uk);
        yk_dash = rk * sin(uk);


        xkdot_dash = rkdot*cos(uk) - yk_dash*ukdot;
        ykdot_dash = rkdot*sin(uk) + xk_dash*ukdot;


        % Compute the longitude of the ascending node (rads) and its
        %   derivative (rads/s).

        OMEGA0 = NavigationData(SVIndex,16);
        OMEGAdot = NavigationData(SVIndex,21);
        
        
        
        OMEGAkdot = OMEGAdot - OMEGAedot;

        OMEGAk = OMEGA0 + OMEGAkdot* tk - OMEGAedot*toe;


      

        % final satellite positions in ecef
        xk = xk_dash * cos(OMEGAk) - yk_dash * cos(ik)*sin(OMEGAk);
        yk = xk_dash * sin(OMEGAk) + yk_dash * cos(ik)*cos(OMEGAk);
        zk = yk_dash * sin(ik);

        % final satellite velocities in ecef

        tmp = ykdot_dash*cos(ik) - yk_dash*sin(ik)*ikdot;

        xkvel = -OMEGAkdot*yk + xkdot_dash*cos(OMEGAk) - tmp*sin(OMEGAk);
        ykvel = OMEGAkdot*xk + xkdot_dash*sin(OMEGAk) + tmp*cos(OMEGAk);
        zkvel = yk_dash*cos(ik)*ikdot + ykdot_dash*sin(ik);

        mu_div_rk3 = - mu/(rk*rk*rk);
        tmp2 = mu_div_rk3 + OMEGAedot*OMEGAedot;

        xkacc = tmp2*xk + 2.0*ykvel*OMEGAedot;
        ykacc = tmp2*yk - 2.0*xkvel*OMEGAedot;
        zkacc = mu_div_rk3*zk;






        %SVPosition(SV,:) = [xk  yk zk];

            SV_X_Data = xk;
            SV_Y_Data = yk;
            SV_Z_Data = zk;
            SV_T_Data = delta_ts;

        SV_Xvel_Data = xkvel;
        SV_Yvel_Data = ykvel;
        SV_Zvel_Data = zkvel;
        SV_Tvel_Data = delta_tsdot;

        SV_Xacc_Data = xkacc;
        SV_Yacc_Data = ykacc;
        SV_Zacc_Data = zkacc;
        SV_Tacc_Data = delta_tsdotdot;


        ValidData = 1;


    else  % DataValid = 0, no ephemeris was found for the satellite within the right time
        
        
            SV_X_Data = 0;
            SV_Y_Data = 0;
            SV_Z_Data = 0;
            SV_T_Data = 0;

        
        
        
        
        SV_Xvel_Data = 0;
        SV_Yvel_Data = 0;
        SV_Zvel_Data = 0;
        SV_Tvel_Data = 0;

        SV_Xacc_Data = 0;
        SV_Yacc_Data = 0;
        SV_Zacc_Data = 0;
        SV_Tacc_Data = 0;
        ValidData = 0;
        % return;

    end

%verified ephemeris and satellite selection from navigation data :
    %6/9/05 TB
    %successfully (1)ignores satellite ephemeris which has bad health flag
    %(2) doesn't use ephemeris more than 2 hours 5 min old
    %(3) ignores satellites which aren't in the ephemeris
    %(4) in all failed cases returns 0 0 0 0 for position and 0 for
    %ValidData



