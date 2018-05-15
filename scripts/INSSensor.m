
%=======================================================
%INERTIAL NAVIGATION SOLUTION
%=======================================================

%this is run at 50Hz, simulating an INS output at 50Hz
for i = startepochHighRate:endepochHighRate

    ins_dt = 1/100;

    if i == startepochHighRate

        %start with the truth

        TMatrix_ECEF2NED = T_ECEF2NED(Lat_truth(i), Lon_truth(i));

        VelocityNED = TMatrix_ECEF2NED*[Xvel_truth(i),Yvel_truth(i),Zvel_truth(i)]';

        INS_state(1,i) = Quaternions_truth(1,i);
        INS_state(2,i) = Quaternions_truth(2,i);
        INS_state(3,i) = Quaternions_truth(3,i);
        INS_state(4,i) = Quaternions_truth(4,i);
        INS_state(5,i) = VelocityNED(1); %vn
        INS_state(6,i) = VelocityNED(2); %ve
        INS_state(7,i) = VelocityNED(3); %vd
        INS_state(8,i) = Lat_truth(i); %lat
        INS_state(9,i) = Lon_truth(i); %lon
        INS_state(10,i) = Hgt_truth(i); %hgt

        %normalise quaternions truth
        omega_x = p_INS_50Hz(i) ;
        omega_y = q_INS_50Hz(i);
        omega_z = r_INS_50Hz(i);


        A_xb = ax_b_INS_50Hz(i)  ;
        A_yb = ay_b_INS_50Hz(i)  ;
        A_zb = az_b_INS_50Hz(i) ;

        A_b = [A_xb, A_yb, A_zb];
        omega_b = [omega_x, omega_y, omega_z];

        
        %gravity = -9.80;

        %Get Gravity estimate (use truth at the moment)
       % g = GravityTruth(i); %note that it is correct to use -ve g i.e.-9.8

        g = 9.81

        %         llh_dot = (INS_state(8:10,i) - INS_state(8:10,i-1))/0.02;
        %
        %            [INS_state(:,i+1)] = INS_Mechanization2(INS_state(:,i), A_b, omega_b, g, ins_dt,llh_dot);
        %

        %INS_state(:,i+1) = INS_state(:,i);



        %just put this in here temporarily to check whether aerosims vn ve and vd calculations are correct (using the lat and lon dots.
        TMatrix_ECEF2NED = T_ECEF2NED( Lat_truth(i+1),Lon_truth(i+1));
        VelocityNED = TMatrix_ECEF2NED*[Xvel_truth(i+1),Yvel_truth(i+1),Zvel_truth(i+1)]';


        INS_state(1,i+1) = Quaternions_truth(1,i+1);
        INS_state(2,i+1) = Quaternions_truth(2,i+1);
        INS_state(3,i+1) = Quaternions_truth(3,i+1);
        INS_state(4,i+1) = Quaternions_truth(4,i+1);
        INS_state(5,i+1) = VelocityNED(1); %vn
        INS_state(6,i+1) = VelocityNED(2); %ve
        INS_state(7,i+1) = VelocityNED(3); %vd
        INS_state(8,i+1) = Lat_truth(i+1); %lat
        INS_state(9,i+1) = Lon_truth(i+1); %lon
        INS_state(10,i+1) = Hgt_truth(i+1); %hgt



        %normalise quaternions
        [INS_state(1,i+1),INS_state(2,i+1),INS_state(3,i+1),INS_state(4,i+1)] = Normalise_Quat(INS_state(1,i+1),INS_state(2,i+1),INS_state(3,i+1),INS_state(4,i+1));





    else

        omega_x = p_INS_50Hz(i) ;
        omega_y = q_INS_50Hz(i);
        omega_z = r_INS_50Hz(i);


        A_xb = ax_b_INS_50Hz(i)  ;
        A_yb = ay_b_INS_50Hz(i)  ;
        A_zb = az_b_INS_50Hz(i) ;

        A_b = [A_xb, A_yb, A_zb];
        omega_b = [omega_x, omega_y, omega_z];


        %normalise quaternions

        %normalise quaternions
        [INS_state(1,i),INS_state(2,i),INS_state(3,i),INS_state(4,i)] = Normalise_Quat(INS_state(1,i),INS_state(2,i),INS_state(3,i),INS_state(4,i));




        %Get Gravity estimate (use truth at the moment)
        g = GravityTruth(i); %note that it is correct to use -ve g i.e.-9.8

        llh_dot = (INS_state(8:10,i) - INS_state(8:10,i-1))/0.01;

        %llh_dot = [0 0 0];

        [INS_state(:,i+1)] = INS_Mechanization2(INS_state(:,i), A_b, omega_b, g, ins_dt,llh_dot);


        %normalise quaternions
        [INS_state(1,i),INS_state(2,i),INS_state(3,i),INS_state(4,i)] = Normalise_Quat(INS_state(1,i),INS_state(2,i),INS_state(3,i),INS_state(4,i));


    end


    Latpos_INS(i) = INS_state(8,i);
    Lonpos_INS(i) = INS_state(9,i);
    Hgtpos_INS(i) = INS_state(10,i);

    Nvel_INS(i) = INS_state(5,i);
    Evel_INS(i) = INS_state(6,i);
    Dvel_INS(i) = INS_state(7,i);

end %end for




%convert to ECEF and convert quaternions

for i = startepochHighRate:endepochHighRate

    PositionECEF_INS = LLH2ECEF(Latpos_INS(i),Lonpos_INS(i),Hgtpos_INS(i));

    Xpos_INS50Hz(i) = PositionECEF_INS(1);
    Ypos_INS50Hz(i) = PositionECEF_INS(2);
    Zpos_INS50Hz(i) = PositionECEF_INS(3);



    quat(1) = INS_state(1,i);
    quat(2) = INS_state(2,i);
    quat(3) = INS_state(3,i);
    quat(4) = INS_state(4,i);

    [euler] = QuatToEuler(quat);


    phi_INS50Hz(i)= euler(1);
    theta_INS50Hz(i)= euler(2);
    psi_INS50Hz(i) = euler(3);




end

%INS at 1 Hz

%aGPS = 1
for i = startepoch:endepoch

    %read in the data and perform necessary conversions
    %this should be the GPS inputs but use the truth for now.
    %ECEF positions
    %Note that first row is simulation time data thats why the index starts at
    %2 in the following

    Xpos_INS(i) = Xpos_INS50Hz(i*100+1);
    Ypos_INS(i) = Ypos_INS50Hz(i*100+1);
    Zpos_INS(i) = Zpos_INS50Hz(i*100+1);


    %INS measured NED velocities
    Vn_INS(i) = INS_state(5,i*100+1);
    Ve_INS(i) = INS_state(6,i*100+1);
    Vd_INS(i) = INS_state(7,i*100+1);


    %convert NED velocity to XYZ vel.

    TMatrix = T_ECEF2NED(Latpos_INS(i*100+1),Lonpos_INS(i*100+1));
    Vel_INS_ecef = TMatrix'*[Vn_INS(i), Ve_INS(i), Vd_INS(i)]';

    Xvel_INS(i) = Vel_INS_ecef(1);
    Yvel_INS(i) = Vel_INS_ecef(2);
    Zvel_INS(i) = Vel_INS_ecef(3);


    % euler estimates
    phi_INS(i) = phi_INS50Hz(i*100+1);
    theta_INS(i) = theta_INS50Hz(i*100+1);
    psi_INS(i) = psi_INS50Hz(i*100+1);


    %raw output from INS
    ax_b_INS1Hz(i) = ax_b_INS_50Hz(i*100+1);
    ay_b_INS1Hz(i) = ay_b_INS_50Hz(i*100+1);
    az_b_INS1Hz(i) = az_b_INS_50Hz(i*100+1);


    omega_x1Hz(i) = p_INS_50Hz(i*100+1);
    omega_y1Hz(i) = q_INS_50Hz(i*100+1);
    omega_z1Hz(i) = r_INS_50Hz(i*100+1);




    q0_INS1Hz(i) = INS_state(1,i*100+1);
    q1_INS1Hz(i) = INS_state(2,i*100+1);
    q2_INS1Hz(i) = INS_state(3,i*100+1);
    q3_INS1Hz(i) = INS_state(4,i*100+1);


    Latpos_INS1Hz(i) =  Latpos_INS(i*100+1);
    Lonpos_INS1Hz(i) = Lonpos_INS(i*100+1);
    Hgtpos_INS1Hz(i) = Hgtpos_INS(i*100+1);

end


%end INS sensor code
