function  [UKF_x_hat_kup,UKF_Px_kup,xs,Ws]  = GARD_SigmaPointPropagation(...
    UKF_x_hat_kminus,...
    UKF_Px_kminus,...
    UKF_Params,...
    Sensor_Params,...
    xs,...
    Ws,...
    RegenerateSPs)


ins_dt = Sensor_Params.dt;
A_xb_0 = Sensor_Params.Axb;
A_yb_0 = Sensor_Params.Ayb;
A_zb_0 = Sensor_Params.Azb;
omega_x_0 = Sensor_Params.p;
omega_y_0 = Sensor_Params.q;
omega_z_0 = Sensor_Params.r;



if RegenerateSPs == 1
    [xs,Ws] = GARD_GenerateSigmaPoints(...
                UKF_x_hat_kminus,...
                UKF_Px_kminus,...
                UKF_Params.Q,...
                UKF_Params.R,...
                UKF_Params.alpha,...
                UKF_Params.beta,...
                UKF_Params.kapa);


end 

if ~isreal(xs)
   disp('warning - imag sigma points'); 
   xs = real(xs);
end


    xs_0 = xs(:,1);
    xs_i = xs(:,2:size(xs,2));
    W_0_m = Ws(1);
    W_0_c = Ws(2);
    W_i_m = Ws(3);
    W_i_c = Ws(4);
    
Na = size(xs,1);
NumberStates = size(UKF_x_hat_kminus,1);
NumberMeasurements = size(UKF_Params.R,1);
ProcessNoiseStates = size(UKF_Params.Q,1);

        % propogate sigma-points through system dynamics
        for i=0:2*Na
            if i==0  % the mean sigma point
                %% add process noise sigma points to sensor inputs
                Acc_in(1) = A_xb_0 - xs_0(11) + xs_0(NumberStates+1);
                Acc_in(2) = A_yb_0 - xs_0(12) + xs_0(NumberStates+2);
                Acc_in(3) = A_zb_0 - xs_0(13) + xs_0(NumberStates+3);
                Omega_in(1) = omega_x_0 - xs_0(14) + xs_0(NumberStates+4);
                Omega_in(2) = omega_y_0 - xs_0(15) + xs_0(NumberStates+5);
                Omega_in(3) = omega_z_0 - xs_0(16) + xs_0(NumberStates+6);
                
                g = GravityModel(xs_0(1:3));
                [xs_0_k,output] = GARD_INSMechanisation(xs_0,ins_dt,Acc_in,Omega_in,g);
                
                % update clock
                xs_0_k(17) = xs_0(17) + xs_0(18)*ins_dt + xs_0(NumberStates+7);
                xs_0_k(18) = xs_0(18) + xs_0(NumberStates+8);
                                                     
%                 % update random walk process for sensor biases
                    
                 % 1st order GM model
%                 xs_0_k(11) = (1.0-x_accel_beta)*xs_0(11) + xs_0(19);
%                 xs_0_k(12) = (1.0-y_accel_beta)*xs_0(12) + xs_0(20);
%                 xs_0_k(13) = (1.0-z_accel_beta)*xs_0(13) + xs_0(21);
%                 xs_0_k(14) = (1.0-x_gyro_beta)*xs_0(14) + xs_0(22);
%                 xs_0_k(15) = (1.0-x_gyro_beta)*xs_0(15) + xs_0(23);
%                 xs_0_k(16) = (1.0-x_gyro_beta)*xs_0(16) + xs_0(24);
%                 
%                 xs_0_k(11) = GaussMarkov_Process2(xs_0(11), x_accel_beta,xs_0(19)/ins_dt,ins_dt);
%                 xs_0_k(12) = GaussMarkov_Process2(xs_0(12), y_accel_beta,xs_0(20)/ins_dt,ins_dt);
%                 xs_0_k(13) = GaussMarkov_Process2(xs_0(13), z_accel_beta,xs_0(21)/ins_dt,ins_dt);
%                 xs_0_k(14) = GaussMarkov_Process2(xs_0(14), x_gyro_beta,xs_0(22)/ins_dt,ins_dt);
%                 xs_0_k(15) = GaussMarkov_Process2(xs_0(15), y_gyro_beta,xs_0(23)/ins_dt,ins_dt);
%                 xs_0_k(16) = GaussMarkov_Process2(xs_0(16), z_gyro_beta,xs_0(24)/ins_dt,ins_dt);
                
                 % random walk model
                xs_0_k(11) = xs_0(11) + xs_0(NumberStates+1)*ins_dt;
                xs_0_k(12) = xs_0(12) + xs_0(NumberStates+2)*ins_dt;
                xs_0_k(13) = xs_0(13) + xs_0(NumberStates+3)*ins_dt;
                xs_0_k(14) = xs_0(14) + xs_0(NumberStates+4)*ins_dt;
                xs_0_k(15) = xs_0(15) + xs_0(NumberStates+5)*ins_dt;
                xs_0_k(16) = xs_0(16) + xs_0(NumberStates+6)*ins_dt;
                % random constant model
%                 xs_0_k(11) = xs_0(11);
%                 xs_0_k(12) = xs_0(12);
%                 xs_0_k(13) = xs_0(13);
%                 xs_0_k(14) = xs_0(14);
%                 xs_0_k(15) = xs_0(15);
%                 xs_0_k(16) = xs_0(16);

                % propogate process noise for sensor bias
                 for m=19:24
                     xs_0_k(m) = xs_0(m);
                 end
                
                % propogate process noise for clock states
                xs_0_k(25) = xs_0(25);
                xs_0_k(26) = xs_0(26);
                
                % propogate measurement noise
                for m=27:Na
                    xs_0_k(m) = xs_0(m);
                end
                
                                                     
            else

                %% add process noise sigma points to sensor inputs
                Acc_in(1) = A_xb_0 - xs_i(11,i) + xs_i(NumberStates+1,i);
                Acc_in(2) = A_yb_0 - xs_i(12,i) + xs_i(NumberStates+2,i);
                Acc_in(3) = A_zb_0 - xs_i(13,i) + xs_i(NumberStates+3,i);
                Omega_in(1) = omega_x_0 - xs_i(14,i) + xs_i(NumberStates+4,i);
                Omega_in(2) = omega_y_0 - xs_i(15,i) + xs_i(NumberStates+5,i);
                Omega_in(3) = omega_z_0 - xs_i(16,i) + xs_i(NumberStates+6,i);
                
                g = GravityModel(xs_i(1:3,i));
                [xs_i_k(1:10,i),output] = GARD_INSMechanisation(xs_i(:,i),ins_dt,Acc_in,Omega_in,g);
                
                
% 
                % update clock
                xs_i_k(17,i) = xs_i(17,i) + xs_i(18,i)*ins_dt + xs_i(NumberStates+7,i);
                xs_i_k(18,i) = xs_i(18,i) + xs_i(NumberStates+8,i);
                                                     
                % update random walk process for sensor biases
%                  xs_i_k(11,i) = (1.0-x_accel_beta)*xs_i(11,i) + xs_i(19,i);
%                  xs_i_k(12,i) = (1.0-y_accel_beta)*xs_i(12,i) + xs_i(20,i);
%                  xs_i_k(13,i) = (1.0-z_accel_beta)*xs_i(13,i) + xs_i(21,i);
%                  xs_i_k(14,i) = (1.0-x_gyro_beta)*xs_i(14,i) + xs_i(22,i);
%                  xs_i_k(15,i) = (1.0-y_gyro_beta)*xs_i(15,i) + xs_i(23,i);
%                  xs_i_k(16,i) = (1.0-z_gyro_beta)*xs_i(16,i) + xs_i(24,i);
% 
%                 xs_i_k(11,i) = GaussMarkov_Process2(xs_i(11,i), x_accel_beta,xs_i(19,i)/ins_dt,ins_dt);
%                 xs_i_k(12,i) = GaussMarkov_Process2(xs_i(12,i), y_accel_beta,xs_i(20,i)/ins_dt,ins_dt);
%                 xs_i_k(13,i) = GaussMarkov_Process2(xs_i(13,i), z_accel_beta,xs_i(21,i)/ins_dt,ins_dt);
%                 xs_i_k(14,i) = GaussMarkov_Process2(xs_i(14,i), x_gyro_beta,xs_i(22,i)/ins_dt,ins_dt);
%                 xs_i_k(15,i) = GaussMarkov_Process2(xs_i(15,i), y_gyro_beta,xs_i(23,i)/ins_dt,ins_dt);
%                 xs_i_k(16,i) = GaussMarkov_Process2(xs_i(16,i), z_gyro_beta,xs_i(24,i)/ins_dt,ins_dt);                 
         
                 xs_i_k(11,i) = xs_i(11,i) + xs_i(NumberStates+1,i)*ins_dt;
                 xs_i_k(12,i) = xs_i(12,i) + xs_i(NumberStates+2,i)*ins_dt;
                 xs_i_k(13,i) = xs_i(13,i) + xs_i(NumberStates+3,i)*ins_dt;
                 xs_i_k(14,i) = xs_i(14,i) + xs_i(NumberStates+4,i)*ins_dt;
                 xs_i_k(15,i) = xs_i(15,i) + xs_i(NumberStates+5,i)*ins_dt;
                 xs_i_k(16,i) = xs_i(16,i) + xs_i(NumberStates+6,i)*ins_dt;

%                  xs_i_k(11,i) = xs_i(11,i);
%                  xs_i_k(12,i) = xs_i(12,i);
%                  xs_i_k(13,i) = xs_i(13,i);
%                  xs_i_k(14,i) = xs_i(14,i);
%                  xs_i_k(15,i) = xs_i(15,i);
%                  xs_i_k(16,i) = xs_i(16,i);

                 % propogate process noise for sensor bias
                 for m=NumberStates+1:NumberStates+6
                     xs_i_k(m,i) = xs_i(m,i);
                 end
                 
                 % propogate process noise for clock states
                 xs_i_k(NumberStates+7,i) = xs_i(NumberStates+7,i);
                 xs_i_k(NumberStates+8,i) = xs_i(NumberStates+8,i);
               
                % propogate measurement noise
                 for m=NumberStates+ProcessNoiseStates:Na 
                     xs_i_k(m,i) = xs_i(m,i);
                 end
            end
        end

        
        % update xs_i
        xs_0 = xs_0_k;
        xs_i = xs_i_k;
        
        % calculate the predicted 
        % mean and state covariance
         xs = [xs_0_k,xs_i_k];
         [UKF_x_hat_kup,UKF_Px_kup] = GARD_GenerateMeanFromSigmaPoints(xs,Ws,NumberStates);
