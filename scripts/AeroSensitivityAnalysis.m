
clear all
feature accel on
warning off

    
 %simulate rate in hz , can't be greater than 50, must be a whole number
 Rate = 100; 
    
  %Use state 0 of the random number generator

randn('state',0);
  
 %startepoch = 10*Rate;  %start at 10 seconds
 
startepoch = 1;  %start epoch is 10 for 1 Hz Rate and 19 for 2 Hz
 %data is at 50 Hz so take every 50th for 1 second  
endepoch = 5000;
% endepoch = 386;
%if you want at 1 hz a = 50/1 = 50, therefore you use i*a - (a-1), ie i*50-49
%if you want at 5 hz you do a = 50/5 = 10 therefore you use i*10-9
%for 2 hz you use a = 50/2 = 25, therefore you use i*25-24
%etc
step_size = 1/Rate;  
 
%=====================================================================
%Gravity Model
%=====================================================================
% Calculate Earth Parameters (WGS-84)

%calculate gravity based on true position

%[Gravity] = Earth_Gravity(x_LLH(1:3)',0.0818191908426^2,9.7803267714,0.00193185138639);

%just advance it by 1 save me getting rid of all the -1's in the code below
%startepoch = startepoch+1;

Gravity = 9.81;

%=====================================================================
%Atmospheric Parameters
%=====================================================================

%Standard Atmosphere in SI units

%Altitude % ASsume GA aircraft operating betweeen 0 and 15,000 feet. 

%Temperature

%[P_atm T_atm p_atm a_atm] = Atmosphere_STDAtmosphere(x_LLH(3),287.0531,1.40,9.80665);

% 
% P_atm = AtmosTruth(1,startepoch); 
% 
% T_atm = AtmosTruth(2,startepoch) ;
% 
% p_atm = AtmosTruth(3,startepoch) ;
% 
% a_atm = AtmosTruth(4,startepoch) ;
% 
% 
% %absolute ambient temperature in Kelvins
% T_amb = T_atm; %K
% R_gas = 287.0531,; %m^2/(K.s^2)
% 
% %Pressure
% 
% %Density
% rho = p_atm ; %kg/m^3 %density
% 
% Atmos = [T_amb R_gas rho];   %just use these parameters for now. 
 
%use truth for pqr

  NumberPoints = 1;        
    
  N = NumberPoints;
  var = 0.001;
        
x = sqrt(var)*randn(N,1);  %created 100 indepndant random variable samples from the distribution.
      
%the constants are the aircraft parameters, gravity, environmental
%parameters

%changing one value at a time. 
            
%changing aircraft parameters

%Initial Conditions, which will depend upon previous epoch
V_T0 = 40;
beta0 = 0;
alpha0 = 0;
phi0 = 0;
theta0 = 0;
psi0 = 0;
p0 = 0;
q0 = 0;
r0 = 0
pn0 = 0; %these arent used as inputs
pe0 = 0;  
h0 = 0;
q00 = 1;
q10 = 0;
q20 = 0;
q30 = 0;
Lat0 = 27*pi/180;
Lon0 = 153*pi/180;
Hgt0 = 2000;
u0 = 40;
v0 = 0;
w0 = 0;


        
    
for N = 1:NumberPoints

%x0(startepoch,:) = [V_T0+N*2, beta0+N*0.5*pi/180 ,alpha0+N*0.5*pi/180, phi0+N*0.5*pi/180, theta0+N*0.5*pi/180, psi0+N*0.5*pi/1802, p0+N*0.03*pi/180,q0+N*0.03*pi/180,r0+N*0.03*pi/180,pn0, pe0, h0, q00, q10, q20, q30, Lat0+N*2/Rnormal, Lon0+N*2/Rnormal, Hgt0+N*2,u,v,w];    
   

x0(startepoch,:) = [V_T0, beta0 ,alpha0, phi0, theta0, psi0, p0,q0,r0,pn0, pe0, h0, q00, q10, q20, q30, Lat0, Lon0, Hgt0,u0,v0,w0,0,0,0];    
    

%x0(startepoch,:) = [2, 2 ,2, 2, 2, 2, 2,2,2,2, 2,2, 2, 2, 2, 2, 2, 2, 2,2,2,2,x(N),0,0];    
    

for i = startepoch:10         
 
 
 delta_x_control(i,:) = [0, 0,0,0]; %this throttle needs to be scaled, not suree what it means at the moment, but probably needs to represent a change in thrust, since it is multiplied by the Xdelta_t term in the model    %control inputs will be the current control inputs at i
   
     

  %calculate MACH number for next time around
      M0 = x0(i,1)/340.29;
      
      LatModel(i) = x0(i,17);
      LonModel(i) = x0(i,18);
      HgtModel(i) = x0(i,19);  
      
      [ECEFModel] = LLH2ECEF(LatModel(i),LonModel(i),HgtModel(i));

        %current position
        XposModel(i) = ECEFModel(1);
        YposModel(i) = ECEFModel(2);
        ZposModel(i) = ECEFModel(3);
      
        
      AeroCoeffCorrections = zeros(1,6);
      Parameter_Noise = zeros(1,6);
      Wind_NED = [0,0,0];
      GravityTruth = 9.81;
      ts = 0.02;;
                                      
                
    [dx] = AerodynamicModelStevensNewSensivity(M0,GravityTruth,x0(i,:),delta_x_control(i,:));
    
          
     x0(i+1,1:22) = x0(i,1:22) + ts*dx(1:22);    
          
     x0(i+1,22:24) =[dx(23),dx(24),dx(25)];       
      
     %Grab the accelerations and roll rates      
      
      u_dot(i,N) = dx(20);
      v_dot(i,N) = dx(21);
      w_dot(i,N) = dx(22);
      
      ps(i,N) = dx(23);
      qs(i,N) = dx(24);
      rs(i,N) = dx(25);            
      

end  %end i

x0ALL(:,:,N) =  x0(:,:) ;

end %end N


%PLOTS to show the effect of different inputs on the calculated
%accelerations and p q r's

%hist(x(N));




%Plot input



%plot outputs for all parameters




























