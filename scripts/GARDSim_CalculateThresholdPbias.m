%$Id: GARDSim_CalculateThresholdPbias.m 2724 2009-07-26 12:48:03Z greerd $
%This code calculates the values for a values (used to calculate) threshold and pbias
%for different degrees of freedom. and PfalseAlarm and Pmiss 
%input the PFalseAlarm and Pmiss and the range of degrees' of freedoms you
%would like to calculate the a and pbias for.
%outputs are stored in the a vector and the lambdatrue vector

feature accel on;




%----------------------------------------------------------------------
%Calculate 'a' and 'pbias 'parameters for parity RAIM scheme
%---------------------------------------------------------------------
%determine the threshold using chi-square statistics.
%Chi-square degrees of freedom , d = n - 4 = 2.

%d = N-4;

%INPUTS
%------------------------
% PFalseAlarm = 2e-5;  % 10e-5/hr over 8 seconds TTA from WAAS MOPS
% Pmiss = 0.04799;
% drange = [1:3]; %this is the range of degrees of freedom to calculate for, note that is must always start at 1 for use with GARDsim_RAIMParity.m
% %------------------------



PFalseAlarm = 0.333e-6;  % 10e-5/hr over 8 seconds TTA from WAAS MOPS
Pmiss = 0.5e-3;
drange = [1:10];


numds = size(drange);
numds2 = numds(2);


for j = 1:numds2


%calculate a
a(j) = chi2inv(1-PFalseAlarm,drange(j));
%normalised Td
Td_norm = sqrt(a(j));
%Unnormalised Td
%Td = SigmaS*sqrt(a(numds));


%-----------------------------
%determine protection radius. 
%-----------------------------
%First calculate pbias



%X = ncx2inv(P,V,DELTA)
%V is degrees of freedom
%P is probabilities
%lambda is noncentrality parameter



%Iterate until the right value for noncentrality parameter is found. 
%This will be when the integral = 0.001

%Initial guess of lambda
e = 0.5;
lambda = 10;
i = 1;
while abs(e) > 5e-5
    
    x = (0:0.001:a(j));
    p1 = ncx2pdf(x,drange(j),lambda);   
    zt = trapz(p1); %is this integration method accurate enough?
    z= zt*0.001;  %multiply by the spacing , i've chosen 0.01, seems to be accurate enough.
    
    e = z - Pmiss
    if abs(e) > 0.5

        if e>0
            lambda = lambda + 3;  %increase if e positive
        else
         lambda = lambda - 3;     %decrease if e negative

        end

    end    
    if abs(e) < 0.5 & abs(2) > 0.0005

        if e>0
            lambda = lambda + 0.1;  %increase if e positive
        else
         lambda = lambda - 0.1;     %decrease if e negative

        end

    end    
    
    if abs(e) < 0.0005 & abs(e) > 0.0002

        if e>0
            lambda = lambda + 0.02;  %increase if e positive
        else
         lambda = lambda - 0.02;     %decrease if e negative

        end
    end


    if abs(e) <= 0.0002 & abs(e) > 0.00001

        if e>0
            lambda = lambda + 0.001;  %increase if e positive
        else
         lambda = lambda - 0.001;     %decrease if e negative

        end      

    end


    if abs(e) <= 0.00001 & abs(e) > 0.000003

        if e>0
            lambda = lambda + 0.0001;  %increase if e positive
        else
         lambda = lambda - 0.0001;     %decrease if e negative

        end      

    end

    if abs(e) <= 0.000003
        if e>0
            lambda = lambda + 0.000001;  %increase if e positive
        else
         lambda = lambda - 0.000001;     %decrease if e negative

        end      
    end        

    
%     if abs(e) <=0.0000001
%         break;
%     end
    
    
    lambda
    
    lambdacheck(i) = lambda;
    echeck(i) = abs(e);
    
    i = i+1;
    
    
    
end  %end while


lambdatrue(j) = lambda;
j

end  %end for j loop.






