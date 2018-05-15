function [CurveFit, BestFit,  normX, Xest] = LeastSquaresBestFitWindGMProcess(Signal, Order, dT, FindBestFit,randU);

feature accel on
%Least Squares Data Fitting First Order Gauss Markov Process

%GM PRocess series of the form X(k) = exp(-Beta*dT)*X(k-1) +
%sqrt((WN_Std_dev^2)*(1-exp(-2*Beta*dT)))*randU




%by Troy Bruggemann 27 Feb 2008


%
% inputs:

% Signal - the data
% Order - the order you want fitted to the data if it isn't finding the  'best fit'
% FindBestFit - if 1, finds the best fit to the data, if 0, just calculates the fourier series with order as specified by 'Order'


%outputs

% CurveFit - the fitted fourier series to the data
% BestFit - the order of the best fitting fourier series to the data
% normX - the tolerance of the residuals
% Xest - the Fourier series parameters :   a0, a1, a2..aN, b1, b2 ...bN



%Give a Fourier Series fit or have it detect which is the best fit order
%Fourier Series to the data

% 1. enter initial guess for fitting parameters
% 2.  Test the initial guess against the data
% 3.  Use least squares to refine the initial guess until is minimised
% 4 . The fit curve is then the best fit of the data in least squares sense
% 5. Repeats the process for each order and finds the best ordered fourier
% series that best fits the data, determined by the sum of squares of
% residuals





BestFitTemp = Order;
clear y Xest f M Minit Resid DeltaX normx;

warning off;%turn off annoying warnings for singular inverse matrices

feature accel on

%find the best fit curve


if FindBestFit == 1

    for N = 1:30

        clear y Xest f M Minit Resid DeltaX normx;

        y = Signal;
        numpoints = size(y);
        sizey = numpoints(2);


        nummeas = numpoints(1);

        x = 1:sizey;
        %order of polynomial to fit to data
        %N = 6;
        
        

        %Vector of initial guess of parameters a0 a1 etc..

        %Xest = zeros(1,N+1);
        %preallocate arrays

        a0 = 0;

        a1 = zeros(1,N);  %first estimate is zero

        b1 = zeros(1,N);



        Xest = [a0, a1, b1];   %state vector

        Minit = zeros(sizey,2*N+1);

        f = zeros(1,sizey);




        %Fourier series


        for i = 1:100

            a0 = Xest(1);
            a1 = Xest(2:N+1);
            b1 = Xest(N+2:2*N+1);



            for x = 1:sizey
                f(x) = a0/2;
                for j = 1:N
                    % f(x) = f(x) + Xest(N-j+2)*x.^(j-1);

                    %f(x) =  Xest(1)*x.^5 + Xest(2)*x.^4 + Xest(3)*x.^3  + Xest(4)*x.^2 + Xest(5)*x.^1 + Xest(6) ;
                    %construct polynomial


                    f(x) = f(x) +  a1(j)*cos(j*x) + b1(j)*sin(j*x);

                end
            end






            for index = 1:nummeas


                Resid((index-1)*length(f)+1:index*length(f)) = y(index,:)-f;

                %count = count+length(f);
            end



            for x = 1:sizey

                Minit(x,1) = 1/2;

                for j = 2:N+1
                    Minit(x,j) =  cos((j-1)*x);
                    Minit(x,j+N) = sin((j-1)*x);

                end

            end


            for index = 1:nummeas

                M((index-1)*length(f)+1:index*length(f),:) = Minit(:,:);
            end




            A = M'*M;
            b = M'*Resid';

            DeltaX = A\b;


            %DeltaX = inv(M'*M)*M'*Resid';
            normX = norm(DeltaX);
            if normX < 1e-3;
                break;
            else
                Xest = Xest +DeltaX';
            end


        end

        for x = 1:sizey
            f(x) = a0/2;
            for j = 1:N
                f(x) = f(x) +  a1(j)*cos(j*x) + b1(j)*sin(j*x);
                %f(x) =  Xest(1)*x.^5 + Xest(2)*x.^4 + Xest(3)*x.^3  + Xest(4)*x.^2 + Xest(5)*x.^1 + Xest(6) ;
                %construct polynomial
            end
        end
        %Order(N) = sum((y-f).^2);

        for index = 1:nummeas

            ydiff(index,:) = y(index,:)-f;

        end

        Order(N) = sum(sum(ydiff(:,:)'.^2));

        %Order(N) = sum( ydiff(1,:).^2+ ydiff(2,:).^2 + ydiff(3,:).^2 + ydiff(4,:).^2);

    end

    %evaluate choice of order

    Minimum = min(Order);
    BestFitTemp = find(Order == Minimum);
    BestFitTemp = BestFitTemp(1);  %in case it returns two or more orders at the same minimum, choose the lowest order one



end






clear y Xest f M Minit Resid DeltaX normx




y = Signal;

numpoints = size(y);
sizey = numpoints(2);

nummeas = numpoints(1);

nummeas = 2; %extra set of data of random numbers
x = 1:sizey;


%generate random numbers

rn = randn(1,sizey);


rnmeas = randn(1,sizey); %generate another set as measurements

rnmeas2 = randn(1,sizey);
 
rnmeas3 = randn(1,sizey);

%rn = randU;

%order of polynomial to fit to data
N = BestFitTemp;


%beta = 1/300;   %these cant be zero otherwise get NaN due to M(x,2);
%sigma = 5;




%use autocorrelation function to get initial value

maxlags = sizey;
[autocorr,lags] = xcorr(y,maxlags);
autocorr_plot = autocorr/length(y);



sigma = sqrt( autocorr_plot(length(y)+1));


%the peak of the autocorrelation function is the variance





value = 0.368*sigma^2; 
A = autocorr_plot(maxlags+1:length(autocorr_plot));
index_Vec = find( value - 0.01*value < A & A <  value + 0.01*value);
ax_tau_check = A(min(index_Vec));
ax_tau = min(index_Vec);

beta = 1/ax_tau;








% beta = 1/100;   %these cant be zero otherwise get NaN due to M(x,2);
% sigma = 343.24245;
%randnum = 0.1;  %this must be between -3 and 3 (according to normal distribution with mean 0 and variance 1

%a1 = zeros(1,N);  %first estimate is zero

%b1 = zeros(1,N);
%Vector of initial guess of parameters a0 a1 etc..

%number of states will be N+1 , the +1 is because of the a0/2 term, since i
%estimate a0 as well.

Xest = [beta, sigma, rn];   %state vector
%Xest = [a0, a1, b1];   %state vector

%preallocate arrays

%Minit = zeros(sizey,2*N+1);


Minit = zeros(sizey,2+sizey);   % 2 states, beta and sigma
%Xest = zeros(1,N+1);
f = zeros(1,sizey);



for i = 1:100

    
    
    
    %beta and sigma must be positive
    
    
    Xest(1) = abs(Xest(1));
     Xest(2) = abs(Xest(2));

    beta = Xest(1);
    sigma = Xest(2);
    rn = Xest(3:length(Xest));  %these dont have to be positive
    
    
    
   % randnum = Xest(3);
%     a0 = Xest(1);
%     a1 = Xest(2:N+1);
%     b1 = Xest(N+2:2*N+1);

% 
% 
% if sigma < 0
%     
%     turd = 1
% 
% end
% 
% 
% if beta < 0
%     
%     turd = 1
% 
% end

% 
% 
% for x = 1:1000
%     
%     n(x) = 1/(sqrt(2*pi))*exp(-x/2);
% 
% end


  f(1) = sigma*rn(1)*dT; 
  
    %f(1) = sigma*randnum*dT; 
    for x = 2:sizey
       % f(x) = a0/2;
        
      
       % for j = 1:N
            % f(x) = f(x) + Xest(N-j+2)*x.^(j-1);

            %f(x) =  Xest(1)*x.^5 + Xest(2)*x.^4 + Xest(3)*x.^3  + Xest(4)*x.^2 + Xest(5)*x.^1 + Xest(6) ;
            %construct polynomial


           % f(x) = f(x) +  a1(j)*cos(j*x) + b1(j)*sin(j*x);
           
          % f(x) = f(x-1)*exp(-beta*dT) + sqrt((sigma^2)*(1-exp(-2*beta*dT)))*rn(x);
           
           f(x) = f(x-1)*exp(-beta*dT) + sigma*sqrt(1-exp(-2*beta*dT))*rn(x);

        %end
        
        
    end
    %Resid = y-f;



    %
    %
    %         Resid = y-f;
    %
    %
    %     for x = 1:sizey
    %
    %
    %         M(x,1) = 1/2;
    %
    %          for j = 2:N+1
    %          M(x,j) =  cos((j-1)*x);
    %
    %          M(x,j+N) = sin((j-1)*x);
    %
    %
    %          end
    %
    %
    %
    %     end


    %     count = 1;
    %     for index = 1:nummeas
    %
    %         Resid(1:length(f)) = y(1,:)-f;
    %         Resid(length(f)+1:2*length(f)) = y(2,:)-f;
    %         Resid(2*length(f)+1:3*length(f)) = y(3,:)-f;
    %         Resid(3*length(f)+1:4*length(f)) = y(4,:)-f;
    %
    %         %count = count+length(f);
    %     end


    
    %get new measurements from normal distribution
    
    
   
    
    
    Resid(1:length(f)) = y(1,:) - f;
    Resid(length(f)+1:2*length(f)) =  rnmeas - rn;
       Resid(2*length(f)+1:3*length(f)) =  rnmeas2 - rn;
        Resid(3*length(f)+1:4*length(f)) =  rnmeas3 - rn;
    
%     for index = 1:nummeas
% 
%         Resid((index-1)*length(f)+1:index*length(f)) = y(index,:)-f;
%         
%         
%         
%         
%         %count = count+length(f);
%     end
    
    
    
   


    %initial values
    Minit(1,1) = 0;
    Minit(1,2) = rn(1)*dT;
    Minit(1,3:sizey+2) = sigma*dT;
    
    df_dbetaprev = 0; 
    df_dsigmaprev = 0;
    df_drnprev(1:sizey) = 0;
    for x = 2:sizey

        Minit(x,1) = exp(-beta*dT)*(df_dbetaprev - f(x-1)*dT) + (sigma*dT*exp(-2*beta*dT))/sqrt(1-exp(-2*beta*dT))*rn(x);
        
         Minit(x,2) = df_dsigmaprev*exp(-beta*dT) + sqrt(1-exp(-2*beta*dT))*rn(x);
         
         
         Minit(x,3:sizey+2) = df_drnprev(1:sizey).*exp(-beta*dT) + sigma*sqrt(1-exp(-2*beta*dT));
         
         df_dbetaprev =Minit(x,1);
         df_dsigmaprev = Minit(x,2);
         
          df_drnprev(1:sizey) = Minit(x,3:sizey+2);
%          %protect from singularity, if is NaN due to denominator being 0 just set to large value, 
%          if Minit(x,2) > 500000000 
%              
%              Minit(x,2) = 100000;
%          end        
        
%         for j = 2:N+1
%             Minit(x,j) =  cos((j-1)*x);
%             Minit(x,j+N) = sin((j-1)*x);
% 
%         end

    end
        
    
    
    %Minit(sizey+1:sizey+1 + sizey, 3:sizey+2) = eye(sizey, sizey);   

%     for pup = sizey+1:4*sizey
%         
%         for mup = 3:sizey+2
%             
%              Minit(pup, mup) = 1; 
%         end
%         
%     end
    
    
    Minit(sizey+1:2*sizey,3:sizey+2) = eye(1000,1000);
    
       Minit(2*sizey+1:3*sizey,3:sizey+2) = eye(1000,1000);
       
         Minit(3*sizey+1:4*sizey,3:sizey+2) = eye(1000,1000);
        
           
    
    
            
    
   % Minit(sizey+1:4*sizey, 3:sizey+2) = eye(sizey*3,sizey*3);
    
    M = Minit;

%     for index = 1:nummeas
% 
%         M((index-1)*length(f)+1:index*length(f),:) = Minit(:,:);
%     end



    A = M'*M;
    b = M'*Resid';
    DeltaX = A\b;
    %DeltaX = inv(M'*M)*M'*Resid';
    normX = norm(DeltaX);
    if normX < 1e-3;
        break;
    else
        Xest = Xest +DeltaX';
        
        
        Xest(1:2) = abs(Xest(1:2));
    end
    
    
     f(1) = sigma*rn(1)*dT; 
    for x = 2:sizey
       % f(x) = a0/2;
        
      
       % for j = 1:N
            % f(x) = f(x) + Xest(N-j+2)*x.^(j-1);

            %f(x) =  Xest(1)*x.^5 + Xest(2)*x.^4 + Xest(3)*x.^3  + Xest(4)*x.^2 + Xest(5)*x.^1 + Xest(6) ;
            %construct polynomial


           % f(x) = f(x) +  a1(j)*cos(j*x) + b1(j)*sin(j*x);
                               
               
           f(x) = f(x-1)*exp(-beta*dT) + sigma*sqrt(1-exp(-2*beta*dT))*rn(x);
        %end
        
        
    end

%     for x = 1:sizey
%         f(x) = a0/2;
%         for j = 1:N
%             f(x) = f(x) +  a1(j)*cos(j*x) + b1(j)*sin(j*x);
%             %f(x) =  Xest(1)*x.^5 + Xest(2)*x.^4 + Xest(3)*x.^3  + Xest(4)*x.^2 + Xest(5)*x.^1 + Xest(6) ;
%             %construct polynomial
%         end
%     end

end

i
CurveFit = f;
BestFit = BestFitTemp



