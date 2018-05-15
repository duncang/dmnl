function [CurveFit, BestFit] = LeastSquaresBestFit(Signal, Order, FindBestFit);


%Least Squares Data Fitting

%by Troy Bruggemann 22 Sept 2005 (c)



%Give a polynomial order to fit or have it detect which is the best fit order polynomial

% 1. enter initial guess for fitting parameters
% 2.  Test the initial guess against the data
% 3.  Use least squares to refine the initial guess until the least squares is minised
% 4 . The fit curve is then the best fit of the data.
% 5. Repeats the process for each order and finds the best ordered polynomial that best fits the data


%~Input~:
%Signal: The input signal which is to have a polynomial fitted to it using least squares.
%Order: The order of the fitted polynomial to use
%FindBestFit: if = 1, user wishes the function to find the best fit curve to the signal, if 0 then the function will use the Order as the order of the polynomial to fit


%~Output~:
%CurveFit: The polynomial of order 'BestFit' to the input signal
%BestFit: The order of the fitted polynomial


%Initial guess for fitting parameters is a vector of 0's.

%go from order 1 to order 20 and find the best fit to the data by using the order which has the best fit to the signal


BestFitTemp = Order; 
clear y Xest f M Resid DeltaX normx;

warning off;%turn off annoying warnings for singular inverse matrices


%find the best fit curve


if FindBestFit == 1    
    
    for N = 1:40
        
        clear y Xest f M Resid DeltaX normx;        
        
        y = Signal;        
        numpoints = size(y);
        sizey = numpoints(2);        
        x = 1:sizey;        
        %order of polynomial to fit to data
        %N = 6;
        
        %Vector of initial guess of parameters a0 a1 etc..
        
        Xest = zeros(1,N+1);        
        %preallocate arrays
        
        M = zeros(sizey,N+1);
        Xest = zeros(1,N+1);
        f = zeros(1,sizey);
        
        
        
        for i = 1:50                    
            for x = 1:sizey                
                f(x) = Xest(1)*x.^N;                
                for j = 1:N
                    
                    f(x) = f(x) + Xest(N-j+2)*x.^(j-1);                    
                    
                    %f(x) =  Xest(1)*x.^5 + Xest(2)*x.^4 + Xest(3)*x.^3  + Xest(4)*x.^2 + Xest(5)*x.^1 + Xest(6) ;
                    %construct polynomial  
                    
                end                
            end
            
            
            
            
            Resid = y-f;
            
            
            
            for x = 1:sizey  
                
                for j = 1:N+1
                    
                    M(x,j) =  x^(N-j+1);
                    
                    %         M(x,1) = x^5;
                    %         M(x,2) =  x^4;  
                    %         M(x,3) =  x^3;
                    %         M(x,4) =  x^2;
                    %         M(x,5) =  x;
                    %         M(x,6) =  1;
                    
                end            
                
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
            f(x) = Xest(1)*x.^N;            
            for j = 1:N
                
                f(x) = f(x) + Xest(N-j+2)*x.^(j-1);                
                
                %f(x) =  Xest(1)*x.^5 + Xest(2)*x.^4 + Xest(3)*x.^3  + Xest(4)*x.^2 + Xest(5)*x.^1 + Xest(6) ;
                %construct polynomial  
                
            end            
        end               
        Order(N) = sum((y-f).^2);         
    end    
    
    %evaluate choice of order 
    
    Minimum = min(Order);
    BestFitTemp = find(Order == Minimum);
    BestFitTemp = BestFitTemp(1);  %in case it returns two or more orders at the same minimum, choose the first one
    
end

clear y Xest f M Resid DeltaX normx;

y = Signal;

numpoints = size(y);
sizey = numpoints(2);
x = 1:sizey;

%order of polynomial to fit to data
N = BestFitTemp;

%Vector of initial guess of parameters a0 a1 etc..

Xest = zeros(1,N+1);

%preallocate arrays

M = zeros(sizey,N+1);
Xest = zeros(1,N+1);
f = zeros(1,sizey);


for i = 1:50       
    for x = 1:sizey        
        f(x) = Xest(1)*x.^N;        
        for j = 1:N            
            f(x) = f(x) + Xest(N-j+2)*x.^(j-1);
            
            %f(x) =  Xest(1)*x.^5 + Xest(2)*x.^4 + Xest(3)*x.^3  + Xest(4)*x.^2 + Xest(5)*x.^1 + Xest(6) ;
            %construct polynomial  
            
        end
        
        
    end
        
    Resid = y-f;
    
    for x = 1:sizey          
        for j = 1:N+1            
            M(x,j) =  x^(N-j+1);                 
        end        
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
    
    for x = 1:sizey        
        f(x) = Xest(1)*x.^N;        
        for j = 1:N            
            f(x) = f(x) + Xest(N-j+2)*x.^(j-1);                        
            %f(x) =  Xest(1)*x.^5 + Xest(2)*x.^4 + Xest(3)*x.^3  + Xest(4)*x.^2 + Xest(5)*x.^1 + Xest(6) ;
            %construct polynomial              
        end        
    end       
     
end

 i
 CurveFit = f;
 BestFit = BestFitTemp  

 

