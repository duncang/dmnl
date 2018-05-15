function [CurveFit, BestFit,  normX, Xest] = LeastSquaresBestFitWindFourierNOrder(Signal, Order, FindBestFit);


%Least Squares Data Fitting Fourier Series

%Fourier series of the form f(x) = a0/2 + a1*cos(x) ..+aN*cos(Nx) + b1*sin(x)..+
%bN*sin(Nx)



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


x = 1:sizey;

%order of polynomial to fit to data
N = BestFitTemp;


a0 = 0;

a1 = zeros(1,N);  %first estimate is zero

b1 = zeros(1,N);
%Vector of initial guess of parameters a0 a1 etc..

%number of states will be N+1 , the +1 is because of the a0/2 term, since i
%estimate a0 as well.


Xest = [a0, a1, b1];   %state vector

%preallocate arrays

Minit = zeros(sizey,2*N+1);
%Xest = zeros(1,N+1);
f = zeros(1,sizey);



for i = 1:1000


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

    for x = 1:sizey
        f(x) = a0/2;
        for j = 1:N
            f(x) = f(x) +  a1(j)*cos(j*x) + b1(j)*sin(j*x);
            %f(x) =  Xest(1)*x.^5 + Xest(2)*x.^4 + Xest(3)*x.^3  + Xest(4)*x.^2 + Xest(5)*x.^1 + Xest(6) ;
            %construct polynomial
        end
    end

end

i
CurveFit = f;
BestFit = BestFitTemp



