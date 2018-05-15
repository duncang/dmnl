% =========================================================================
% M-FILE NAME:  
% VERSION:      $Version$
% AUTHOR:       Lennon Cork
% LAST UPDATE:  $Id: FilterEquations.m 785 2007-03-06 07:59:15Z greerd $
% HISTORY:      
%               $Log$
%               Revision 1.1  2007/01/09 12:27:26  greerd
%               Added lennons UKF equations
% 
% ABSTRACT:     
% REFERENCES:   
% NOTATIONS:    
% =========================================================================
% State Data: xdata = [initial ; enable ; max ; min ; limit mode];
% Output Data: xdata = [initial ; enable ; max ; min ; limit mode];

function [x,P,z,S,e,K,eError] = FilterEquations(Mx,Mz,xdata,zdata,u,y,v,w,ts,P0,Q,R,form,ftype);
if (ftype == 1)
     [x,P,z,S,e,K,eError] = FullStateUKF(Mx,Mz,xdata,zdata,u,y,v,w,ts,P0,Q,R,form);
else
     error('Unrecognised Filter Type');
     x = zeros(length(xdata),1); P = zeros(length(xdata)); z = zeros(length(zdata),1);
     S = zeros(length(zdata)); e = zeros(length(zdata),1); K = zeros(length(xdata),length(zdata));
     eError = 1;
end

% -------------------------------------------------------------------------
% FULLSTATEUKF: Full State Uncented Kalman Filter (UKF)
% -------------------------------------------------------------------------
function [x,P,z,S,e,K,eError] = FullStateUKF(Mx,Mz,xdata,zdata,u,y,v,w,ts,P0,Q,R,form);
% Variable Initalisation
eError = 0;
x0 = xdata(:,1); se = xdata(:,2); z0 = zdata(:,1); oe = zdata(:,2); 
nx = length(x0); nz = length(z0); nv = length(v); nw = length(w); 
fs = find(se); fo = find(oe); na = nx+nv+nw;
X = zeros(nx,2*na+1); Z = zeros(nz,2*na+1);
% Augmented State (zero mean)
xa = [x0' zeros(nv,1)' zeros(nw,1)']';
% Augmented Covariance
Pa = [P0            zeros(nx,nv)    zeros(nx,nw); ...
      zeros(nv,nx)  Q               zeros(nv,nw); ...
      zeros(nw,nx)  zeros(nw,nv)    R           ];
% Augmented Sigma Points (Start With k = 100)
[Xa G] = nut(Pa,xa,1);
% Check if Data is Imaginary
i = 1;  % Set i - counter to stop infinite loop
while (~isreal(Xa)||~all(all((isfinite(Xa))))) && (i <= 5)
   % Generate Augmented Sigma Points with new k
   [Xa G] = nut(Pa,xa,10^i); 
   i = i + 1;   % Incriment i
   if (i == 5); 
       warning('error calculating sigma points'); 
       eError = 1;
   end
end
% SIGMA POINT SEPERATION - Remove reference to data being available
% Seperated Sigma Points
if (eError == 1)
    X0 = x0*ones(1,2*na+1); V = zeros(nv,2*na+1); W = zeros(nw,2*na+1);
else
    X0 = Xa(1:nx,:); V = Xa(nx+1:nx+nv,:);  W = Xa(nx+nv+1:nx+nv+nw,:);
end
% Sigma Point Calculations
for i = 1:(2*na+1)
    % State Equations and Output
    [X(:,i) Z(:,i) dError] = SystemEquations(Mx,Mz,[X0(:,i), xdata(:,2:5)],zdata,u,v+V(:,i),w+W(:,i),ts,form);
    if (dError == 1); eErro = 1; break; end
    % Zero limit for correct calculations
    X(:,i) = limit(X(:,i),-(xdata(:,3)-xdata(:,4))/2,(xdata(:,3)-xdata(:,4))/2,xdata(:,5));
    Z(:,i) = limit(Z(:,i),-(zdata(:,3)-zdata(:,4))/2,(zdata(:,3)-zdata(:,4))/2,zdata(:,5));
end
if (~isreal(X)||~all(all((isfinite(X))))); warning('se error'); eError = 1; end;
% Initalisation
x = x0; z = z0; P = P0;
e = zeros(nz,1); S = zeros(nz,nz); PS = zeros(nx,nz); K = zeros(nx,nz);
% State Mean, Output Mean and Output Error
x(fs) = umean(X(fs,:),G); z(fo) = umean(Z(fo,:),G); e(fo) = y(fo)-z(fo);
% State, Output and State-Output Covariance
P(fs,fs) = ucov(X(fs,:),x(fs),X(fs,:),x(fs),G); 
S = ucov(Z,z,Z,z,G); 
PS(fs,fo) = ucov(X(fs,:),x(fs),Z(fo,:),z(fo),G);
% Kalman Gain
if (eError == 0)
    if (cond(S(fo,fo))>1e25)
        warning('conditioning error'); 
        eError = 1;
    else
        K(fs,fo) = PS(fs,fo)*pinv(S(fo,fo));
    end
end;
x(fs) = x(fs) + K(fs,fo)*e(fo); P(fs,fs) = P(fs,fs) - K(fs,fo)*S(fo,fo)*K(fs,fo)';
% Relimit for consistancy
x = limit(x, xdata(:,4), xdata(:,3), xdata(:,5));
z = limit(z, zdata(:,4), zdata(:,3), zdata(:,5));





