%%% Implementation of J. Farrell IMU coast: Not a Silver bullet, ION 1999
%% 
%% $Id: FarrellIMUCoast.m 1875 2008-07-15 04:42:43Z n2523710 $
%%
grav=9.8;
N = input('Enter # of holding pattern cycles   ');
M = 4*N;
tleg = 60.0;
delt = 1.0;
fps = (6076./3600)*input('Enter speed in Kts   ');
speed = fps*12./39.37;
Iq = round(tleg/delt);
I = M*Iq;
f1 = (1.e-5)*input('Enter ROLL Misalignment Multiplier (ref to 10ppm)  ');
f2 = (1.e-5)*input('Enter PITCH Misalignment Multiplier (ref to 10ppm)  ');
rerr = zeros(3,1);
verr = zeros(3,1);
Y = zeros(3,1);
vx = zeros(I,1);
vy = zeros(I,1);
rx = zeros(I,1);
ry = zeros(I,1);
qd = zeros(I,1);
hdg = zeros(I,1);
xe = zeros(I,1);
ye =  zeros(I,1);
vex =  zeros(I,1);
vey =  zeros(I,1);
Yx =  zeros(I,1);
Yy =  zeros(I,1);


for I=1:I,
    k = 1+rem(I-1,Iq);
    iq = 1+fix((I-1)/Iq);
    qd(I) = 1+rem(iq-1,4);
    j = 1+rem(I-1,Iq*4);
    ir = 1+fix((I-1)/(Iq*4));
    
    if(qd(I)==1)
        cR = 1./sqrt(1.+(speed*(pi/tleg)/grav)^2);
        sR = sqrt(1.-cR^2);
        hdg(I) = pi*k*delt/tleg;
        cH = cos(hdg(I));
        sH = sin(hdg(I));
        TAP = [1 0 0;0 cR sR; 0 -sR cR]*[cH sH 0;-sH cH 0;0 0 1];
        Y = Y + TAP' * (pi * delt/tleg) * [f1;f2;0]; 
        if(I <= Iq*4)
            vx(j) = speed*cH;
            vy(j) = speed*sH;
            rx(j) = (tleg * speed/pi)*sH;
            ry(j) = (tleg * speed/pi)*(1.-cH);
        end
    elseif(qd(I) == 2)
        cR=  1;
        sR = 0;
        hdg(I) = pi;
        cH = -1;
        sH = 0;
        TAP = [1 0 0;0 cR sR; 0 -sR cR]*[cH sH 0;-sH cH 0;0 0 1];
        
        if(I<=Iq*4)
            vx(j) = -speed;
            vy(j) = 0;
            rx(j) = -speed*k*delt;
            ry(j) = speed*2*tleg/pi;
        end
    elseif(qd(I)==3)
        cR = 1./sqrt(1.+(speed*(pi/tleg)/grav)^2);
        sR = sqrt(1-cR^2);
        hdg(I) = pi*(k*delt/tleg - 1);
        cH = cos(hdg(I));
        sH = sin(hdg(I));
        TAP = [1 0 0;0 cR sR; 0 -sR cR]*[cH sH 0;-sH cH 0;0 0 1];
        Y = Y + TAP' * (pi * delt/tleg) * [f1;f2;0]; 
        if(I <= Iq*4)
            vx(j) = speed*cH;
            vy(j) = speed*sH;
            rx(j) = (tleg * speed)*(sH/pi-1.);
            ry(j) = (tleg * speed/pi)*(1.-cH);
        end
    elseif(qd(I) ==4)
                cR=  1;
        sR = 0;
        hdg(I) = pi;
        cH = -1;
        sH = 0;
        TAP = [1 0 0;0 cR sR; 0 -sR cR];
        if(I<=Iq*4)
            vx(j) = speed;
            vy(j) = 0;
            rx(j) = speed*(k*delt-tleg);
            ry(j) = 0;
        end
    end
    Y = Y + (delt/6378137) * [-verr(2);verr(1);0];
    Yx(I) = Y(1);
    Yy(I) = Y(2);
    A = TAP' * [0;0;-grav/cR];
    YxA = [Y(2) * A(3) - Y(3) * A(2);Y(3)*A(1) - Y(1)*A(3);Y(1)*A(2) - Y(2)*A(1)];
    rerr = rerr + delt * verr/2;
    verr = verr + delt * YxA;
    vex(I) = verr(1);
    vey(I) = verr(2);
    rerr = rerr + delt * verr /2 ;
    x(ir,j) = rx(j) - rerr(1);
    y(ir,j) = ry(j) - rerr(2);
    xe(I) = rerr(1);
    ye(I) = rerr(2);
end

figure();
subplot(2,2,1),plot(Yx);
ylabel('Yx');
subplot(2,2,2),plot(Yy);
ylabel('Yy');
subplot(2,2,3),plot(vex);
ylabel('vex');
subplot(2,2,4),plot(vey);
ylabel('vey');
ry(Iq) = 10000; % ficticious point

if N==8, ry(Iq) = 12000;,end

figure();
axis('square');
subplot(111),plot(ry,rx,'.r',y(ir,:),x(ir,:),'.g');
grid off;
title(['RACE TRACK HOLDING PATTERN - ',int2str(N),' CYCLES']);
xlabel('East Excursion - Metres');
ylabel('North Excursion - Metres');

figure();
subplot(211),plot(xe)
grid on;
title(['Position Error - ',int2str(N),' CYCLES']);
ylabel('North Position Error - Metres');
subplot(212),plot(ye);
grid on
xlabel('Number of Seconds Since Start');
ylabel('East Position Error - Metres');

        