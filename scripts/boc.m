


f_L1 = 1575.42 * 1e6; % GPS L1 in MHz

f_BOC_Carrier = 1 * 1e6;
f_BOC_Data = 1 * 1e6;



t = [0:0.01/f_L1:10/f_BOC_Carrier];
N = length(t);

C = sin (2*pi*f_L1*t);


t1 = [0:0.01/f_BOC_Carrier:10/f_BOC_Carrier];
BOC_C = square(2*pi*f_BOC_Carrier*t);

BOC_D = ones(1,length(BOC_C));

% modulate carrier with BOC

C1 = conv(C,BOC_C);

f1 = fft(C1);

%plot(t,BOC_C);
%plot(t,C);


