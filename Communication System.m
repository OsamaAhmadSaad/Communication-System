clear;
clc;

%Reading the first message
[m1,f1]=audioread('Short_RussianVoice.wav');
m1=m1(:,1)+m1(:,2);

%Reading the second message
[m2,f2]=audioread('Short_BBCArabic2.wav');
m2=m2(:,1)+m2(:,2);

%to make m1 of the same size of m2
z=zeros(37184,1); 
m1=[m1; z];

f1=10*f1;
f2=10*f2;
CH1=interp(m1,10);
CH2=interp(m2,10);

%Generating the carriers
fc1=100*10^3;
fc2=150*10^3;
t=0:(length(CH1)-1);
Carrier1=transpose(cos(2*pi*fc1*t*(1/f1)));
Carrier2=transpose(cos(2*pi*fc2*t*(1/f2)));

%The transmitted signals
S1=2*CH1.*Carrier1; 
S2=2*CH2.*Carrier2;
S_total=S1+S2;
n=length(S_total);

%Recieving process
fprintf("Select your Desired Channel\nFor First Channel press 1 \nand for Second channel press 2\n");
x=input('');
while ~(x==1||x==2) 
disp('Error')
x = input(' ');
end

if (x==1)
w=fc1; %Carrier Frequency 1
elseif(x==2)
w=fc2; %Carrier Frequency 2
end

%RF Stage (BPF)

BW=22*10^3;
RF_fpass1=w-BW; %Lower Cut off Frequency 
RF_fpass2=w+BW; %Upper Cut off Frequency but here we must insure that w+BW is less than 2WIF
RF_Lower_Coefficient= RF_fpass1/(f1/2);
RF_Upper__Coefficient= RF_fpass2/(f1/2);
[b,a]=butter(5,[RF_Lower_Coefficient 
RF_Upper__Coefficient]);
RF_out=filter(b,a,S_total);


%oscillator
offset=0;
WIF=25*10^3;
IF_carrier=transpose(cos(2*pi*(w+WIF+offset)*t*(1/f1)));
Mixer1_out=2.*RF_out.*IF_carrier;

%IF Stage (BPF)
IF_fpass1=WIF-BW; %Lower Cut off Frequency
IF_fpass2=WIF+BW; %Upper Cut off Frequency
IF_Lower_Coefficient= IF_fpass1/(f1/2);
IF_Upper__Coefficient= IF_fpass2/(f1/2);
[bb,aa]=butter(6,[IF_Lower_Coefficient IF_Upper__Coefficient]);
IF_out=filter(bb,aa,Mixer1_out);

%Getting the signal at baseband
Baseband_carrier=transpose(cos(2*pi*(WIF)*t*(1/f1)));
Mixer2_out=2.*Baseband_carrier.*IF_out;

%LPF 
Cut_off_freq=BW; %Cut off Frequency
Coffiecent=Cut_off_freq/(f1/2);
[b1, a1] = butter(6,Coffiecent);
LPF_out = 2*filter(b1,a1,Mixer2_out); 

%The Recieved Message
audiowrite('Recieved Message_test.wav',LPF_out,f1);

%Plotting
f=-length(S_total)/2:length(S_total)/2-1;

m1_f=fft(CH1,length(CH1)); %fourier transform of the first message
m2_f=fft(CH2,length(CH2)); %fourier transform of the second message
S_total_f=fft(S_total,n); %fourier transform of the transmitted signal
RF_out_f=fft(RF_out,length(RF_out)); %fourier transform of the signal after RF Stage
Mixer1_out_f=fft(Mixer1_out,length(Mixer1_out));%fourier transform of the signal after first mixer
IF_out_f=fft(IF_out,length(IF_out)); %fourier transform of the signal after the BPF
Mixer2_out_f=fft(Mixer2_out,length(Mixer2_out)); %fourier transform of the signal after second mixer
LPF_out_f=fft(LPF_out,n);

%plot(f*f1/n,fftshift(abs(m1_f)));
%plot(f*f1/n,fftshift(abs(m2_f)));
%plot(f*f1/n,fftshift(abs(S_total_f)));
%plot(f*f1/length(RF_out_f),fftshift(abs(RF_out_f)));
%plot(f*f1/length(Mixer1_out_f),fftshift(abs(Mixer1_out_f)));
%plot(f*f1/length(IF_out_f),fftshift(abs(IF_out_f)));
%plot(f*f1/length(Mixer2_out_f),fftshift(abs(Mixer2_out_f)));
plot((f*f1/length(LPF_out_f)),fftshift(abs(LPF_out_f)));
