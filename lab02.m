clear all;										
close all;										
msg=round(rand(1,1000));										
%1/2 rated convolutional Encoder										
trellis=poly2trellis(3,[6 7]);										
user=convenc(msg,trellis);										
% Convolutionally encoded data(0,1) are  mapping into +1/1										
%% To convert the binary sequences to bipolar NRZ format										
length_user=length(user);										
for i=1:length_user										
if user(i)==0										
user(i)=-1;										
end										
end										
fc=5000; %%carrier frequency, %KHz										
eb=.5;     %% energy per bit										
bitrate=1000;% 1KHz										
tb=1/bitrate; %% time per bit of message sequence										
chiprate=10000;										
tc=1/chiprate;										
%%% CDMA transmitter for a single user										
t=tc:tc:tb*length_user; 										
%%plotting base band signal for user										
basebandsig=[];										
for i=1:length_user										
for j=tc:tc:tb										
if user(i)==1 										
basebandsig=[basebandsig 1];										
else 										
basebandsig=[basebandsig -1];										
end										
end										
end										
figure(1)										
stairs(t(1:800),basebandsig(1:800))										
xlabel('Time(sec)')										
ylabel('Binary value')										
set(gca,'ytick',[ -1  1 ])										
title('A segment of original binary sequence for a single user')										
%%%% BPSK Modulation 										
bpskmod=[];										
for i=1:length_user										
for j=tc:tc:tb										
bpskmod=[bpskmod sqrt(2*eb)*user(i)*cos(2*pi*fc*j)];										
end										
end										
%length(bpskmod)										
number=length(t); %Total number of time segments										
spectrum=abs(fft(bpskmod));										
sampling_frequency=2*fc;										
sampling_interval=(1.0/sampling_frequency);										
nyquest_frequency=1.0/(2.0*sampling_interval);										
for i=1:number										
frequency(i)=(1.0/(number*sampling_interval)).*i;										
end										
figure(2)										
plot(frequency,spectrum)										
title('Frequency Domain analysis of BPSK modulated signal for a single user')										
xlabel('Frequency (Hz)')										
ylabel('Magnitude')										
grid on										
%% PN generator for  a single user										
%% let initial seed for a single user  is 1000										
seed=[1 -1 1 -1];  %convert it into bipolar NRZ format 										
spreadspectrum=[];										
pn=[];										
for i=1:length_user										
for j=1:10 %chip rate is 10 times the bit rate										
pn=[pn seed(4)];  										
if seed (4)==seed(3) temp=-1;										
else temp=1;										
end										
seed(4)=seed(3);										
seed(3)=seed(2);										
seed(2)=seed(1);										
seed(1)=temp;										
end										
end										
% each bit has 100 samples. and each pn chip has 10 samples. there r 										
% 10 chip per bit there fore size of pn samples and original bit is same										
pnupsampled=[];										
len_pn=length(pn);										
for i=1:len_pn										
for j=10*tc:10*tc:tb										
if pn(i)==1 										
pnupsampled=[pnupsampled 1];										
else 										
pnupsampled=[pnupsampled -1];										
end										
end										
end										
length_pnupsampled=length(pnupsampled);										
sigtx=bpskmod.*pnupsampled;										
figure(3)										
plot(t(1:200), sigtx(1:200))										
title('A segment of Transmitted  DS CDMA signal')										
xlabel('Time(sec)')										
ylabel('Amplitude')										
grid on										
%%%%%%%%%%%%%Adding fadig channel effect%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%										
chan=rayleighchan(1/chiprate,100);										
chan.ResetBeforeFiltering=0;										
fad=abs(filter(chan,ones(size(sigtx))));										
fadedsig=fad.*sigtx;										
snr_in_dBs=0:1.0:10;										
for m=1:length(snr_in_dBs)										
ber(m)=0.0;										
composite_signal=awgn(fadedsig,snr_in_dBs(m),'measured');  %% SNR of % dbs										
%%%%%%%%%%%%%DEMODULATION FOR USER 1%%%%%%%%%%%%%%%%%%%%										
rx=composite_signal.*pnupsampled;										
%%%% BPSK demodulation for a single user										
demodcar=[];										
for i=1:length_user										
for j=tc:tc:tb										
demodcar=[demodcar sqrt(2*eb)*cos(2*pi*fc*j)];										
end										
end										
bpskdemod=rx.*demodcar;										
len_dmod=length(bpskdemod);										
sum=zeros(1,len_dmod/10);										
for i=1:len_dmod/10										
for j=(i-1)*10+1:i*10										
sum(i)=sum(i)+bpskdemod(j);										
end										
end										
sum;										
rxbits=[];										
for i=1:length_user										
if sum(i)>0										
rxbits=[rxbits 1];										
else										
rxbits=[rxbits 0];										
end										
end										
tblen = 3; delay = tblen; % Traceback length										
decoded = vitdec(rxbits,trellis,tblen,'cont','hard');										
[number,rat] = biterr(decoded(delay+1:end),msg(1:end-delay));										
ber(m)=rat;										
end % for m										
figure(4)										
plot(snr_in_dBs,ber);										
xlabel('Signal to noise ratio(dB)');  										
ylabel('BER');										
legend('BER simulation for a single user');										
title(' Coded BER simulation under AWGN and Rayleigh fading channel ')										
grid on 										
										
