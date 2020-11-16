%plot data quick from trics folder in SciNet on Anna
clear;clc; %close all

h=1;
%mac:
filename = '/Volumes/scinet_iqoqi/measurements/trics_data/2017-11-02/173059/PMT1_1.txt';%long spectra
PMT{h}= importfile(filename);
%
h=2;
%mac:
filename = '/Volumes/scinet_iqoqi/measurements/trics_data/2017-11-02/164530/PMT1_1.txt';%long spectra
PMT{h}= importfile(filename);
%

save PMT;

load PMT;

for h=1:length(PMT);
%figure(2);
%subplot(2,1,1);
%plot(PMT(:,3)-PMT(7,3), PMT(:,5));
%plot(PMT{h}(:,3), PMT{h}(:,5));%plots the absolute double-pass AOM frequency vs mean excitation
%xlabel('AOM Frequency, MHz')
%hold on;
%for gg=0:10;
%plot(297.635*ones(10,1)+gg*0.700,0:0.1:0.9,'k');
%end
%subplot(2,1,2);
figure(1);

for gg=-3:2;
plot(0*ones(10,1)+gg*5.593,0:0.1:0.9,'--r','LineWidth',1);
plot(0*ones(10,1)+gg*5.593,0:0.1:0.9,'--r','LineWidth',1);
end

plot(2*(PMT{h}(:,3))-2*297.6325, PMT{h}(:,5),'-s'); %plots the real frequency difference seen by the ion
hold on
xlabel('Real Frequency, MHz')
%title('5000ms pulse, pi-time on carrier ~50 mus')
for gg=[-1,1];
plot(0*ones(10,1)+gg*1.155,0:0.1:0.9,'k');
plot(0*ones(10,1)+gg*1.706,0:0.1:0.9,'r');
end

for gg=[-1,1];
plot(-11.18*ones(10,1)+gg*1.155,0:0.1:0.9,'g');
plot(-11.18*ones(10,1)+gg*1.706,0:0.1:0.9,'b');
end

for gg=[-1,1];
plot(5.593*ones(10,1)+gg*1.155,0:0.1:0.9,'c');
plot(5.593*ones(10,1)+gg*1.706,0:0.1:0.9,'m');
end

for gg=[-1,1];
plot(-3*5.593*ones(10,1)+gg*1.155,0:0.1:0.9,'c');
plot(-3*5.593*ones(10,1)+gg*1.706,0:0.1:0.9,'m');
end

end

%subplot(2,1,1);
grid on; grid minor; %ylim([0 0.5]); xlim([-4 4])

%subplot(2,1,2); grid on; grid minor; ylim([0 0.5])