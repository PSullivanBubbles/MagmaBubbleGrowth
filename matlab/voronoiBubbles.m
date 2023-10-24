clear;

%%%  Script completing part 2 of the numerical modelling exercise



volumes = 10*[0.1;0.15; 0.2; 1];

%% Run shell model for each point in the array
radii = {};
time={};

parfor l=1:size(volumes,1)
    disp("Shell model "+ num2str(l)+" of " + num2str(size(volumes,1)))
    [radii(l), time(l)]=Bubble_Growth_Modelling_Function(volumes(l));
end

t= 0:1:(200*3600);
phi = zeros(size(t,2),1);
phiT = zeros(size(t,2),size(volumes,1));

for i=1:size(t,2)

    tau=t(i);
    liquidVolume =0;
    gasVolume=0;
    for j=1:size(volumes,1)
        R = findR(radii{j},time{j},tau);

        liquidVolume = liquidVolume+ volumes(j)/(1e11);
        gasVolume = gasVolume + (4/3).*pi()*R.^3;

    end
    phi(i) = gasVolume/(gasVolume+liquidVolume);

end


R=zeros(size(volumes,1),size(t,2));


for j=1:size(volumes,1)
    R(j,:) = spline(time{j},radii{j},t);
end


figure();
hold on;
for i=1:size(volumes,1)
    semilogx(t./(volumes(i).^(2/3)),R(i,:)./(volumes(i).^(1/3)));
end
axis([0 1e4 0 2e-4])
hold off;

DP= 1;
mu=1e2;
R0=1e-5;

Vexp = DP*R0/4/mu * exp(DP*t/4/mu);

Vscale = 1.9e-6;
Vsqrt = Vscale./(t.^0.5);

V = diff(R,1,2);
figure(8);
%v =[0; 0; 0; 0; 0;0;0];
v = zeros(size(volumes));
V=[v V];

semilogx (t/3600/volumes(1)^(2/3),V(1,:)./Vsqrt);
hold on;
for i=2:size(V,1)
    semilogx (t/3600/volumes(i)^(2/3),V(i,:)./Vsqrt);
end
semilogx(t/3600,Vexp./Vsqrt);
semilogx(t/3600,Vsqrt./Vsqrt);
legend( 'Bubble 1', 'Bubble 2', 'Bubble 3', 'Bubble 4','initial fit', 'late fit');
axis([1e-4 50 0 1.2])
hold off;





liquidVolume1 =volumes(1)/(1e11);
liquidVolume2 =volumes(2)/(1e11);
liquidVolume3 =volumes(3)/(1e11);
liquidVolume4 =volumes(4)/(1e11);
% liquidVolume5 =volumes(5)/(1e11);
% liquidVolume6 =volumes(6)/(1e11);
% liquidVolume7 =volumes(7)/(1e11);

phi1 = zeros(size(t,2),1);
phi2=phi1;
phi3=phi1;
phi4=phi1;
% phi5=phi1;
% phi6=phi1;
% phi7=phi1;

for i=1:size(t,2)
    gasVolume1=(4/3).*pi()*R(1,i).^3;
    gasVolume2=(4/3).*pi()*R(2,i).^3;
    gasVolume3=(4/3).*pi()*R(3,i).^3;
    gasVolume4=(4/3).*pi()*R(4,i).^3;
    %     gasVolume5=(4/3).*pi()*R(5,i).^3;
    %     gasVolume6=(4/3).*pi()*R(6,i).^3;
    %     gasVolume7=(4/3).*pi()*R(7,i).^3;

    phi1(i) = gasVolume1/(gasVolume1+liquidVolume1);
    phi2(i) = gasVolume2/(gasVolume2+liquidVolume2);
    phi3(i) = gasVolume3/(gasVolume3+liquidVolume3);
    phi4(i) = gasVolume4/(gasVolume4+liquidVolume4);
    %     phi5(i) = gasVolume5/(gasVolume5+liquidVolume5);
    %     phi6(i) = gasVolume6/(gasVolume6+liquidVolume6);
    %     phi7(i) = gasVolume7/(gasVolume7+liquidVolume7);

end
figure(6);
plot(t/3600,phi,'--');
hold on;
plot(t/3600,phi1);
plot(t/3600,phi2);
plot(t/3600,phi3);
plot(t/3600,phi4);
% plot(t/3600,phi5);
% plot(t/3600,phi6);
% plot(t/3600,phi7);
axis([0 10 0 1]);
legend('Global Porosity', 'Porosity Bubble 1', 'Porosity Bubble 2', 'Porosity Bubble 3', 'Porosity Bubble 4', 'Porosity Bubble 5', 'Porosity Bubble 6', 'Porosity Bubble 7');
hold off;


%% FFT of porosities

Fs = 1/t(2);            % Sampling frequency
T = 1/Fs;             % Sampling period
L=size(t,2);          % Length of signal
f = Fs*(0:(L/2))/L;

Y=fft(phi);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);

Y1=fft(phi1);
P21 = abs(Y1/L);
P11 = P21(1:L/2+1);
P11(2:end-1) = 2*P11(2:end-1);

Y2=fft(phi2);
P22 = abs(Y2/L);
P12 = P22(1:L/2+1);
P12(2:end-1) = 2*P12(2:end-1);

Y3=fft(phi3);
P23 = abs(Y3/L);
P13 = P23(1:L/2+1);
P13(2:end-1) = 2*P13(2:end-1);

Y4=fft(phi4);
P24 = abs(Y4/L);
P14 = P24(1:L/2+1);
P14(2:end-1) = 2*P14(2:end-1);
%
% Y5=fft(phi5);
% P25 = abs(Y5/L);
% P15 = P25(1:L/2+1);
% P15(2:end-1) = 2*P15(2:end-1);
%
% Y6=fft(phi6);
% P26 = abs(Y6/L);
% P16 = P26(1:L/2+1);
% P16(2:end-1) = 2*P16(2:end-1);
%
% Y7=fft(phi7);
% P27 = abs(Y7/L);
% P17 = P27(1:L/2+1);
% P17(2:end-1) = 2*P17(2:end-1);

%Pmean = mean([phi1,phi2,phi3,phi4,phi5,phi6,phi7],2);
Ymean = volumes(1)/sum(volumes)*Y1 + volumes(2)/sum(volumes)*Y2 + volumes(3)/sum(volumes)*Y3 + volumes(4)/sum(volumes)*Y4;
Pmean = volumes(1)/sum(volumes)*P11 + volumes(2)/sum(volumes)*P12 + volumes(3)/sum(volumes)*P13 + volumes(4)/sum(volumes)*P14;



figure(20);
hold off;
loglog(f,P1,'--');
hold on;
loglog(f,P11,'.');
loglog(f,P12,'.');
loglog(f,P13,'.');
loglog(f,P14,'.');
% loglog(f,P15,'.');
% loglog(f,P16,'.');
% loglog(f,P17,'.');
loglog(f,Pmean,'--');
%loglog(f,3e-7 * f.^-1);
legend('Global Porosity', 'Porosity Bubble 1', 'Porosity Bubble 2', 'Porosity Bubble 3', 'Porosity Bubble 4','Mean Porosity');

%legend('Global Porosity', 'Porosity Bubble 1', 'Porosity Bubble 2', 'Porosity Bubble 3', 'Porosity Bubble 4', 'Porosity Bubble 5', 'Porosity Bubble 6', 'Porosity Bubble 7','Mean Porosity');
hold off;



figure(21);
hold off;
loglog(f,P1,'--');
hold on;
loglog(f*(volumes(1)^(2/3)),P11/(volumes(1)^(2/3)),'.');
loglog(f*(volumes(2)^(2/3)),P12/(volumes(2)^(2/3)),'.');
loglog(f*(volumes(3)^(2/3)),P13/(volumes(3)^(2/3)),'.');
loglog(f*(volumes(4)^(2/3)),P14/(volumes(4)^(2/3)),'.');
% loglog(f,P15,'.');
% loglog(f,P16,'.');
% loglog(f,P17,'.');
loglog(f,Pmean,'--');
%loglog(f,3e-7 * f.^-1);
legend('Global Porosity', 'Porosity Bubble 1', 'Porosity Bubble 2', 'Porosity Bubble 3', 'Porosity Bubble 4','Mean Porosity');

%legend('Global Porosity', 'Porosity Bubble 1', 'Porosity Bubble 2', 'Porosity Bubble 3', 'Porosity Bubble 4', 'Porosity Bubble 5', 'Porosity Bubble 6', 'Porosity Bubble 7','Mean Porosity');
hold off;




weightedMean = ifft(Ymean);

figure(6);
plot(t/3600,phi,'--');
hold on;
plot(t/3600,phi1);
plot(t/3600,phi2);
plot(t/3600,phi3);
plot(t/3600,phi4);
plot(t/3600,weightedMean);

axis([0 10 0 1]);
legend('Global Porosity', 'Porosity Bubble 1', 'Porosity Bubble 2', 'Porosity Bubble 3', 'Porosity Bubble 4', 'Weighted Mean');
hold off;
