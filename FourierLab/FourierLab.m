%% PHYS 423 Fourier Lab Jason Pruitt

%% Ex 2 1-D Double slit
clear all; close all; clc;

sze = 5000;
slts = zeros(1, sze);
a = 100; % slit sep in microns
b = 20; % slit width in microns

xprime = linspace(-.5, .5, sze);
z = 2; % distance from slits to screen in meters
lam = 500e-9; % meters

k = (2*pi)/lam;
sinTheta = xprime./z; 
beta = (1/2)*k*b*(10^-6)*sinTheta; 
alpha = (1/2)*k*a*(10^-6)*sinTheta;

I1 = ((sin(beta)./beta).^2).*(cos(alpha)).^2; % Analytical expression
                                              % for double slit irradiance

slts(sze/2-a/2-b/2+1:sze/2-a/2+b/2) = ones(1,b); % left slit
slts(sze/2+a/2-b/2+1:sze/2+a/2+b/2) = ones(1,b); % right slit

f = fft(slts); 
f = fftshift(f);
I = f.*conj(f); 
I = I/max(I);

% plots superimposed 
fig1 = figure(1);
plot(xprime, I1, 'k')
hold on 
plot(xprime, I, 'b.')
xlabel('Screen Position [${m}$]', 'interpreter', 'latex')
ylabel('Irradiance [arb]')
legend('Analytical', 'FFT')
saveas(fig1, '1dDouble.png')

%% Ex Single Slit Diffraction with phase shift

lam = 500e-9; 
k = 2*pi/lam;
d = 50e-6; % slit width in meters
b = 50; % slith width in microns
sze = 5000; 
theta = linspace(-.35, .35, sze);

% analytical routine
s = sin(theta);
p = (sin(k*s*d/2)./(k*s*d/2) - sin(k*s*d/4)./(k*s*d/4)); % carried out in 
                                                         % meters
p = p.*p;
m = max(p);

% fourier transform routine
slt = zeros(1, sze);
slt(sze/2 - b/2 + 1:sze/2 + b/2 ) = ones(1, b); % slit
slt(sze/2 - b/4 + 1:sze/2 + b/4 ) = -1*ones(1, b/2); % phase element

% to visualize 
% figure(20)
% plot(slt)

f = fft(slt); 
f = fftshift(f);
I = f.*conj(f); 
I = I/max(I);

% Plots superimposed 
fig2 = figure(2);
plot(theta, p/m);
xlabel('angle [rad]');
ylabel('Irradiance [arb]')
hold on
plot(theta, I)
legend('Analytical', 'FFT')
saveas(fig2, '1dSingle.png')

%% Ex 3 2-D diffraction 
clear all; close all; clc;

sze = 1024;
slot_sze_x = 20;
slot_sze_y = 10;
aper = zeros(sze,sze);
slot = ones(slot_sze_x, slot_sze_y);
aper(sze/2-slot_sze_x/2:sze/2+slot_sze_x/2-1, ...
    sze/2-slot_sze_y/2:sze/2+slot_sze_y/2-1) = slot;
figure(35)
imagesc(aper);

f = fft2(aper);
lp = f.*conj(f);
lp = fftshift(lp);

fig3 = figure(3);
surf(log(1 + lp));
shading interp; 
view(125, 90)
gr = zeros(256,3);
gr(:,1) = (0:255)/(255*2);
gr(:,2) = (0:255)/(255*1.0);
gr(:,3) = (0:255)/(255*5);
colormap(gr)
xlabel('x')
ylabel('y')
zlabel('log(Irradiance)')
saveas(fig3, '2dDouble.png')
colorbar
