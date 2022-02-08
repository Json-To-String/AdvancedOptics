%% Modeling Optical Waveguides Lab - Jason Pruitt
clear all; close all; clc;
format long

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% n_eff vs d graph for waveguide TE modes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lam = .633e-6;  % wavelength
k = 2*pi/lam;   % wavenum
n1 = 1.0;       % cover/air
n2 = 1.5095;    % guide
n3 = 1.4711;    % substrate
thick = 0:.01e-6:10e-6;
thickM = 0:.01:10;

figure()
modes = linspace(0, 9, 10);
nArrTE = zeros(1, length(thick));

for j = 1:length(modes)
    for i = 1:length(thick)
        [nArrTE(i), ~] = binarySplit(n1, n2, n3, ...
                                     thick(i), k, modes(j));
    end
    plot(thickM, nArrTE, '.')
    hold on
end
xlabel('Guide thickness [\mum]')
ylabel('Effective index')
title('TE')
ylim([1.4712 1.51])

figure()
nArrTM = zeros(1, length(thick));
for j = 1:length(modes)
    for i = 1:length(thick)
        [~, nArrTM(i)] = binarySplit(n1, n2, n3, ...
                                     thick(i), k, modes(j));
    end
    plot(thickM, nArrTM, '.')
    hold on
end
xlabel('Guide thickness [\mum]')
ylabel('Effective index')
title('TM')
ylim([1.4712 1.51])
xlim([0.29 10])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Electric field in a guide
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

d = 2e-6;       % thickness in m                   
mode0 = 0;
mode1 = 1;

[neffTE0, neffTM0] = binarySplit(n1, n2, n3, d, k, mode0);
[neffTE1, neffTM1] = binarySplit(n1, n2, n3, d, k, mode1);

[x0, E0] = eField(n1, n2, n3, neffTE0, d, k, mode0);
[x1, E1] = eField(n1, n2, n3, neffTE1, d, k, mode1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Extra credit
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

makeGif(E0, mode0)
makeGif(E1, mode1)


%% functions

function [midTE, midTM] = binarySplit(n1, n2, n3, d, k, mode)

    funcTE = @(x) 2*d*sqrt(k*k*(n2*n2 - x*x)) ... % extra path
                 - 2*atan(sqrt(x*x - n1*n1)./sqrt(n2*n2 - x*x)) ... 
                 - 2*atan(sqrt(x*x - n3*n3)./sqrt(n2*n2 - x*x)) ... 
                 - 2*mode*pi;
    
    funcTM = @(x) 2*d*sqrt(k*k*(n2*n2 - x*x)) ... % extra path
                 - 2*atan(sqrt(x*x - n1*n1)./(((n1/n2)^2)*sqrt(n2*n2 - x*x))) ... 
                 - 2*atan(sqrt(x*x - n3*n3)./(((n2/n3)^2)*sqrt(n2*n2 - x*x))) ... 
                 - 2*mode*pi;
    
    % TE Loop
    tol = 1e-5;
    nup = n2 - tol;
    ndo = n3 + tol;
    midTE = (nup + ndo)/2;
    while abs(midTE - nup) > tol
       prod = funcTE(nup)*funcTE(midTE);
       if prod > 0
           nup = midTE;
       else
           ndo = midTE;
       end
       midTE = (nup + ndo)/2;
    end
    
    % TM Loop
    tol = 1e-5;
    nup = n2 - tol;
    ndo = n3 + tol;
    midTM = (nup + ndo)/2;
    while abs(midTM - nup) > tol
       prod = funcTM(nup)*funcTM(midTM);
       if prod > 0
           nup = midTM;
       else
           ndo = midTM;
       end
       midTM = (nup + ndo)/2;
    end
end

function [x, E] = eField(n1, n2, n3, neff, d, k, mode)


    % for within the guide
    k2x = k*sqrt(n2^2 - neff^2);

    % for outside the guide
    a1x = k*sqrt(neff^2 - n1^2);
    a3x = k*sqrt(neff^2 - n3^2);

    phi = atan(a1x/k2x) - k2x*d/2;

    A1 = 1;
    A2 = exp(-a1x*d/2)/cos((k2x*d/2) + phi);
    A3 = A2*cos((-k2x*d/2) + phi)/exp(-a3x*d/2);

    x = linspace(-d, d, 100);
    E = zeros(1, length(x));

    for i = 1:length(x)

        if x(i) > d/2
            E(i) = A1*exp(-a1x*x(i));
    %         fprintf('\nAbove')

        elseif x(i) < d/2 && x(i) > 0
            E(i) = A2*cos(k2x*x(i) + phi);
    %         fprintf('\nWithin')

        elseif x(i) < 0 && x(i) > -d/2
            E(i) = A2*cos(k2x*x(i) + phi);
    %         fprintf('\nWithin')

        elseif x(i) < -d/2
            E(i) = A3*exp(a3x*x(i));
    %         fprintf('\nBelow')
        end
    end
    
    figure()
    plot(x, E)
    title("Electric field in Guide for TE mode " + num2str(mode))
    xline(d/2, 'k--')
    xline(-d/2, 'k--')
    xlabel('x [m]')
    ylabel('E [arb]')

end

function makeGif(E, mode)

    axis tight manual
    filename = "oscE" + num2str(mode) + ".gif";
    f = figure();
    for t = 1:100
        plot(E*cos(t/10))
        title('Oscillating electric field')
        xlabel('Position [m]')
        ylabel('E-field [arb]')
        ylim([-1.5e-4 1.5e-4])
        drawnow

        % Capture plot as image
        frame = getframe(f);
        im = frame2im(frame);
        [imind, cm] = rgb2ind(im, 256);

        % Write to gif file
        if t == 1
            imwrite(imind, cm, filename, 'gif', 'Loopcount', inf);
        else
            imwrite(imind, cm, filename, 'gif', 'WriteMode', 'append');
        end

    end
end