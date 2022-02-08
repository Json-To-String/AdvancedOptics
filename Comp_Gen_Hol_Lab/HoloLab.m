%% Holography Lab Jason Pruitt

clear; close all; clc;

rng('shuffle');
sze0 = 64;

hol = rand(sze0, sze0); % a set of numbers between 0 and 1
                        % change to a random array of -1 or 1
for i = 1:sze0
    for j = 1:sze0
        if hol(i,j)<0.5
            hol(i,j) = -1;
        else
            hol(i,j) = 1;
        end
    end
end

target1 = target4by4(sze0); % 4x4 16 pixel spot target
rect16 = targetRect(sze0);  % 16x16 rectangle pixel

[targetA, holf2A, holA, err_arrA] = annealing(sze0, hol, target1);
[targetB, holf2B, holB, err_arrB] = binSearch(sze0, hol, target1);
[targetC, holf2C, holC, err_arrC] = binSearch(sze0, hol, rect16);

figure(10)
subplot(2,2,1)
imagesc(targetB)
title('target pattern')

subplot(2,2,2)
imagesc(holf2B)
title('Hologram output')

subplot(2,2,3)
imagesc(hol);
title('Original hologram')

subplot(2,2,4)
plot(err_arrB);
title('error')
xlabel('Iterations')
ylabel('Error')
sgtitle('Binary Search')

figure(20)
subplot(2,2,1)
imagesc(targetA)
title('target pattern')

subplot(2,2,2)
imagesc(holf2A)
title('Hologram output')

subplot(2,2,3)
imagesc(hol);
title('Original hologram')

subplot(2,2,4)
plot(err_arrA);
title('error')
xlabel('Iterations')
ylabel('Error')
sgtitle('Simulated Annealing')

figure(30)
subplot(2,2,1)
imagesc(rect16)
title('target pattern')

subplot(2,2,2)
imagesc(holf2C)
title('Hologram output')

subplot(2,2,3)
imagesc(hol);
title('Original hologram')

subplot(2,2,4)
plot(err_arrC);
title('error')
xlabel('Iterations')
ylabel('Error')
sgtitle('Binary Search')

%% functions
function target = target4by4(sze)
    % 4x4 target spot pattern
    target = zeros(sze, sze);
    cent = sze/2+1; 

    target(cent + 8, cent + 8) = 1;
    target(cent + 8, cent + 24) = 1;
    target(cent + 24, cent + 24) = 1;
    target(cent + 24, cent + 8) = 1;

    target(cent - 8, cent - 8) = 1;
    target(cent - 8, cent - 24) = 1;
    target(cent - 24, cent - 24) = 1;
    target(cent - 24, cent - 8) = 1;

    target(cent + 8, cent - 8) = 1;
    target(cent + 8, cent - 24) = 1;
    target(cent + 24, cent - 24) = 1;
    target(cent + 24, cent - 8) = 1;

    target(cent - 8, cent + 8) = 1;
    target(cent - 8, cent + 24) = 1;
    target(cent - 24, cent + 24) = 1;
    target(cent - 24, cent + 8) = 1;
end

function target = targetRect(sze)
    % 16x16 rectangular target pattern
    target = zeros(sze, sze);
    cent = sze/2+1; 
    rect = ones(17,17);   
    target(cent - 8 :cent + 8, cent - 8:cent + 8) = rect;
end

function [target, holf2, hol, err_arr] = binSearch(sze, hol, target)
    
    s = sqrt(sum(sum(target.^2)));
    target = target/s; % making the average "energy" 
                       % a number we can compare with
                       % the eventual hologram output

    % do the first pass through the hologram
    holf = fft2(hol);
    holf = fftshift(holf);
    holf2 = holf.*conj(holf);
    holf2 = holf2./(sqrt(sum(sum(holf2.^2))));

    err1 = sum(sum((target-holf2).^2)); % the first difference between the 
                                        % target and the hologram output.
    iter = 40000;
    err_arr = zeros(1,iter);            % an array so we can plot the error

    tol = .025; % tolerance to terminate when error is small enough
    iterCount = 0;
    for i = 1:iter
        rpix_x = round(1 + rand*(sze-1));             % pick a random x
        rpix_y = round(1 + rand*(sze-1));             % pick a random y

        hol(rpix_x, rpix_y) = -1*hol(rpix_x, rpix_y); % flip that pixel

        holf = fft2(hol);
        holf = fftshift(holf);
        holf2 = holf.*conj(holf);

        holf2 = holf2/(sqrt(sum(sum(holf2.^2))));
        err2 = sum(sum((target-holf2).^2));

        if err2 > err1 % if error increased, flip pixel back/reject change

            hol(rpix_x, rpix_y) = -1*hol(rpix_x, rpix_y);
        else
            err1 = err2; % if error hasn't increased, keep change
        end

        err_arr(i) = err1; % record the error

        if err2 < tol % if error gets low enough
            break
        end

        iterCount = iterCount + 1; % keep track of how many iterates

    end
end

function [target, holf2, hol, err_arr] = annealing(sze, hol, target)
    
    s = sqrt(sum(sum(target.^2)));
    target = target/s; % making the average "energy" 
                       % a number we can compare with
                       % the eventual hologram output

    % do the first pass through the hologram
    holf = fft2(hol);
    holf = fftshift(holf);
    holf2 = holf.*conj(holf);
    holf2 = holf2./(sqrt(sum(sum(holf2.^2))));

    err1 = sum(sum((target-holf2).^2)); % the first difference between the 
                                        % target and the hologram output.
    iter = 40000;
    err_arr = zeros(1,iter);            % an array so we can plot the error

    T = 1;
    tol = .025; % tolerance to terminate when error is small enough
    iterCount = 0;
    for i = 1:iter
        rpix_x = round(1 + rand*(sze-1)); %pick a random x
        rpix_y = round(1 + rand*(sze-1)); %pick a random y
        hol(rpix_x, rpix_y) = -1*hol(rpix_x, rpix_y); %flip that pixel
        holf = fft2(hol);
        holf = fftshift(holf);
        holf2 = holf.*conj(holf);
        holf2 = holf2/(sqrt(sum(sum(holf2.^2))));
        err2 = sum(sum((target-holf2).^2));

        if err2 > err1 % if error increased
            if exp(-1*(err2 - err1)/T) > rand % And if obj func greater than 
                                              % some random number
                hol(rpix_x, rpix_y) = -1*hol(rpix_x, rpix_y); % condition for 
                                                              % flip met
            end
        else
            err1 = err2; % keep flip if error decreased

        end
        err_arr(i) = err1; %record the error

        if err2 < tol
            break
        end

        iterCount = iterCount + 1; % keep track of how many iterates
        T = .9999*T; % lower the temperature every iter

    end
end
