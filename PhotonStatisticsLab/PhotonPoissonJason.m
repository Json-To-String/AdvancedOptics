%% PHYS 423 - Advanced Optics Lab 1: Photon Stats - Poisson
clear all; close all; clc;

rng('shuffle')

poissonRoutine(1, 1000, 7, 1)
poissonRoutine(4, 1000, 25, 2)
poissonRoutine(12, 1000, 60, 3)
poissonRoutine(19, 1000, 100, 4)

function poissonRoutine(givenMean, trials, countNum, figNum)
    nPhotons = [];
    tol = givenMean/trials;

    % A given trial is defined within the for loop
    for i = 1:trials

        % an array of random numbers at a large length is made
        ranArray = rand(1, 1000);

        % a photon event is detected if the element in the array is below the 
        % tolerance

        % an array of logical values is made to reflect which indices yield a
        % detection event (True = 1, False = 0)
        photoArray = (ranArray < tol);

        % the number of nonzero elements from the array of logical values to is
        % reflect the number of photons detected from a given trial
        nPhotons(i) = nnz(photoArray);
    end

    nPhotons = sort(nPhotons, 'descend');

    % For plotting, we want to take the number of occurrences of a given count
    % number and divide it by the number of trials

    % array of zeros to reflect each given count number, with the count number
    % at each index starting from 1
    occur = zeros(1, countNum);
    for i = 1:length(nPhotons)
        % we use the value from each trial as the index to increment the count
        % number
        ind = nPhotons(i);
        occur(ind+1) = occur(ind+1) + 1;
    end
    % division happens after the loop
    occur = occur/trials;

    % This will start at zero because zero counts is a reasonable occurrence
    interval = linspace(0, length(occur), length(occur));
    n0 = [0:countNum-1];
    
    poisson = @(n) ((givenMean.^n).*exp(-givenMean))./(factorial(n));
    fig = figure(figNum);
    hold on
    plot(interval, poisson(n0), 'k*')
    plot(interval, occur, 'bo')
    xlabel('Count Number')
    ylabel('Probability')
    lgd = legend('Poisson', 'Experiment');
    s = 'Average Number of Counts = ' + string(givenMean);
    title(lgd, s)
    hold off
    filename = 'PoissonAvg' + string(givenMean) + '.png';
    saveas(fig, filename)
%     close all;

end
