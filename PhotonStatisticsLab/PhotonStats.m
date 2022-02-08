%% PHYS 423 - Advanced Optics Lab 1: Photon Stats
clear all; close all; clc;

rng('shuffle')

% These can be commented out or in to run individually

poissonRoutine(1, 1000, 7, 1)
poissonRoutine(4, 1000, 25, 2)
poissonRoutine(12, 1000, 60, 3)
poissonRoutine(19, 1000, 100, 4)

boseERoutine(1, 1000, 15, 1, [0, 10])
boseERoutine(4, 2000, 50, 2, [0, 25])
boseERoutine(12, 10000, 150, 3, [0, 60])
boseERoutine(19, 10000, 300, 4, [0, 150])

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

end

function boseERoutine(givenMean, trials, countNum, figNum, xLims)
    x = rand(1, trials);
   
    % array of random meanVals following negative exponential distribution
    meanVals = -givenMean*log(x);
    nPhotons = [];
    for i = 1:trials
        ranArray = rand(1, trials);
        
        % use each meanval to determine tolerance for seeing photon
        tol = meanVals(i)/trials;
        photoArray = (ranArray < tol);
        nPhotons(i) = nnz(photoArray);
    end

    nPhotons = sort(nPhotons, 'descend');
    occur = zeros(1, countNum);

    for i = 1:length(nPhotons)
    ind = nPhotons(i);
        occur(ind+1) = occur(ind+1) + 1;
    end
    occur = occur/trials;

    interval = linspace(0, length(occur), length(occur));
    n0 = [0:countNum-1];
    bose = @(K) (1/(1 + givenMean))*(givenMean/(1 + givenMean)).^K;

    fig = figure(figNum+10);
    hold on
    plot(interval, bose(n0), 'k*')
    plot(interval, occur, 'bo')
    xlabel('Count Number')
    ylabel('Probability')
    xlim(xLims)
    lgd = legend('Bose-Einstein', 'Experiment');
    s = 'Average Number of Counts = ' + string(givenMean);
    title(lgd, s)
    hold off
    filename = 'boseEAvg' + string(givenMean) + '.png';
    saveas(fig, filename)

end
