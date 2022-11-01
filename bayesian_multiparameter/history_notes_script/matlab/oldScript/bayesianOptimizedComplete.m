a = 0.9;
N = 500;
resample_threshold = 0.5;
numGuesses = 1;

debug = 1;

%Reproduction of the old results

%{
linNum = 10;
s = linspace(1, linNum, linNum);
%}

%s = [1 2 4 8 16 32 64];

%s = [1 2 11 51];

L = 32;

phases = (2*pi)/(L)*linspace(0, L, L);
%phases = [0 pi/2];

k = 10;

%numParticles = max(floor(k*max(s)*sqrt(N)), 3000);

numParticles = 1000;

sizeSpace = 5;

Q = diag(zeros(sizeSpace, 1));
Q(1, 1) = 1;

%adaptiveBayesian(numParticles, sizeSpace, N, a, resample_threshold, Q, numGuesses, 2, s, phases, debug);

Num1 = N;
Num2 = 50;

results = zeros(Num1, Num2, sizeSpace+1);

step = 10;

theta = 2;

for N=step:step:Num1
    j=N/step;
    for counter = 1:1:Num2
        numParticles = 1000;
        [results(j, counter, 1:sizeSpace), results(j, counter, sizeSpace+1)] = adaptiveBayesian(numParticles, sizeSpace, N, a, resample_threshold, Q, numGuesses, theta, s, phases, 0);
        save('results', 'results', '-v7.3');
    end
end

%}

%{
D = 1;

quantumResource = 2;

controlPhase = 0;

posteriorSingleMeasurement(D, particleLocations(:, 1), particleLocations(:, 2:sizeSpace), quantumResource, controlPhase, s)
%}

%Mean of the positions of the particles
function muFin = estMean(weights, particleLocations, sizeSpace)

muFin = zeros(sizeSpace, 1);

angleMu = real(mod(angle(exp(1i*particleLocations(:, 1)')*weights(:, 1)), 2*pi));
muVis = particleLocations(:, 2:sizeSpace)'*weights(:, 1);

muFin(1, 1) = angleMu;
muFin(2:sizeSpace, 1) = muVis;

end

%Covariance matrix
function Sigma = estCovariance(weights, particleLocations, numParticles, sizeSpace)

muMatrix = repmat(estMean(weights, particleLocations, sizeSpace), 1, numParticles);

locMinusMean = zeros(sizeSpace, numParticles);

locMinusMean(1, :) = angularDist(particleLocations(:, 1)', muMatrix(1, :));

locMinusMean(2:sizeSpace, :) = particleLocations(:, 2:sizeSpace)'-muMatrix(2:sizeSpace, :);

repWeights = repmat(weights', sizeSpace, 1);

Sigma = (locMinusMean.*repWeights)*locMinusMean';

end

%Weights update algorithm
function newWeights = weightsUpdate(weights, particleLocations, numParticles, D, quantumResource, controlPhase, sizeSpace, s)
newWeights = weights(:, 1).*posteriorSingleMeasurement(D, particleLocations(:, 1), particleLocations(:, 2:sizeSpace), quantumResource, controlPhase, s);
newWeights = newWeights/sum(newWeights);

end

%Posterior probability distribution
%TODO -> Va riscrittta per considerare anche la visibilità
function prob = posteriorSingleMeasurement(D, particlePosition, visibilities, quantumResource, controlPhase, s)

visArray = visibilities(:, find(s == quantumResource));

if(D == -1)
    prob = (1+visArray.*cos(particlePosition*quantumResource+controlPhase))/2;
    %prob = (1+cos(particlePosition*quantumResource+controlPhase))/2;
elseif(D == +1)
    prob = (1-visArray.*cos(particlePosition*quantumResource+controlPhase))/2;
    %prob = (1-cos(particlePosition*quantumResource+controlPhase))/2;
else
    prob = -1;
    disp("An error in 'posteriorSingleMeasurement' has occourred.");
end
    
end

%Circular resampling algorithm with concentration towards the mean
function [newWeights, newParticleLocations] = resamplingCircular(weights, particleLocations, a, numParticles, numParticles2, sizeSpace, N, s)

newParticleLocations = zeros(numParticles2, sizeSpace);

mu = estMean(weights, particleLocations, sizeSpace);

h = sqrt(1-a^2);

newSigma = h^2*estCovariance(weights, particleLocations, numParticles, sizeSpace);

pd = makedist('Multinomial','probabilities', weights);

vectRandom = random(pd, numParticles2, 1);

meanTotal = zeros(sizeSpace, 1);

%TODO: Si può ancora ottimizzare
for i=1:1:numParticles2
    
    if((particleLocations(vectRandom(i), 1) - mu(1, 1)) > pi)
       temp = particleLocations(vectRandom(i), 1) - 2*pi;
       mean = mod(a*temp+(1-a)*mu(1, 1), 2*pi);
    elseif((mu(1, 1) - particleLocations(vectRandom(i), 1)) > pi)
       temp = mu(1, 1) - 2*pi;
       mean = mod(a*particleLocations(vectRandom(i), 1)+(1-a)*temp, 2*pi);
    else
       mean = mod(a*particleLocations(vectRandom(i), 1)+(1-a)*mu(1, 1), 2*pi);
    end
    
    meanVis = a*particleLocations(vectRandom(i), 2:sizeSpace)'+(1-a)*mu(2:sizeSpace, 1);
    
    meanTotal(1, 1) = mean;
    meanTotal(2:sizeSpace, 1) = meanVis;
    
    newParticleLocations(i, :) = mvnrnd(meanTotal, newSigma);
    newParticleLocations(i, 1) = mod(newParticleLocations(i, 1), 2*pi);
    newParticleLocations(i, 2:sizeSpace) = min(1, newParticleLocations(i, 2:sizeSpace));
    newParticleLocations(i, 2:sizeSpace) = max(0, newParticleLocations(i, 2:sizeSpace));
end

newWeights = 1/numParticles2*ones(numParticles2, 1);

%newParticleLocations = randomizeLocations(newParticleLocations, numParticles2, sizeSpace, N, s);

end

%TODO: adding visibilities
%Simplified version of the resamping
function [newWeights, newParticleLocations] = resamplingSimplified(weights, particleLocations, a, numParticles, numParticles2, sizeSpace, N, s)

pd = makedist('Multinomial','probabilities', weights);

vectRandom = random(pd, numParticles2, 1);

newParticleLocations = particleLocations(vectRandom(:), :);
newWeights = 1/numParticles2*ones(numParticles2, 1);

%newParticleLocations = randomizeLocations(newParticleLocations, numParticles, sizeSpace, N, s);

end

function dist = angularDist(ang1, ang2)
    dist = pi - abs(mod((ang1-ang2), 2*pi)-pi);
end

%TODO: adding visibilities
%Importance sampling
function [newWeights, newParticleLocations] = resamplingSimplifiedImportance(weights, particleLocations, a, numParticles, sizeSpace, N, s)

mu = estMean(weights, particleLocations, sizeSpace);
Sigma = estCovariance(weights, particleLocations, numParticles, sizeSpace);

%Multiplication for the factor matrix
piWeights = weights(:, 1).*(1+(ones(numParticles, 1)*(1/(0.1*Sigma(1, 1)^2))).*(angularDist(particleLocations(:, 1), mu(1, 1))).^2);

piWeights(:, 1) = piWeights(:, 1)/sum(piWeights(:, 1));

pd = makedist('Multinomial','probabilities', piWeights);

vectRandom = random(pd, numParticles, 1);

newWeights = weights(vectRandom(:, 1), 1)./piWeights(vectRandom(:, 1), 1);

newParticleLocations = particleLocations(vectRandom(:, 1), :);

newWeights(:, 1) = newWeights(:, 1)/sum(newWeights(:, 1));

newParticleLocations = randomizeLocations(newParticleLocations, numParticles, sizeSpace, N, s);

end

function newParticleLocations = randomizeLocations(particleLocations, numParticles, sizeSpace, N, s)
newParticleLocations = mod(particleLocations(:, :) + normrnd(0, (1/(sqrt(N)*max(s))), [numParticles, sizeSpace]), 2*pi);
end

%Circular utility function
function utility = utilityNVCircular(weights, particleLocations, quantumResource, quantumResourceP, controlPhase, Q, numParticles, sizeSpace, s)

%Possible outcomes of the measurement
numOutcomes = 2;

varArray = zeros(numOutcomes, 1);

counter = 0;

for D=-1:2:1
    
    counter = counter + 1;
    
    hypotheticalWeights = weightsUpdate(weights, particleLocations, numParticles, D, quantumResource, controlPhase, sizeSpace, s);
    
    %weights or hypotheticalWeights ???
    variance = -trace(estCovariance(hypotheticalWeights, particleLocations, numParticles, sizeSpace)*Q);
    
    %Compute the marginalied likelyhood
    margLike = weights(:, 1)'*posteriorSingleMeasurement(D, particleLocations(:, 1), particleLocations(:, 2:sizeSpace), quantumResource, controlPhase, s);
    
    varArray(counter) = variance*margLike;
end

mu = 0;

utility = sum(varArray)-mu*(1-deltaFunc(quantumResource, quantumResourceP));

end

function ris = deltaFunc(quantumResource, quantumResourceP)
    if(quantumResource == quantumResourceP)
        ris = 1;
    else
        ris = 0;
    end
end

%Complete adapive Bayesian experiment design algorithm, using SMC
%approximation. The prior is uniform
function [est, totResources] = adaptiveBayesian(numParticles, sizeSpace, N, a, resample_threshold, Q, numGuesses, trueTheta, s, phases, debug)

totResources = 0;

weights = ones(numParticles, 1)/numParticles;

particleLocations(:, 1) = 2*pi*rand(numParticles, 1);
particleLocations(:, 2:sizeSpace) = rand(numParticles, sizeSpace-1);

if(debug == 1)
    screenShoots = zeros(numParticles, sizeSpace+1, N);
    listEsperiment = zeros(N, 2);
end

quantumResourceP = 1;

%N is the number of experiments

quantResource_opt = 0;
controlPhase_opt = 0;

for iexp=1:1:N
    
    optControls = zeros(numGuesses, 3);

    %If the number of guesses are ever going to big large, then we could
    %think of optimizing this procedure
    for iguess=1:1:numGuesses
        [quantResource_init, controlPhase_init] = guessExperiment2(iguess, s, phases);
        [optControls(iguess, 1), optControls(iguess, 2), optControls(iguess, 3)] = localOptimize(quantResource_init, quantumResourceP, controlPhase_init, weights, particleLocations, numParticles, sizeSpace, Q, s, phases);
    end
    
    quantumResourceP = quantResource_opt;
    
    utilityMax = -9999;
    
    for iguess=1:1:numGuesses
        
        if(optControls(iguess, 3) > utilityMax)
            
            quantResource_opt = optControls(iguess, 1);
            controlPhase_opt = optControls(iguess, 2);
            utilityMax = optControls(iguess, 3);
            
        end
    end
    
    %ExperimentSimulation
    Dexp = experimentSimulation(quantResource_opt, controlPhase_opt, trueTheta);
    
    totResources = totResources + quantResource_opt;
    
    weights = weightsUpdate(weights, particleLocations, numParticles, Dexp, quantResource_opt, controlPhase_opt, sizeSpace, s);
    
    squareWeights = weights.^2;
    
    if(1/(sum(squareWeights)) < resample_threshold*numParticles)
        [weights, particleLocations] = resamplingCircular(weights, particleLocations, a, numParticles, numParticles, sizeSpace, N, s);
    end
    
    if(debug == 1)
        for c=1:1:numParticles
            screenShoots(c, 1, iexp) = weights(c, 1);
            for p=1:1:sizeSpace
                screenShoots(c, p+1, iexp) = particleLocations(c, p);
            end
        end
        
        listEsperiment(iexp, 1) = quantResource_opt;
        listEsperiment(iexp, 2) = controlPhase_opt;
        
    end

end

est = estMean(weights, particleLocations, sizeSpace);

if(debug == 1)
    
    numPlot = 1;
    div = 10000;
    
    listPositions = zeros(numPlot, 1);
    
    hFig = zeros(numPlot, 1);
    
    for j=1:1:numPlot
        listPositions(j) = round(N/numPlot*j);
        %listPositions(j) = 2*j;
    end
    
    sumWeights = zeros(div, numPlot);
    
    for j=1:1:numPlot
        
        hFig(j) = figure('Visible', 'on');
        
        for c=1:1:numParticles
            %Costruisco l'istogramma
            pos = (floor((screenShoots(c, 2, listPositions(j)))/(2*pi/div)))+1;
            sumWeights(pos, j) = sumWeights(pos, j) + screenShoots(c, 1, listPositions(j));
        end
        
        coeff = 1;
        
        for c=1:1:div
            sumWeights(c, j) = (sumWeights(c, j))^(1/coeff);
        end
        
        plot(sumWeights(:, j));

    end
    
    counterQuantum = zeros(max(s), 1);
    
    for c=1:1:N
        counterQuantum(listEsperiment(c, 1)) = counterQuantum(listEsperiment(c, 1)) + 1;
    end
    
    tot = 0;
    
    for c=1:1:size(s, 2)
        disp(s(c) + ": " + counterQuantum(s(c)));
        tot = tot + s(c)*counterQuantum(s(c));
    end
    
    disp("Total number of used resources: " + tot);
    
    ratio = (trace(Q*estCovariance(weights, particleLocations, numParticles, sizeSpace)))*tot/pi;
    
    %disp("Ratio error/piHS: " + ratio);
    
    ratio = (trace(Q*estCovariance(weights, particleLocations, numParticles, sizeSpace)))*sqrt(tot);
    
    %disp("Ratio error/SQL: " + ratio);
    
    counterPhase = zeros(size(phases, 2), 1);
    
    for c=1:1:N
        counterPhase(floor(listEsperiment(c, 2))+1) = counterPhase(floor(listEsperiment(c, 2))+1) + 1;
    end
    
    disp("---------------------------");
    
    disp("cos: " + counterPhase(1));
    disp("sin: " + counterPhase(2));
    
    disp("---------------------------");
    
    disp("Estimator: " + est);
    disp("True Value: " + trueTheta);
    
    %disp("---------------------------")
    %disp(listEsperiment(:, 1)');
    
end

end

function [quantResource_init, controlPhase_init] = guessExperiment1(iexp, s, phases)

r1 = randi([1 size(phases, 2)],1,1);
r2 = randi([1 size(s, 2)],1,1);

quantResource_init = s(r2);
controlPhase_init = phases(r1);

end

function [quantResource_init, controlPhase_init] = guessExperiment2(iexp, s, phases)

quantResource_init = 1;
controlPhase_init = 0;

end

function [quantResource_opt, controlPhase_opt, utilityValue] = localOptimize(quantResource_init, quantumResourceP, controlPhase_init, weights, particleLocations, numParticles, sizeSpace, Q, s, phases)

utilityValue = -9999;

quantResource_opt = 1;
controlPhase_opt = 0;

for i=1:1:size(s, 2)
    for j=1:1:size(phases, 2)
        tmpUtility = utilityNVCircular(weights, particleLocations, s(i), quantumResourceP, phases(j), Q, numParticles, sizeSpace, s);
        if(tmpUtility > utilityValue)
            utilityValue = tmpUtility;
            quantResource_opt = s(i);
            controlPhase_opt = phases(j);
        end
    end
end
end

%Simulation of the experiment
function Dexp = experimentSimulation(quantumResource, controlPhase, trueTheta)

if(quantumResource == 1)
    p = (1+0.900*cos(trueTheta*quantumResource+controlPhase))/2;
end

if(quantumResource == 2)
    p = (1+0.875*cos(trueTheta*quantumResource+controlPhase))/2;
end

if(quantumResource == 11)
    p = (1+0.850*cos(trueTheta*quantumResource+controlPhase))/2;
end

if(quantumResource == 51)
    p = (1+0.765*cos(trueTheta*quantumResource+controlPhase))/2;
end


if(rand(1) < p)
    Dexp = -1;
else
    Dexp = +1;
end

end