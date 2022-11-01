a = 0.9;
N = 200;
resample_threshold = 0.5;
numGuesses = 1;

debug = 1;

%Reproduction of the old results

%{
linNum = 10;
s = linspace(1, linNum, linNum);
%}

%s = [1 2 4 8 16 32 64];

s = [1 2 11 51];

phases = [0 pi/2];
%phases = [0 0];

k = 10;

%numParticles = max(floor(k*max(s)*sqrt(N)), 3000);

numParticles = max(k*max(s)*N, 3000);

sizeSpace = 1;

weights = ones(numParticles, 1)/numParticles;

Q = eye(sizeSpace, sizeSpace);

adaptiveBayesian(numParticles, sizeSpace, N, a, resample_threshold, Q, numGuesses, 2, s, phases, debug);

%{
Num1 = N;
Num2 = 30;

results = zeros(Num1, Num2, 2);

step = 10;

theta = 1;

for N=step:step:Num1
    j=N/step;
    for counter = 1:1:Num2
        [results(j, counter, 1), results(j, counter, 2)] = adaptiveBayesian(numParticles, sizeSpace, N, a, resample_threshold, Q, numGuesses, theta, s, phases, 0);
    end
end

save('results', 'results', '-v7.3');

%}

%{
%Mean of the positions of the particles
function mu = estMean(weights, particleLocations, numParticles, sizeSpace)

mu = zeros(sizeSpace, 1);

for i=1:1:numParticles
    mu(:, 1) = mu(:, 1) + weights(i, 1)*particleLocations(i, :);
end

end

%Covariance matrix
function Sigma = estCovariance(weights, particleLocations, numParticles, sizeSpace)

mu = estMean(weights, particleLocations, numParticles, sizeSpace);

Sigma = zeros(sizeSpace, sizeSpace);


for i=1:1:numParticles
    Sigma(:, :) = Sigma(:, :) + weights(i, 1)*particleLocations(i, :)*particleLocations(i, :)';
end

Sigma(:, :) = Sigma(:, :) - mu(:, 1)*mu(:, 1)';

end
%}

%Circular mean of the positions of the particles
function angleMu = estMeanCircular(weights, particleLocations, numParticles, sizeSpace)
angleMu = real(mod(angle(exp(1i*particleLocations(:, :)')*weights(:, 1)), 2*pi));
end

%Circular variance
function Sigma = estCovarianceCircular(weights, particleLocations, numParticles, sizeSpace)
Sigma = real(sqrt(-2*log(sqrt(norm(exp(1i*particleLocations(:, :)')*weights(:, 1))))));
end

%Weights update algorithm
function newWeights = weightsUpdate(weights, particleLocations, numParticles, D, quantumResource, controlPhase)
newWeights = weights(:, 1).*posteriorSingleMeasurement(D, particleLocations(:, :), quantumResource, controlPhase);
newWeights = newWeights/sum(newWeights);
end

%Posterior probability distribution
%TODO -> Va riscrittta per considerare anche la visibilità
function prob = posteriorSingleMeasurement(D, particlePosition, quantumResource, controlPhase)

if(D == -1)
    prob = (1+cos(particlePosition*quantumResource+controlPhase))/2;
elseif(D == +1)
    prob = (1-cos(particlePosition*quantumResource+controlPhase))/2;
else
    prob = -1;
    disp("An error in 'posteriorSingleMeasurement' has occourred.");
end

end

%Circular resampling algorithm with concentration towards the mean
function [newWeights, newParticleLocations] = resamplingCircular(weights, particleLocations, a, numParticles, numParticles2, sizeSpace, N, s)

newParticleLocations = zeros(numParticles2, sizeSpace);

mu = estMeanCircular(weights, particleLocations, numParticles, sizeSpace);

h = sqrt(1-a^2);

newSigma = h^2*estCovarianceCircular(weights, particleLocations, numParticles, sizeSpace);

pd = makedist('Multinomial','probabilities', weights);

vectRandom = random(pd, numParticles2, 1);

%TODO: Si può ancora ottimizzare
for i=1:1:numParticles2
    if((particleLocations(vectRandom(i), :) - mu) > pi)
       temp = particleLocations(vectRandom(i), :) - 2*pi;
       mean = mod(a*temp+(1-a)*mu, 2*pi);
    elseif((mu - particleLocations(vectRandom(i), :)) > pi)
       temp = mu - 2*pi;
       mean = mod(a*particleLocations(vectRandom(i), :)+(1-a)*temp, 2*pi);
    else
       mean = mod(a*particleLocations(vectRandom(i), :)+(1-a)*mu, 2*pi);
    end
    newParticleLocations(i, :) = mod(mvnrnd(mean, newSigma), 2*pi);
end

newWeights = 1/numParticles2*ones(numParticles2, 1);

%newParticleLocations = randomizeLocations(newParticleLocations, numParticles2, sizeSpace, N, s);

end

%Simplified version of the resamping
function [newWeights, newParticleLocations] = resamplingSimplified(weights, particleLocations, a, numParticles, numParticles2, sizeSpace, N, s)

pd = makedist('Multinomial','probabilities', weights);

vectRandom = random(pd, numParticles2, 1);

newParticleLocations = particleLocations(vectRandom(:), :);
newWeights = 1/numParticles2*ones(numParticles2, 1);

%newParticleLocations = randomizeLocations(newParticleLocations, numParticles, sizeSpace, N, s);

end

%Simplified version of the resamping

function dist = angularDist(ang1, ang2)
    dist = pi - abs(mod((ang1-ang2), 2*pi)-pi);
end

%Importance sampling
function [newWeights, newParticleLocations] = resamplingSimplifiedImportance(weights, particleLocations, a, numParticles, sizeSpace, N, s)

mu = estMeanCircular(weights, particleLocations, numParticles, sizeSpace);
Sigma = estCovarianceCircular(weights, particleLocations, numParticles, sizeSpace);

%Multiplication for the factor matrix
piWeights = weights(:, 1).*(1+(ones(numParticles, 1)*(1/(0.1*Sigma^2))).*(angularDist(particleLocations(:, 1), mu)).^2);

piWeights(:, 1) = piWeights(:, 1)/sum(piWeights(:, 1));

pd = makedist('Multinomial','probabilities', piWeights);

vectRandom = random(pd, numParticles, 1);

newWeights = weights(vectRandom(:, 1), 1)./piWeights(vectRandom(:, 1), 1);

newParticleLocations = particleLocations(vectRandom(:, 1), :);

newWeights(:, 1) = newWeights(:, 1)/sum(newWeights(:, 1));

%newParticleLocations = randomizeLocations(newParticleLocations, numParticles, sizeSpace, N, s);

end

function [newWeights, newParticleLocations] = resamplingSimplifiedImportance2(weights, particleLocations, a, numParticles, sizeSpace, N, s)

[orderedWeights, index] = sort(weights, 'descend');

somma = 0;

%This fraction is arbitrary
f = 0.9;

counter = 1;
while ((somma < f) && (counter < numParticles))
   somma = somma + orderedWeights(counter); 
   counter = counter + 1;
end

extractionWeights1 = orderedWeights(1:(counter-1), 1)/sum(orderedWeights(1:(counter-1), 1));
extractionWeights2 = orderedWeights(counter:numParticles, 1)/sum(orderedWeights(counter:numParticles, 1));

extractionLocations1 = zeros(counter-1, sizeSpace);

extractionLocations2 = zeros(numParticles - counter+1, sizeSpace);

extractionLocations1(1:(counter-1), :) = particleLocations(index(1:(counter-1)), :);

extractionLocations2(1:(numParticles-counter+1), :) = particleLocations(index(counter:numParticles), :);

%[newW1, newLoc1] = resamplingCircular(extractionWeights1, extractionLocations1, a, counter-1, numParticles/2, sizeSpace, N, s);
[newW1, newLoc1] = resamplingSimplified(extractionWeights1, extractionLocations1, a, counter-1, numParticles/2, sizeSpace, N, s);
[newW2, newLoc2] = resamplingSimplified(extractionWeights2, extractionLocations2, a, numParticles-counter+1, numParticles/2, sizeSpace, N, s);

newW1 = (somma)*newW1;
newW2 = (1-somma)*newW2;

%newLoc2 = randomizeLocations(newLoc2, numParticles/2, sizeSpace, N, s);

newWeights = [newW1(:, 1); newW2(:, 1)];
newParticleLocations = [newLoc1(:, :); newLoc2(:, :)];

end

function newParticleLocations = randomizeLocations(particleLocations, numParticles, sizeSpace, N, s)
newParticleLocations = mod(particleLocations(:, :) + normrnd(0, (1/(sqrt(N)*max(s))), [numParticles, sizeSpace]), 2*pi);
end

%Circular utility function
function utility = utilityNVCircular(weights, particleLocations, quantumResource, quantumResourceP, controlPhase, Q, numParticles, sizeSpace)

%Possible outcomes of the measurement
numOutcomes = 2;

varArray = zeros(numOutcomes, 1);

counter = 0;

for D=-1:2:1
    
    counter = counter + 1;
    
    hypotheticalWeights = weightsUpdate(weights, particleLocations, numParticles, D, quantumResource, controlPhase);
    
    %weights or hypotheticalWeights ???
    variance = -estCovarianceCircular(hypotheticalWeights, particleLocations, numParticles, sizeSpace);
    
    %Compute the marginalied likelyhood
    margLike = weights(:, 1)'*posteriorSingleMeasurement(D, particleLocations(:, :), quantumResource, controlPhase);
    
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

particleLocations = rand(numParticles, sizeSpace)*2*pi;

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
    
    weights = weightsUpdate(weights, particleLocations, numParticles, Dexp, quantResource_opt, controlPhase_opt);
    
    squareWeights = weights.^2;
    
    if(1/(sum(squareWeights)) < resample_threshold*numParticles)
        %[weights, particleLocations] = resamplingSimplified(weights, particleLocations, a, numParticles, numParticles, sizeSpace, N, s);
        [weights, particleLocations] = resamplingSimplifiedImportance(weights, particleLocations, a, numParticles, sizeSpace, N, s);
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

est = estMeanCircular(weights, particleLocations, numParticles, sizeSpace);

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
    
    ratio = (estCovarianceCircular(weights, particleLocations, numParticles, sizeSpace))*tot/pi;
    
    disp("Ratio error/piHS: " + ratio);
    
    ratio = (estCovarianceCircular(weights, particleLocations, numParticles, sizeSpace))*sqrt(tot);
    
    disp("Ratio error/SQL: " + ratio);
    
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
    
    disp("---------------------------")
    disp(listEsperiment(:, 1)');
    
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
        tmpUtility = utilityNVCircular(weights, particleLocations, s(i), quantumResourceP, phases(j), Q, numParticles, sizeSpace);
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

p = (1+cos(trueTheta*quantumResource+controlPhase))/2;

if(rand(1) < p)
    Dexp = -1;
else
    Dexp = +1;
end

end