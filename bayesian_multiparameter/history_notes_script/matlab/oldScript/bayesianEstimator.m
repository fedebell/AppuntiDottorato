a = 0.9;
N = 200;
resample_threshold = 0.5;
numGuesses = 1;

debug = 1;

%Reproduction of the old results
%{
linNum = 100;
s = linspace(1, linNum, linNum);
%}

s = [1 2 11 51];

phases = [0 pi/2];
%phases = [0 0];

k = 3;

numParticles = max(k*max(s)*N, 5000);

sizeSpace = 1;

weights = ones(numParticles, 1)/numParticles;

Q = eye(sizeSpace, sizeSpace);

%adaptiveBayesian(numParticles, sizeSpace, N, a, resample_threshold, Q, numGuesses, 1, s, phases, debug);

Num1 = N;
Num2 = 20;

results = zeros(Num1, Num2, 2);

step = 5;

for N=step:step:Num1
    j=N/step;
    for counter = 1:1:Num2
        [results(j, counter, 1), results(j, counter, 2)] = adaptiveBayesian(numParticles, sizeSpace, N, a, resample_threshold, Q, numGuesses, 1, s, phases, 0);
    end
end

save('results', 'results', '-v7.3');

theta = 1.8;

[est, res] = adaptiveBayesian(numParticles, sizeSpace, N, a, resample_threshold, Q, numGuesses, theta, s, phases, 1);

abs(est-theta)*res/pi

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

%Resampling algorithm with concentration towards the mean
function [newWeights, newParticleLocations] = resampling(weights, particleLocations, a, numParticles, sizeSpace)

newParticleLocations = zeros(numParticles, sizeSpace);
newWeights = zeros(numParticles, 1);

mu = estMeanCircular(weights, particleLocations, numParticles, sizeSpace);

h = sqrt(1-a^2);

newSigma = h^2*estCovarianceCircular(weights, particleLocations, numParticles, sizeSpace);

pd = makedist('Multinomial','probabilities', weights);
rng('default');

vectRandom = random(pd, numParticles, 1);

%TODO: Si può ancora ottimizzare
for i=1:1:numParticles
    newParticleLocations(i, :) = mod(mvnrnd(a*particleLocations(vectRandom(i), :)+(1-a)*mu, newSigma), 2*pi);
    newWeights(i, 1) = 1/numParticles;
end

end

%Circular resampling algorithm with concentration towards the mean
function [newWeights, newParticleLocations] = resamplingCircular(weights, particleLocations, a, numParticles, numParticles2, sizeSpace, N, s)

newParticleLocations = zeros(numParticles2, sizeSpace);

mu = estMeanCircular(weights, particleLocations, numParticles, sizeSpace);

h = sqrt(1-a^2);

newSigma = h^2*estCovarianceCircular(weights, particleLocations, numParticles, sizeSpace);

rng('default');
pd = makedist('Multinomial','probabilities', weights);

vectRandom = random(pd, numParticles2, 1);

%TODO: Si può ancora ottimizzare
for i=1:1:numParticles2
    newParticleLocations(i, :) = mod(mvnrnd(mod(angle(a*exp(1i*particleLocations(vectRandom(i), :))+(1-a)*exp(1i*mu)), 2*pi), newSigma), 2*pi);
end

newWeights = 1/numParticles2*ones(numParticles2, 1);

%newParticleLocations = randomizeLocations(newParticleLocations, numParticles2, sizeSpace, N, s);

end

%Simplified version of the resamping
function [newWeights, newParticleLocations] = resamplingSimplified(weights, particleLocations, a, numParticles, numParticles2, sizeSpace, N, s)

pd = makedist('Multinomial','probabilities', weights);
rng('default');

vectRandom = random(pd, numParticles2, 1);

newParticleLocations = particleLocations(vectRandom(:), :);
newWeights = 1/numParticles2*ones(numParticles2, 1);

%newParticleLocations = randomizeLocations(newParticleLocations, numParticles, sizeSpace, N, s);

end

%Simplified version of the resamping
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

function dist = angularDist(ang1, ang2)
    dist = pi - abs(mod((ang1-ang2), 2*pi)-pi);
end

%Importance sampling
function [newWeights, newParticleLocations] = resamplingSimplifiedImportance(weights, particleLocations, a, numParticles, sizeSpace, N, s)

mu = estMeanCircular(weights, particleLocations, numParticles, sizeSpace);
Sigma = estCovarianceCircular(weights, particleLocations, numParticles, sizeSpace);

%Multiplication for the factor matrix
piWeights = weights(:, 1).*(1+(ones(numParticles, 1)*(1/Sigma)).*(angularDist(particleLocations(:, 1), mu)).^2);

piWeights(:, 1) = piWeights(:, 1)/sum(piWeights(:, 1));

pd = makedist('Multinomial','probabilities', piWeights);
rng('default');

vectRandom = random(pd, numParticles, 1);

newWeights = weights(vectRandom(:, 1), 1)./piWeights(vectRandom(:, 1), 1);

newParticleLocations = particleLocations(vectRandom(:, 1), :);

newWeights(:, 1) = newWeights(:, 1)/sum(newWeights(:, 1));

%newParticleLocations = randomizeLocations(newParticleLocations, numParticles, sizeSpace, N, s);

end

function newParticleLocations = randomizeLocations(particleLocations, numParticles, sizeSpace, N, s)
newParticleLocations = particleLocations(:, :) + normrnd(0, (1/(N*max(s))), [numParticles, sizeSpace]);
end

%TODO: Still to optimize
function utility = utilityNV(weights, particleLocations, quantumResource, controlPhase, Q, numParticles, sizeSpace)

%Possible outcomes of the measurement

numOutcomes = 2;

varArray = zeros(numOutcomes, 1);

counter = 0;

for D=-1:2:1
    
    counter = counter + 1;
    
    hypotheticalWeights = weightsUpdate(weights, particleLocations, numParticles, D, quantumResource, controlPhase);
    
    hypotheticalMu = estMean(hypotheticalWeights, particleLocations, numParticles, sizeSpace);
    
    variance = 0;
    
    %FIXME -> In questa sommatoria ci vogliono i pesi vecchi o nuovi?
    for i=1:1:numParticles
        variance = variance - hypotheticalWeights(i, 1)*(particleLocations(i, :)-hypotheticalMu)'*Q*(particleLocations(i, :)-hypotheticalMu);
    end
    
    %Compute the marginalied likelyhood
    
    margLike = 0;
    
    for i=1:1:numParticles
        margLike = margLike + weights(i, 1)*posteriorSingleMeasurement(D, particleLocations(:, i), quantumResource, controlPhase);
    end
    
    variance = variance*margLike;
    
    varArray(counter) = variance;
end

utility = sum(varArray);

end

%Circular utility function
function utility = utilityNVCircular(weights, particleLocations, quantumResource, controlPhase, Q, numParticles, sizeSpace)

%Possible outcomes of the measurement
numOutcomes = 2;

varArray = zeros(numOutcomes, 1);

counter = 0;

for D=-1:2:1
    
    counter = counter + 1;
    
    hypotheticalWeights = weightsUpdate(weights, particleLocations, numParticles, D, quantumResource, controlPhase);
    
    variance = -estCovarianceCircular(hypotheticalWeights, particleLocations, numParticles, sizeSpace);
    
    %Compute the marginalied likelyhood
    margLike = weights(:, 1)'*posteriorSingleMeasurement(D, particleLocations(:, :), quantumResource, controlPhase);
    
    varArray(counter) = variance*margLike;
end

utility = sum(varArray);

end

%TODO: STill to optimize
function entr = utilityIG(weights, particleLocations, quantumResource, controlPhase, Q, numParticles, sizeSpace)

numOutcomes = 2;

postArray = zeros(numParticles, numOutcomes);

probDArray = zeros(numOutcomes, 1);

counter = 0;

for D=-1:2:1
    
    counter = counter + 1; 
    
    for i=1:1:numParticles
       postArray(i, counter) = posteriorSingleMeasurement(D, particleLocations(i, :), quantumResource, controlPhase); 
    end
    
    for i=1:1:numParticles
        probDArray(counter) =  probDArray(counter) + weights(i)*postArray(i, counter);
    end
    
    entr = shannonEntropy(probDArray);
    
    for i=1:1:numParticles
        entr = entr - weights(i, 1)*shannonEntropy(postArray(i, :));
    end
    
end

end

function entropy = shannonEntropy(p)
    entropy = -sum(p.*log2(p));
end

%Complete adapive Bayesian experiment design algorithm, using SMC
%approximation
%Scelgo il priore uniforme
function [est, totResources] = adaptiveBayesian(numParticles, sizeSpace, N, a, resample_threshold, Q, numGuesses, trueTheta, s, phases, debug)

totResources = 0;

weights = ones(numParticles, 1)/numParticles;

particleLocations = rand(numParticles, sizeSpace)*2*pi;

if(debug == 1)
    screenShoots = zeros(numParticles, sizeSpace+1, N);
    listEsperiment = zeros(N, 2);
end

%N is the number of experiments
for iexp=1:1:N
    
    optControls = zeros(numGuesses, 3);

    %If the number of guesses are ever going to big large, then we could
    %think of optimizing this procedure
    for iguess=1:1:numGuesses
        [quantResource_init, controlPhase_init] = guessExperiment2(iguess, s, phases);
        [optControls(iguess, 1), optControls(iguess, 2), optControls(iguess, 3)] = localOptimize(quantResource_init, controlPhase_init, weights, particleLocations, numParticles, sizeSpace, Q, s, phases);
    end
    
    quantResource_opt = 0;
    controlPhase_opt = 0;
    
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
        [weights, particleLocations] = resampling(weights, particleLocations, a, numParticles, sizeSpace);
        %[weights, particleLocations] = resamplingSimplifiedImportance(weights, particleLocations, a, numParticles, sizeSpace, N, s);
    end
    
    %ACHTUNG: Il numero di particelle cambia nel corso dell'esecuzione
    %del programma, ma il numero all'inizio è sempre numParticle
    
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
    
    %Threeshold probability for weights
    %{
    if(mod(N, 10) == 0)
        
        p_th = 1/(10*N*max(s)*numParticles);
        %p_th = 10^(-11);
        %p_th = 0;
        counter = 0;
        
        for i=1:1:numParticles
            if(weights(i, 1) > p_th)
                counter = counter + 1;
            end
        end
        
        newParticleLocations = zeros(counter, sizeSpace);
        newWeights = zeros(counter, 1);
        
        counter = 0;
        
        for i=1:1:numParticles
            if(weights(i, 1) > p_th)
                counter = counter + 1;
                newParticleLocations(counter, 1) = particleLocations(i, 1);
                newWeights(counter) = weights(i);
            end
        end
        
        tot = 0;
        for i=1:1:counter
            tot = tot + newWeights(i);
        end
        
        for i=1:1:counter
            newWeights(i) = newWeights(i)/tot;
        end
        
        weights = newWeights;
        particleLocations = newParticleLocations;
        numParticles = counter;
    end   
    %}
    
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
    
end

end

function [quantResource_init, controlPhase_init] = guessExperiment1(iexp, s, phases)

%s = [1 2 11 51];
%phases = [0 pi/2];

r1 = randi([1 size(phases, 2)],1,1);
r2 = randi([1 size(s, 2)],1,1);

quantResource_init = s(r2);
controlPhase_init = phases(r1);

end

function [quantResource_init, controlPhase_init] = guessExperiment2(iexp, s, phases)

quantResource_init = 1;
controlPhase_init = 0;

end

function [quantResource_opt, controlPhase_opt, utilityValue] = localOptimize(quantResource_init, controlPhase_init, weights, particleLocations, numParticles, sizeSpace, Q, s, phases)

%s = [1 2 11 51];
%phases = [0 pi/2];

utilityValue = -9999;

quantResource_opt = 1;
controlPhase_opt = 0;

for i=1:1:size(s, 2)
    for j=1:1:size(phases, 2)
        tmpUtility = utilityNVCircular(weights, particleLocations, s(i), phases(j), Q, numParticles, sizeSpace);
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