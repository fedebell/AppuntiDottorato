numParticles = 1000;
sizeSpace = 1;

particleLocations = rand(sizeSpace, numParticles)*2*pi;
weights = ones(numParticles, 1)/numParticles;

mu = estMean(weights, particleLocations, numParticles, sizeSpace);
Sigma = estCovariance(weights, particleLocations, numParticles, sizeSpace);
f = @funct;
square = estFunction(weights, particleLocations, numParticles, sizeSpace, f);

a = 0.95;
Q = eye(sizeSpace, sizeSpace);
approx_ratio = 0.95;
N = 10000;
resample_threshold = 0.5;
numGuesses = 1;

debug = 1;

est = adaptiveBayesian(numParticles, sizeSpace, N, a, resample_threshold, approx_ratio, Q, numGuesses, 1, debug);


%Mean of the positions of the particles
function mu = estMean(weights, particleLocations, numParticles, sizeSpace)

mu = zeros(sizeSpace, 1);

for i=1:1:numParticles
    mu(:, 1) = mu(:, 1) + weights(i, 1)*particleLocations(:, i);
end

end

%Circular mean of the positions of the particles
function angleMu = estMeanCircular(weights, particleLocations, numParticles, sizeSpace)

mu = zeros(sizeSpace, 1);

for i=1:1:numParticles
    mu(:, 1) = mu(:, 1) + weights(i, 1)*exp(1i*particleLocations(:, i));
end

angleMu = mod(angle(mu), 2*pi);

end

%CCircular variance
function Sigma = estCovarianceCircular(weights, particleLocations, numParticles, sizeSpace)

mu = zeros(sizeSpace, 1);

for i=1:1:numParticles
    mu(:, 1) = mu(:, 1) + weights(i, 1)*exp(1i*particleLocations(:, i));
end

R = sqrt(norm(mu));

Sigma = sqrt(-2*loog(R));

end

%Espectation value of a generic function
function est = estFunction(weights, particleLocations, numParticles, sizeSpace, funct)

est = zeros(sizeSpace, 1);

for i=1:1:numParticles
    est(:, 1) = est(:, 1) + weights(i, 1)*exp(funct(particleLocations(:, i));
end

end

%Test function
function square = funct(particlePosition)

square = particlePosition.^2;

end

%Weight update algorithm
function newWeights = weightsUpdate(weights, particleLocations, numParticles, D, quantumResource, controlPhase)

newWeights = zeros(numParticles, 1);

for i=1:1:numParticles
    newWeights(i, 1) = weights(i, 1)*posteriorSingleMeasurement(D, particleLocations(:, i), quantumResource, controlPhase);
end

normConstant = sum(newWeights);

for i=1:1:numParticles
    newWeights(i, 1) = newWeights(i, 1)/normConstant;
end

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

%Resampling algorithm
function [newWeights, newParticleLocations] = resampling(weights, particleLocations, a, numParticles, sizeSpace)

newParticleLocations = zeros(sizeSpace, numParticles);
newWeights = zeros(numParticles, 1);

mu = estMean(weights, particleLocations, numParticles, sizeSpace);

h = sqrt(1-a^2);

newSigma = h^2*estCovariance(weights, particleLocations, numParticles, sizeSpace);

pd = makedist('Multinomial','probabilities', weights);
rng('default');

vectRandom = random(pd, numParticles, 1);

for i=1:1:numParticles
    j = vectRandom(i);
    modMu = a*particleLocations(:, j)+(1-a)*mu;
    newParticleLocations(:, i) = mvnrnd(modMu, newSigma);
    newWeights(i, 1) = 1/numParticles;
end

end

function [newWeights, newParticleLocations] = resamplingSimplified(weights, particleLocations, a, numParticles, sizeSpace)

newParticleLocations = zeros(sizeSpace, numParticles);
newWeights = zeros(numParticles, 1);

pd = makedist('Multinomial','probabilities', weights);
rng('default');

vectRandom = random(pd, numParticles, 1);

for i=1:1:numParticles
    newWeights(i, 1) = 1/numParticles;
    newParticleLocations(:, i) = particleLocations(:, vectRandom(i));
end

end

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
        variance = variance - hypotheticalWeights(i, 1)*(particleLocations(:, i)-hypotheticalMu)'*Q*(particleLocations(:, i)-hypotheticalMu);
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

function entr = utilityIG(weights, particleLocations, quantumResource, controlPhase, Q, numParticles, sizeSpace)

numOutcomes = 2;

postArray = zeros(numParticles, numOutcomes);

probDArray = zeros(numOutcomes, 1);

counter = 0;

for D=-1:2:1
    
    counter = counter + 1; 
    
    for i=1:1:numParticles
       postArray(i, counter) = posteriorSingleMeasurement(D, particleLocations(:, i), quantumResource, controlPhase); 
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


%Reduces the set of particles
function [newWeights, newParticleLocations] = reducedApproximation(weights, particleLocations, approx_ratio, numParticles, sizeSpace)

newWeightsparz = zeros(numParticles, 1);

newParticleLocationsparz = zeros(sizeSpace, numParticles);

numTilde = ceil(numParticles*approx_ratio);

newWeights = zeros(numTilde, 1);

newParticleLocations = zeros(sizeSpace, numTilde);

randPermutation = randperm(numParticles);

for i=1:1:numParticles
    
    newWeightsparz(i, 1) = weights(randPermutation(i));
    newParticleLocationsparz(:, i) = particleLocations(:, randPermutation(i));
    
end

[~, index] = sort(newWeightsparz, 'descend');

for i=1:1:numTilde
    newWeights(i, 1) = newWeightsparz(index(i));
    newParticleLocations(:, i) = newParticleLocationsparz(:, index(i));
end

%FIXME -> Aggiungo qui una normalizzazione dei pesi non segnata
%nell'articolo

norm = sum(newWeights);

for i=1:1:numTilde
    newWeights(i, 1) = newWeights(i, 1)/norm;
end

end

%Complete adapive Bayesian experiment design algorithm, using SMC
%approximation
%Scelgo il priore uniforme
function est = adaptiveBayesian(numParticles, sizeSpace, N, a, resample_threshold, approx_ratio, Q, numGuesses, trueTheta, debug)

weights = ones(numParticles, 1)/numParticles;

particleLocations = zeros(sizeSpace, numParticles);

if(debug == 1)
    screenShoots = zeros(numParticles, sizeSpace+1, N);
    listEsperiment = zeros(N, 2);
end

for i=1:1:numParticles
    %TODO: Da modificare
    particleLocations(:, i) = 2*pi*rand(sizeSpace);
end

%N is the number of experiments
for iexp=1:1:N
    
    if(approx_ratio ~= 1)
        
        [weights, particleLocations] = reducedApproximation(weights, particleLocations, approx_ratio, numParticles, sizeSpace);
        numParticles = size(particleLocations, 2);
        
    end
    
    optControls = zeros(numGuesses, 3);
    
    for iguess=1:1:numGuesses
        [quantResource_init, controlPhase_init] = guessExperiment2(iguess);
        [optControls(iguess, 1), optControls(iguess, 2), optControls(iguess, 3)] = localOptimize(quantResource_init, controlPhase_init, weights, particleLocations, numParticles, sizeSpace, Q);
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
    
    weights = weightsUpdate(weights, particleLocations, numParticles, Dexp, quantResource_opt, controlPhase_opt);
    
    squareWeights = zeros(numParticles, 1);
    
    for i=1:1:numParticles
        squareWeights(i, 1) = weights(i, 1)^2;
    end
    
    if(1/(sum(squareWeights)) < resample_threshold)
        [weights, particleLocations] = resamplingSimplified(weights, particleLocations, a, numParticles, sizeSpace);
    end
    
    %ACHTUNG: Il numero di particelle cambia nel corso dell'esecuzione
    %del programma, ma il numero all'inizio è sempre numParticle
    
    if(debug == 1)
        for c=1:1:numParticles
            screenShoots(c, 1, iexp) = weights(c, 1);
            for p=1:1:sizeSpace
                screenShoots(c, p+1, iexp) = particleLocations(p, c);
            end
        end
        
        listEsperiment(iexp, 1) = quantResource_opt;
        listEsperiment(iexp, 2) = controlPhase_opt;
        
    end
    
end

est = estMean(weights, particleLocations, numParticles, sizeSpace);



if(debug == 1)
    
    numPlot = 1;
    div = 1000;
    
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
            sumWeights(pos, j) = sumWeights(pos) + screenShoots(c, 1, listPositions(j));
        end
        
        plot(sumWeights(:, j));
    end
    
    s = [1 2 4 8 16 32 64 128 256 512];
    
    counterQuantum = zeros(max(s), 1);
    
    for c=1:1:N
        counterQuantum(listEsperiment(c, 1)) = counterQuantum(listEsperiment(c, 1)) + 1;
    end
    
    for c=1:1:size(s, 2)
        disp(s(c) + ": " + counterQuantum(s(c)));
    end
    
    counterPhase = zeros(2, 1);
    
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

function [quantResource_init, controlPhase_init] = guessExperiment1(iexp)

%s = [1 2 11 51];
s = [1 2 4 8 16 32 64 128 256 512];
phases = [0 pi/2];

r1 = randi([1 2],1,1);
r2 = randi([1 size(s, 2)],1,1);

quantResource_init = s(r2);
controlPhase_init = phases(r1);

end

function [quantResource_init, controlPhase_init] = guessExperiment2(iexp)

quantResource_init = 1;
controlPhase_init = 0;

end

function [quantResource_opt, controlPhase_opt, utilityValue] = localOptimize(quantResource_init, controlPhase_init, weights, particleLocations, numParticles, sizeSpace, Q)

s = [1 2 4 8 16 32 64 128 256 512];
%s = [1 2 11 51];
phases = [0 pi/2];

utilityValue = -9999;

quantResource_opt = 1;
controlPhase_opt = 0;

for i=1:1:size(s, 2)
    for j=1:1:2
        tmpUtility = utilityNV(weights, particleLocations, s(i), phases(j), Q, numParticles, sizeSpace);
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
