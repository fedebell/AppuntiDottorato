%Application of the Bayesian algorithm of Granade to the estimation of a phase with unknown visibilities

%maximum number of resources of the simulation
N = 30000;

%Parameters of the resampling algorithm
a = 0.9;

resample_threshold = 0.5;
numGuesses = 1;

%Debug mode
debug = 1;

s = [1 2 11 51];

%Possible polarization measurement
phases = [0 pi/2];

%Extendend polarization measurements
%L = 32;
%phases = (2*pi)/(L)*linspace(0, L, L);

%Number of particles of the particle filter
numParticles = 1000;

%Size of the space of parameters in which the estimation is performed
sizeSpace = 5;

%Weight of the covariance matrix appearing in the utility function
%We take into account only the variance of the phase, this means that the visibilities are nuisance parameters
Q = zeros(sizeSpace);
Q(1, 1) = 1;
Q(3, 3) = 1;

%True phase that we need to estimate
theta = 0.38000;

%Single adaptive Bayesian estimation
%array = adaptiveBayesian(numParticles, sizeSpace, N, a, resample_threshold, Q, numGuesses, theta, s, phases, debug);

%The following code is used to produce a simulation and save the results in a variable results

%Maximum number of resources used
Num1 = N;

%Number of repetition of each experiment with a fixed resource number
Num2 = 50;

NumTheta = 17;

thetaArray = [0.00235
0.06145
0.29345
0.21635
0.38000
1.13975
1.26430
0.49620
0.46280
1.44890
1.66450
1.75520
1.87500
2.12710
2.58995
2.74000
2.96000];

%The first index is the number of resource used in the estimation, the second index the repetition of the experiment, the third index conains the estimatimated phase and visibilities + the used number of total resources, which will be higher but very close to Num1.
%results(N, j, :) = (theta, V1, V2, V3, V4, totResources) refers to the jth execution of an experiment with N resources (more precisely totResources) that gave (theta, V1, V2, V3, V4) as output.
results = zeros(NumTheta, Num1, Num2, sizeSpace);


numCores = feature('numcores');
p = parpool(numCores);

parfor k=1:1:17
    for i=1:1:Num2
        array = adaptiveBayesian(numParticles, sizeSpace, N, a, resample_threshold, Q, numGuesses, thetaArray(k), s, phases, 0);
        for j=1:1:Num1
            results(k, j, i, :) = array(:, j);
        end
        disp(k + ":" + i);
    end
end

save('results', 'results', '-v7.3');

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
function prob = posteriorSingleMeasurement(D, particlePosition, visibilities, quantumResource, controlPhase, s)

visArray = visibilities(:, find(s == quantumResource));

if(D == -1)
    prob = (1+visArray.*cos(particlePosition*quantumResource+controlPhase))/2;
elseif(D == +1)
    prob = (1-visArray.*cos(particlePosition*quantumResource+controlPhase))/2;
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

%Theoretically we could chose a differet number of particles after the resamplin. In our simualtion we choose always numParticles2 = numParticles1
for i=1:1:numParticles2

    %This code compute the linear combination between two angles
    if((particleLocations(vectRandom(i), 1) - mu(1, 1)) > pi)
       temp = particleLocations(vectRandom(i), 1) - 2*pi;
       mean = mod(a*temp+(1-a)*mu(1, 1), 2*pi);
    elseif((mu(1, 1) - particleLocations(vectRandom(i), 1)) > pi)
       temp = mu(1, 1) - 2*pi;
       mean = mod(a*particleLocations(vectRandom(i), 1)+(1-a)*temp, 2*pi);
    else
       mean = mod(a*particleLocations(vectRandom(i), 1)+(1-a)*mu(1, 1), 2*pi);
    end
    
    %Linear combination of visibilities
    meanVis = a*particleLocations(vectRandom(i), 2:sizeSpace)'+(1-a)*mu(2:sizeSpace, 1);
    
    meanTotal(1, 1) = mean;
    meanTotal(2:sizeSpace, 1) = meanVis;
    
    %Extraction of the new particles location from a Gaussian
    newParticleLocations(i, :) = mvnrnd(meanTotal, newSigma);
    
    %The angle is cast in [0, 2 \pi)
    newParticleLocations(i, 1) = mod(newParticleLocations(i, 1), 2*pi);
    
    %The visibilities extracted from a Gaussian are casted in [0, 1]
    newParticleLocations(i, 2:sizeSpace) = min(1, newParticleLocations(i, 2:sizeSpace));
    newParticleLocations(i, 2:sizeSpace) = max(0, newParticleLocations(i, 2:sizeSpace));
end

%The new particles are weighted equally
newWeights = 1/numParticles2*ones(numParticles2, 1);

end

%Function that computes the angular distance 
function dist = angularDist(ang1, ang2)
    dist = pi - abs(mod((ang1-ang2), 2*pi)-pi);
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
    
    %This is the hypothetical variance
    variance = -trace(estCovariance(hypotheticalWeights, particleLocations, numParticles, sizeSpace)*Q);
    
    %Compute the marginalied likelyhood
    margLike = weights(:, 1)'*posteriorSingleMeasurement(D, particleLocations(:, 1), particleLocations(:, 2:sizeSpace), quantumResource, controlPhase, s);
    
    varArray(counter) = variance*margLike;
end

%We tried to force the measurements with the same q-plate to be performed in sequence, but we didn't obtain good results.
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

%Complete adaptive Bayesian experiment design algorithm, using sequential Montecarlo. The prior is uniform.
function estArray = adaptiveBayesian(numParticles, sizeSpace, N, a, resample_threshold, Q, numGuesses, trueTheta, s, phases, debug)

%Maximum number of experiment to be performed, that is maximum number of photons to be used.
maxExp = 10000;

%Total number of resources
totResources = 0;
estArray = zeros(sizeSpace, N);

%The weights are initially uniform.
weights = ones(numParticles, 1)/numParticles;

%The particle locations are initially chosen at random
particleLocations(:, 1) = 2*pi*rand(numParticles, 1);
particleLocations(:, 2:sizeSpace) = rand(numParticles, sizeSpace-1);

%0.900, 0.875, 0.850, 0.765,

%{
for p=1:1:numParticles
    particleLocations(:, 2) = 0.900;
    particleLocations(:, 3) = 0.875;
    particleLocations(:, 4) = 0.850;
    particleLocations(:, 5) = 0.765;
end
%}

%Quantum resource (q-plate) to be used at each experiment
quantumResourceP = 1;


quantResource_opt = 0;
controlPhase_opt = 0;

%Counter of the number of experiments (total photons used)
iexp = 0;

%While the total number of consumed resources is still less than N and the total number of experiments has not been reaches 
while((totResources < N) && (iexp < maxExp))
    
    iexp = iexp + 1;
    
    optControls = zeros(numGuesses, 3);

    %We compute here the optimal control, that is the optimal q-plate to be used in the next experiment and the polarization measurement to be performed
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
    
    %Update of the total number of used resources
    totResources = totResources + quantResource_opt;
    
    %Update of the weights.
    weights = weightsUpdate(weights, particleLocations, numParticles, Dexp, quantResource_opt, controlPhase_opt, sizeSpace, s);
    
    squareWeights = weights.^2;
    
    %If the particles are not used optimally, then we resample
    if(1/(sum(squareWeights)) < resample_threshold*numParticles)
        [weights, particleLocations] = resamplingCircular(weights, particleLocations, a, numParticles, numParticles, sizeSpace, N, s);
    end
    
    %Debug instructions
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
    
    %Mean of the particle distributions
    estArray(:, totResources) = estMean(weights, particleLocations, sizeSpace);


end

%Number of photons that were used to reach the total number of resources totResoiurces
numExp = iexp;

%Mean of the particle distributions
est = estMean(weights, particleLocations, sizeSpace);

%Debug code
if(debug == 1)
    
    numPlot = 1;
    div = 10000;
    
    listPositions = zeros(numPlot, 1);
    
    hFig = zeros(numPlot, 1);
    
    for j=1:1:numPlot
        listPositions(j) = round(numExp/numPlot*j);
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
    
    for c=1:1:numExp
        counterQuantum(listEsperiment(c, 1)) = counterQuantum(listEsperiment(c, 1)) + 1;
    end
    
    tot = 0;
    
    for c=1:1:size(s, 2)
        disp(s(c) + ": " + counterQuantum(s(c)));
        tot = tot + s(c)*counterQuantum(s(c));
    end
    
    disp("Total number of used resources: " + tot);
    
    ratio = sqrt((trace(Q*estCovariance(weights, particleLocations, numParticles, sizeSpace))))*tot/pi;
    
    disp("Ratio error/piHS: " + ratio);
    
    ratio = sqrt((trace(Q*estCovariance(weights, particleLocations, numParticles, sizeSpace))))*sqrt(tot);
    
    disp("Ratio error/SQL: " + ratio);
    
    counterPhase = zeros(size(phases, 2), 1);
    
    for c=1:1:numExp
        counterPhase(floor(listEsperiment(c, 2))+1) = counterPhase(floor(listEsperiment(c, 2))+1) + 1;
    end
    
    %This representation works only if the possible phases are just 0 and pi/2
    
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

%Random experiment choice
function [quantResource_init, controlPhase_init] = guessExperiment1(iexp, s, phases)

r1 = randi([1 size(phases, 2)],1,1);
r2 = randi([1 size(s, 2)],1,1);

quantResource_init = s(r2);
controlPhase_init = phases(r1);

end

%Trivial choice of the experiment
function [quantResource_init, controlPhase_init] = guessExperiment2(iexp, s, phases)

quantResource_init = 1;
controlPhase_init = 0;

end

%Computation of the utility value, that is the expected variance of the posterior probability distribution
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

%The visibilities are those of the experiment
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
