%Copyright (c) 2024, Iclal Gor

%All rights reserved.

%Developer: Iclal Gor (Aydın Adnan Menderes University, Department of Mathematics)

%Contact Info: iclal@adu.edu.tr

clear; close; clc;

T0 = getBaseTime;


nVar = 30;           % Number of Decision Variables

VarSize = [1 nVar];   % Size of Decision Variables Matrix

VarMin = -1.28; % Lower Bound of Variables
VarMax = 1.28; % Upper Bound of Variables

%VarMin=[12 12 12 12]; %GearTrainDesign nVar = 4;
%VarMax=[60 60 60 60];

%% PSO Parameters

MaxIt = 1;          % Maximum Number of Iterations
nPop = 100;           % Population Size (Swarm Size)

% PSO Parameters
w = 1;                % Inertia Weight
wdamp = 0.99;         % Inertia Weight Damping Ratio
c1 = 1.5;             % Personal Learning Coefficient
c2 = 2.0;             % Global Learning Coefficient

% Velocity Limits
VelMax = 0.1 * (VarMax - VarMin);
VelMin = -VelMax;



%% Initialization
empty_particle.Position = [];
empty_particle.Cost = [];
empty_particle.Velocity = [];
empty_particle.Best.Position = [];
empty_particle.Best.Cost = [];

particle = repmat(empty_particle, nPop, 1);

GlobalBest.Cost = inf;

method = 'PSO';
maxTrials = 25;

setMutation = 1;

%% PSO without Hypersphere Dynamics and Mutation (PSO-original)
Sol_0 = repmat(struct('fmin', [], 'nFeval', [], 'nIter', [], 'elapsedTime', [], 'zeroCostCount', []), maxTrials, 1);
fmin = zeros(1, maxTrials);
nFeval = zeros(1, maxTrials);
nIter = zeros(1, maxTrials);
T = zeros(1, maxTrials);

fprintf('PSO Original \n');
zeroCostCount = 0;  % Initialize zero cost count

Func_name = 'F13'; % Fitness Function
[VarMin,VarMax,dim,fobj]=Get_Functions_details(Func_name);

%Func_name = 'GearTrainDesign';
%[VarMin,VarMax,dim,fobj]=Get_Engineering_Problems_details(Func_name);

for k = 1:maxTrials
    fprintf("\nTrial Id: %d\n",k);
    tic;
    [BestSol, BestCost, nIter, funccount] = pso(fobj, nVar, VarMin, VarMax, MaxIt, nPop, setMutation);
    T(k) = toc;
    Sol_0(k).fmin = BestSol.Cost;
    Sol_0(k).BestCost = BestCost;
    Sol_0(k).nFeval = funccount;
    Sol_0(k).nIter = nIter;
    Sol_0(k).elapsedTime = T(k);
end

T1 = mean(T);
T2 = (T1 - min(T)) / T0;

%% PSO with Hypersphere Dynamics and Mutation (PSO - HDM)
fprintf('PSO-HDM \n');
Sol_1 = repmat(struct('fmin', [], 'nFeval', [], 'nIter', [], 'elapsedTime', [], 'zeroCostCount', []), maxTrials, 1);
fmin_hypersphere = zeros(1, maxTrials);
nFeval_hypersphere = zeros(1, maxTrials);
nIter_hypersphere = zeros(1, maxTrials);
T_hypersphere = zeros(1, maxTrials);
zeroCostCount_hypersphere = 0;  % Initialize zero cost count

for k = 1:maxTrials
    tic;
   [BestSol_hypersphere, BestCost_hypersphere, nIter, funccount] = pso_hypersphere_dinamik_h_mutate_particles(fobj, nVar, VarMin, VarMax, MaxIt, nPop, setMutation);
     T_hypersphere(k) = toc;
    Sol_1(k).fmin_hypersphere = BestSol_hypersphere.Cost;
    Sol_1(k).BestCost = BestCost_hypersphere;
    Sol_1(k).nFeval = funccount;
    Sol_1(k).nIter = nIter;
    Sol_1(k).elapsedTime = T_hypersphere(k);
end

T1_hypersphere = mean(T_hypersphere);
T2_hypersphere = (T1_hypersphere - min(T_hypersphere)) / T0;



fprintf('\nPSO Original \n');
% Display results
fmins = [Sol_0(:).fmin];
[~, minInd_0] = min(fmins);
[~, maxInd_0] = max(fmins);
fprintf('%s after %d iterations: %d\n\n', method, Sol_0(minInd_0).nIter);

% Algorithm complexity measures
fprintf('Algorithm complexity measures:\n');
fprintf('T1: %8.6e\n', T1);
fprintf('T2 = (mean(T) - min(T)) / T0: %8.6e depending on base time T0: %8.6e seconds\n\n', T2, T0);

% Display Performance Scores
fprintf('Best of Score obtained from %s: %8.6e\n', method, min(fmins));
fprintf('Worst of Score obtained from %s: %8.6e\n', method, max(fmins));
fprintf('Mean of Scores obtained from %s: %8.6e +/- %8.6e\n', method, mean(fmins), std(fmins));


fprintf('\nPSO WITH HYPERSPHERE DYNAMICS AND MUTATION (PSO-HDM)\n');
% Display results
fmins_hypersphere = [Sol_1(:).fmin_hypersphere];
[~, minInd_1] = min(fmins_hypersphere);
[~, maxInd_1] = max(fmins_hypersphere);
fprintf('%s after %d iterations: %d\n\n', method, Sol_1(minInd_1).nIter);

% Algorithm complexity measures
fprintf('Algorithm complexity measures:\n');
fprintf('T1: %8.6e\n', T1_hypersphere);
fprintf('T2 = (mean(T) - min(T)) / T0: %8.6e depending on base time T0: %8.6e seconds\n\n', T2_hypersphere, T0);

% Display Performance Scores
fprintf('Best of Score obtained from %s: %8.6e\n', method, min(fmins_hypersphere));
fprintf('Worst of Score obtained from %s: %8.6e\n', method, max(fmins_hypersphere));
fprintf('Mean of Scores obtained from %s: %8.6e +/- %8.6e\n', method, mean(fmins_hypersphere), std(fmins_hypersphere));

figure;
plot(Sol_0(minInd_0).BestCost, 'LineWidth', 2);
hold on;
plot(Sol_1(minInd_1).BestCost, 'LineWidth', 2);
xlabel('Trial');
ylabel('Best Cost');
legend({'PSO', 'PSO-HDM'});
title('Best Cost of Algorithms',Func_name);
grid on;

figure;
plot(1:maxTrials, T, 'r-', 'LineWidth', 2);
hold on;
plot(1:maxTrials, T_hypersphere, 'b-', 'LineWidth', 2);
yline(mean(T), 'r--', 'LineWidth', 2);  % Ortalama zaman yatay çizgi
yline(mean(T_hypersphere), 'b--', 'LineWidth', 2);  % Ortalama zaman yatay çizgi
xlabel('Trial');
ylabel('Elapsed Time (seconds)');
title('Elapsed Time Comparison');
legend({'Without Hypersphere Dynamics and Mutation', 'With Hypersphere Dynamics and Mutation', 'Mean Without Hypersphere Dynamics and Mutation', 'Mean With Hypersphere Dynamics and Mutation'});
grid on;


