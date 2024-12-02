function [BestSol, BestCost, nIter, funccount] = pso_hypersphere_dinamik_h_mutate_particles(CostFunction, nVar, VarMin, VarMax, MaxIt, nPop, setMutation)

%% PSO Parameters
VarSize = [1 nVar];    % Decision Variables Matrix Size
VelMax = 0.1*(VarMax-VarMin);
VelMin = -VelMax;

% Inertia Weight
w = 1;
wdamp = 0.99;
c1 = 1.5;
c2 = 2.0;
selected_h = pi/2;
hdamp = 0.9;  % H-damping factor (reduction in step size)


% Initialize Population
empty_particle.Position = [];
empty_particle.Cost = [];
empty_particle.Velocity = [];
empty_particle.Best.Position = [];
empty_particle.Best.Cost = [];

particle = repmat(empty_particle, nPop, 1);
alpha = zeros(1, nVar - 1);

GlobalBest.Cost = inf;
funccount = 0;

% Initialize Population using Sobol Sequence
sobolSeq = sobolset(nVar); % Sobol set for nVar dimensions
sobolPoints = net(sobolSeq, nPop); % Generate nPop points

%%%
%%%aynı noktalar üretilsin rassal olmasın - unifrnd ile
%rng(0, 'twister');

%%%
for i = 1:nPop
    % Normalize the Sobol points to [VarMin, VarMax]
    %particle(i).Position = VarMin + (VarMax - VarMin) .* sobolPoints(i, :);
    particle(i).Position = unifrnd(VarMin, VarMax, VarSize);
    %particle(i) = jitterMutation(particle(i), 1, nVar, 0.5, VarMin, VarMax, CostFunction);

    % Initialize Velocity
    particle(i).Velocity = zeros(VarSize);

    % Evaluation
    particle(i).Cost = CostFunction(particle(i).Position);
    funccount = funccount + 1;

    % Update Personal Best
    particle(i).Best.Position = particle(i).Position;
    particle(i).Best.Cost = particle(i).Cost;

    % Update Global Best
    if particle(i).Best.Cost < GlobalBest.Cost
        GlobalBest = particle(i).Best;
    end
end


BestCost = zeros(MaxIt, 1);

%% PSO Main Loop with Hypersphere Update
for it = 1:MaxIt
    % Dinamik olarak c1 ve c2 değerlerini güncelle
    c1 = 1.2 + 0.3 * (1 - it / MaxIt); % İlk başta 1.5, iterasyon arttıkça azalır
    c2 = 1.5 + 0.5 * (1 - it / MaxIt); % İlk başta 2.0, iterasyon arttıkça azalır

    for i = 1:nPop
        if norm(particle(i).Best.Position - particle(i).Position) == 0
            % Update Velocity
            particle(i).Velocity = w * particle(i).Velocity ...
                + c1 * rand(VarSize) .* (particle(i).Best.Position - particle(i).Position) ...
                + c2 * rand(VarSize) .* (GlobalBest.Position - particle(i).Position);

            % Apply Velocity Limits
            particle(i).Velocity = max(particle(i).Velocity, VelMin);
            particle(i).Velocity = min(particle(i).Velocity, VelMax);

            % Update Position
            particle(i).Position = particle(i).Position + particle(i).Velocity;
        else
            % Hypersphere 1
            center1 = (particle(i).Best.Position + particle(i).Position) / 2;
            radius1 = norm(particle(i).Position - particle(i).Best.Position) / 2;

            % Hypersphere 2
            center2 = GlobalBest.Position;
            radius2 = norm(particle(i).Position - GlobalBest.Position);

            % Angular coordinates calculation
            for ii = 1:nVar - 2
                d = sqrt(sum(particle(i).Position(ii+1:nVar).^2));
                alpha(ii) = acot(particle(i).Position(ii) / d);
            end
            alpha(nVar - 1) = 2 * acot((sqrt(sum(particle(i).Position(nVar-1:nVar).^2)) + particle(i).Position(nVar-1)) / particle(i).Position(nVar));
            rho = norm(particle(i).Position, 2);


            % Noktalar doğrusal mı?
            yanit = dogrusalmi(particle(i).Position,particle(i).Best.Position,GlobalBest.Position);
            if (yanit == 1)
                % I. Durum
                if (norm(GlobalBest.Position - particle(i).Best.Position,2)  > norm(GlobalBest.Position - particle(i).Position,2))
                    particle(i).Velocity = w * particle(i).Velocity + c2 * rand(VarSize) .* (GlobalBest.Position - particle(i).Position);
                    % Apply Velocity Limits
                    particle(i).Velocity = max(particle(i).Velocity, VelMin);
                    particle(i).Velocity = min(particle(i).Velocity, VelMax);

                    % Update Position
                    particle(i).Position = particle(i).Position + particle(i).Velocity;
                else
                    % II. Durum
                    newPos = calculateNewPosition(rho, alpha, pi*rand, nVar);
                    newPos2 = calculateNewPosition(rho, alpha, -pi*rand, nVar);
                    dE_dalpha1 = (CostFunction(newPos) - CostFunction(newPos2)) / (2 * h); % Açısal Türev kullanalım
                    yon = [-1 1];
                    newPosition = hiperkureUzerindeOtelemenew(center1, radius1, dE_dalpha1, particle(i).Position, yon(randi(2)));
                    particle(i).Velocity = particle(i).Position  - newPosition;
                    particle(i).Position = newPosition;
                end
            else
                % III. durum
                % dE_dalpha1 and dynamic step size adjustment
                initial_h = pi / 2;  % Starting step size
                step_adjustment_factor = 0.5;  % Factor for dynamic adjustment
                h = initial_h;
                flag=0; % Nokta hiper küre dışında
                while (h > 10^-36)
                    newPos = calculateNewPosition(rho, alpha, h, nVar);
                    newPos2 = calculateNewPosition(rho, alpha, -h, nVar);
                    dE_dalpha1 = (CostFunction(newPos) - CostFunction(newPos2)) / (2 * h); % Açısal Türev kullanalım
                    newPosition = hiperkureUzerindeOtelemenew(center1, radius1, dE_dalpha1, particle(i).Position, -1);

                    chck = isInHypersphere(newPosition, center2, radius2);
                    if chck
                        flag=1;
                        newCost = CostFunction(newPosition);
                        funccount = funccount + 1;
                        if particle(i).Cost > newCost
                            particle(i).Position = newPosition;
                            particle(i).Cost = newCost;
                            selected_h = h;
                        end
                        %break;  %Nokta hiper küre içinde
                    end
                    h=h*step_adjustment_factor;
                end


                if ~flag
                    particle(i).Velocity = w * particle(i).Velocity ...
                        + c1 * rand(VarSize) .* (particle(i).Best.Position - particle(i).Position) ...
                        + c2 * rand(VarSize) .* (GlobalBest.Position - particle(i).Position);

                    % Apply Velocity Limits
                    particle(i).Velocity = max(particle(i).Velocity, VelMin);
                    particle(i).Velocity = min(particle(i).Velocity, VelMax);

                    % Update Position
                    particle(i).Position = particle(i).Position + particle(i).Velocity;
                end

            end
        end

        % Apply Position Limits
        particle(i).Position = max(particle(i).Position, VarMin);
        particle(i).Position = min(particle(i).Position, VarMax);

        % Evaluation
        particle(i).Cost = CostFunction(particle(i).Position);
        funccount = funccount + 1;

        % Update Personal Best
        if particle(i).Cost < particle(i).Best.Cost
            particle(i).Best.Position = particle(i).Position;
            particle(i).Best.Cost = particle(i).Cost;

            % Update Global Best
            if particle(i).Best.Cost < GlobalBest.Cost
                GlobalBest = particle(i).Best;
            end
        end
    end

    % Apply Mutation if setMutation is true
    if (setMutation) %&& (rand<0.5)
        if rand < 0.5
            particle = gaussianMutation(particle, nPop, nVar, VarMin, VarMax, CostFunction);
        else
            %jitterFactor = 0.05;  % Jitter faktörünü ayarla (küçük bir değer)
            jitterFactor = 0.05 +0.45 * (1 - it / MaxIt);
            particle = jitterMutation(particle, nPop, nVar, jitterFactor, VarMin, VarMax, CostFunction);
        end
        for i=1:nPop
            if particle(i).Best.Cost < GlobalBest.Cost
                GlobalBest = particle(i).Best;
            end
        end
    end

    BestCost(it) = GlobalBest.Cost;

    % Display Iteration Information
    disp(['Iteration ' num2str(it) ': Best Cost = ' num2str(BestCost(it))]);

    % Damping Inertia Weight
    w = w * wdamp;

    % Break if BestCost is 0
    if BestCost(it) == 0
        break;
    end
end

BestSol = GlobalBest;
nIter = it;
end

function particle = jitterMutation(particle, nPop, nVar, jitterFactor, VarMin, VarMax, CostFunction)
for i = 1:nPop
    % 100 yeni pozisyon oluştur (her satır bir yeni pozisyon)
    noise = (rand(1000, nVar) - 0.5) * jitterFactor; % jitterFactor ile gürültü ekle
    newPositions = particle(i).Position + noise; % Eski pozisyona ekle

    % Pozisyonları sınırlar içinde tut
    newPositions = max(newPositions, VarMin);
    newPositions = min(newPositions, VarMax);

    % Yeni pozisyonları değerlendir
    costs = zeros(40, 1);
    for j = 1:40
        costs(j) = CostFunction(newPositions(j, :)); % Her yeni pozisyonun maliyetini hesapla
    end

    % En iyi pozisyonu ve maliyeti seç
    [~, bestIdx] = min(costs); % En düşük maliyete sahip olanın indeksini al
    bestCost = costs(bestIdx);
    bestPosition = newPositions(bestIdx, :);

    % % En iyi pozisyonu kaydet
    if bestCost<particle(i).Cost
        particle(i).Position = bestPosition;
        particle(i).Cost = bestCost;
        % Personal best'i de güncelle
        if bestCost < particle(i).Best.Cost
            particle(i).Best.Position = bestPosition;
            particle(i).Best.Cost = bestCost;
        end
    end
end
end



function particle = gaussianMutation(particle, nPop, nVar, VarMin, VarMax, CostFunction)
for i = 1:nPop
    % 100 yeni pozisyon oluştur (her satır bir yeni pozisyon)
    noise = randn(1000, nVar) * 0.1; % Gaussian gürültü (standart sapma 0.1)
    newPositions = particle(i).Position + noise; % Eski pozisyona ekle

    % Pozisyonları sınırlar içinde tut
    newPositions = max(newPositions, VarMin);
    newPositions = min(newPositions, VarMax);

    % Yeni pozisyonları değerlendir
    costs = zeros(40, 1);
    for j = 1:40
        costs(j) = CostFunction(newPositions(j, :)); % Her yeni pozisyonun maliyetini hesapla
    end

    % En iyi pozisyonu ve maliyeti seç
    [~, bestIdx] = min(costs); % En düşük maliyete sahip olanın indeksini al
    bestCost = costs(bestIdx);
    bestPosition = newPositions(bestIdx, :);

    % % En iyi pozisyonu kaydet
    if bestCost<particle(i).Cost
        particle(i).Position = bestPosition;
        particle(i).Cost = bestCost;
        % Personal best'i de güncelle
        if bestCost < particle(i).Best.Cost
            particle(i).Best.Position = bestPosition;
            particle(i).Best.Cost = bestCost;
        end
    end
    
    
end
end




%DENEDİM 10 tane oluştur
% function particle = jitterMutation(particle, nPop, nVar, jitterFactor, VarMin, VarMax, CostFunction)
% for i = 1:nPop
%     % Birkaç jitter noktası oluştur
%     nMutate = 10; % Kaç yeni nokta üretilecek
%     bestCost = particle(i).Cost;
%     bestPosition = particle(i).Position;
%     
%     for j = 1:nMutate
%         % Jitter noise: küçük bir rastgele değer
%         noise = (rand(1, nVar) - 0.5) * jitterFactor; % jitterFactor jitter büyüklüğü için
%         newPosition = particle(i).Position + noise;
% 
%         % Pozisyonun sınırlar içinde olduğundan emin ol
%         newPosition = max(newPosition, VarMin);
%         newPosition = min(newPosition, VarMax);
% 
%         % Yeni pozisyonu değerlendir
%         newCost = CostFunction(newPosition);
% 
%         % Yeni pozisyon daha iyi ise güncelle
%         if newCost < bestCost
%             bestCost = newCost;
%             bestPosition = newPosition;
%         end
%     end
% 
%     % En iyi pozisyonu kaydet
%     particle(i).Position = bestPosition;
%     particle(i).Cost = bestCost;
% 
%     % Personal best'i de güncelle
%     if bestCost < particle(i).Best.Cost
%         particle(i).Best.Position = bestPosition;
%         particle(i).Best.Cost = bestCost;
%     end
% end
% end
% 
% 
% function particle = gaussianMutation(particle, nPop, nVar, VarMin, VarMax, CostFunction)
% for i = 1:nPop
%     % Birkaç Gaussian noktası oluştur
%     nMutate = 10; % Kaç yeni nokta üretilecek
%     bestCost = particle(i).Cost;
%     bestPosition = particle(i).Position;
%     
%     for j = 1:nMutate
%         % Gaussian noise
%         noise = randn(1, nVar) * 0.1; % Standart sapma 0.1
%         newPosition = particle(i).Position + noise;
% 
%         % Pozisyonun sınırlar içinde olduğundan emin ol
%         newPosition = max(newPosition, VarMin);
%         newPosition = min(newPosition, VarMax);
% 
%         % Yeni pozisyonu değerlendir
%         newCost = CostFunction(newPosition);
% 
%         % Yeni pozisyon daha iyi ise güncelle
%         if newCost < bestCost
%             bestCost = newCost;
%             bestPosition = newPosition;
%         end
%     end
% 
%     % En iyi pozisyonu kaydet
%     particle(i).Position = bestPosition;
%     particle(i).Cost = bestCost;
% 
%     % Personal best'i de güncelle
%     if bestCost < particle(i).Best.Cost
%         particle(i).Best.Position = bestPosition;
%         particle(i).Best.Cost = bestCost;
%     end
% end
% end




%ORJİNAL

% function particle = jitterMutation(particle, nPop, nVar, jitterFactor, VarMin, VarMax, CostFunction)
% for i = 1:nPop
%     % Jitter noise: küçük bir rastgele değer
%     noise = (rand(1, nVar) - 0.5) * jitterFactor; % jitterFactor jitter büyüklüğü için
%     newPosition = particle(i).Position + noise;
% 
%     % Pozisyonun sınırlar içinde olduğundan emin ol
%     newPosition = max(newPosition, VarMin);
%     newPosition = min(newPosition, VarMax);
% 
%     % Yeni pozisyonu değerlendir
%     newCost = CostFunction(newPosition);
% 
%     % Eğer yeni pozisyon daha iyiyse parçacığı güncelle
%     if newCost < particle(i).Cost
%         particle(i).Position = newPosition;
%         particle(i).Cost = newCost;
%         % Personal best'i de güncelle
%         if newCost <particle(i).Best.Cost
%             particle(i).Best.Position = newPosition;
%             particle(i).Best.Cost = newCost;
%         end
%     end
% end
% end
% 
% 
% 
% function particle = gaussianMutation(particle, nPop, nVar, VarMin, VarMax, CostFunction)
% for i = 1:nPop
%     % Gaussian noise
%     noise = randn(1, nVar) * 0.1; % Standart sapma 0.1
%     %noise = randn(floor(nPop*0.4), nVar) * 0.1; % Standart sapma 0.1
%     newPosition = particle(i).Position + noise;
% 
%     % Ensure position is within bounds
%     newPosition = max(newPosition, VarMin);
%     newPosition = min(newPosition, VarMax);
% 
%     % Evaluate new position
%     newCost = CostFunction(newPosition);
% 
%     % Update particle position and cost if the new position is better
%     if newCost < particle(i).Cost
%         particle(i).Position = newPosition;
%         particle(i).Cost = newCost;
%         % Update personal best
%         if newCost <particle(i).Best.Cost
%             particle(i).Best.Position = newPosition;
%             particle(i).Best.Cost = newCost;
%         end
%     end
% end
% end
