function [BestSol, BestCost, Iteration, funccount] = pso(CostFunction, nVar, VarMin, VarMax, MaxIt, nPop, setMutation)

    w = 1;
    wdamp = 0.99;
    c1 = 1.5;
    c2 = 2.0;

    VarSize = [1 nVar];
    VelMax = 0.1*(VarMax-VarMin);
    VelMin = -VelMax;

    empty_particle.Position = [];
    empty_particle.Cost = [];
    empty_particle.Velocity = [];
    empty_particle.Best.Position = [];
    empty_particle.Best.Cost = [];

    particle = repmat(empty_particle, nPop, 1);

    GlobalBest.Cost = inf;

    % Initialize function count
    funccount = 0;

    %%%
    %%%aynı noktalar üretilsin rassal olmasın- unifrnd ile
    %rng(0, 'twister');

    %%%
    for i = 1:nPop
        particle(i).Position = unifrnd(VarMin, VarMax, VarSize);
        particle(i).Velocity = zeros(VarSize);
        particle(i).Cost = CostFunction(particle(i).Position);
        funccount = funccount + 1;  % Increment function count
        particle(i).Best.Position = particle(i).Position;
        particle(i).Best.Cost = particle(i).Cost;

        if particle(i).Best.Cost < GlobalBest.Cost
            GlobalBest = particle(i).Best;
        end
    end

    BestCost = zeros(MaxIt, 1);
    zeroCostCount = 0;  % Initialize zero cost count

    for it = 1:MaxIt
        for i = 1:nPop
            particle(i).Velocity = w*particle(i).Velocity ...
                + c1*rand(VarSize).*(particle(i).Best.Position-particle(i).Position) ...
                + c2*rand(VarSize).*(GlobalBest.Position-particle(i).Position);

            particle(i).Velocity = max(particle(i).Velocity, VelMin);
            particle(i).Velocity = min(particle(i).Velocity, VelMax);

            particle(i).Position = particle(i).Position + particle(i).Velocity;

            IsOutside = (particle(i).Position<VarMin | particle(i).Position>VarMax);
            particle(i).Velocity(IsOutside) = -particle(i).Velocity(IsOutside);

            particle(i).Position = max(particle(i).Position, VarMin);
            particle(i).Position = min(particle(i).Position, VarMax);

            particle(i).Cost = CostFunction(particle(i).Position);
            funccount = funccount + 1;  % Increment function count

            % Mutation - Regional Domination Policy
            if rand<0.2
                if (setMutation) && (numel(particle) > 3)
                    % Implement your mutation function here if needed
                end
            end

            if particle(i).Cost < particle(i).Best.Cost
                particle(i).Best.Position = particle(i).Position;
                particle(i).Best.Cost = particle(i).Cost;

                if particle(i).Best.Cost < GlobalBest.Cost
                    GlobalBest = particle(i).Best;
                end
            end
        end

        BestCost(it) = GlobalBest.Cost;
         disp(['Iteration ' num2str(it) ': Best Cost = ' num2str(BestCost(it))]);
        w = w*wdamp;

        % Count zero costs
        if BestCost(it) == 0
            zeroCostCount = zeroCostCount + 1;
        end

        % Break if BestCost is 0
        if BestCost(it) == 0
            break;
        end
    end

    BestSol = GlobalBest;
    Iteration = it;
    output.funccount = funccount;  % Store function count in output struct
    output.BestCost = BestCost;
    output.Iteration = Iteration;
    
    % Store zero cost count in output struct
    output.zeroCostCount = zeroCostCount;
end
