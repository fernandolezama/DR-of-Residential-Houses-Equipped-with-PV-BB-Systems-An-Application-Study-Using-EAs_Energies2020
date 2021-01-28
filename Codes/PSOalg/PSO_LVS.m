function [Fit_and_p,FVr_bestmemit, fitMaxVector,Best_xOpt]=PSO_LVS(PSOparameters,Data,lowerB,upperB,id,iRuns)
%% Algorithm Information 
% Created by Ricardo Faia on 19/12/2019 %
% Original PSO provided from https://yarpiz.com/ and adpated by Ricardo Faia %
% Original VS provided from https://web.itu.edu.tr/~bdogan/VortexSearch/VS.htm %
% PSO-LVS (Particle Swarm Optimization with Local Vortex Search) %

%% Problem Definition
lowerlimit =    lowerB; %lower limit of the problem
upperlimit =    upperB; %upper limit of the problem

nVar=           length(upperB);            % Number of Decision Variables

VarMin=         lowerlimit;         % Lower Bound of Variables
VarMax=         upperlimit;         % Upper Bound of Variables

%% PSO Parameters

MaxIt=                      PSOparameters.MaxIt;       % Maximum Number of Iterations
nPop=                       PSOparameters.nPop;        % Population Size (Swarm Size)
wmin=                       PSOparameters.wmin;        % Inertia Weight
wmax=                       PSOparameters.wmax;
c1=                         PSOparameters.c1;          % Personal Learning Coefficient
c2=                         PSOparameters.c2;          % Global Learning Coefficient
BRM=                        PSOparameters.BRM;
pPSOinit=                   1; 
alpha=                      0.999;

%% Initialization
Cost=[];%vectores
LocalBest.Position=[];
LocalBest.Cost=[];
GlobalBest.Position=[];
GlobalBest.Cost=[];
if length(VarMin)==1
    minPositionsMatrix=repmat(VarMin,nPop,nVar);
    maxPositionsMatrix=repmat(VarMax,nPop,nVar);
else
    minPositionsMatrix=repmat(VarMin,nPop,1); 
    maxPositionsMatrix=repmat(VarMax,nPop,1);
end
% generate initial population

rand('state',iRuns) %Guarantee same initial population

Position=genpop(nPop,nVar,minPositionsMatrix,maxPositionsMatrix); 

Velocity=0.1*Position;  

% Evaluation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FM_pop=Position;
% load('Tomlab_sol');
% FM_pop(1,:)=Results.x;
%------Evaluate the best member after initialization----------------------
for i=1:nPop
    [xOpt(i)] = Create_Results_str(Data);
    [xOpt(i)] = decoding_sol(id,FM_pop(i,:),Data,xOpt(i));
    [S_val(i,1),xOpt(i)] = Fit_Evaluation(xOpt(i),Data,id);
%    [FM_pop(i,:)] = encoding_repaired_sol(xOpt(i),FM_pop(i,:),id,Data);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[S_bestval,I_best_index] = min(S_val); % This mean that the best individual correspond to the best worst performance
Best_xOpt=xOpt(I_best_index);
FVr_bestmemit = FM_pop(I_best_index,:); % best member of current iteration
fitMaxVector(1,1) = S_bestval;
Position=FM_pop;
Cost=S_val;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

LocalBest.Position=Position; 
LocalBest.Cost(:,1)=Cost; 
inx=Cost<LocalBest.Cost; 
LocalBest.Cost(inx)=Cost(inx); 
[~,idk]=min(LocalBest.Cost); 
GlobalBest.Position=LocalBest.Position(idk,:); 
GlobalBest.Cost=LocalBest.Cost(idk,1); 
BestSol.Result(1,1)=GlobalBest.Cost;  
BestCost(1)=GlobalBest(1).Cost; 
pPSO(1)=pPSOinit;
evals=nPop;
fprintf('Fitness value: %f\n',BestCost(1))
fprintf('Iteration: %d\n',1)
%% PSO Main Loop
for it=2:MaxIt
    %Inertia update - updata at each iteration
    w=wmax-((wmax-wmin)/MaxIt)*it;
    % Velocity Limits - updata at each iteration - optional
    f=0.5-((0.5-0.01)/MaxIt)*it;
    VelMax=repmat(f*(VarMax-VarMin),nPop,1);
    VelMin=-VelMax;
    %random criation
    r1=rand;
    r2=rand;
    % Update Velocity
    %Velocity=w.*Velocity+...
    %            c1*r1.*(LocalBest.Position-Position)+...
    %            c2*r2.*(repmat(GlobalBest(it-1).Position,[nPop,1])-Position);
    Velocity=w.*Velocity+...
                c1*r1.*(LocalBest.Position-Position)+...
                c2*r2.*(repmat(GlobalBest(it-1).Position,[nPop,1]).*(1+randn([nPop, nVar]))-Position);          
    % Apply Velocity Limits
    Velocity =max(Velocity,VelMin);
    Velocity =min(Velocity,VelMax);
        
    %Roulette Wheel Selection 
    %pPSO(it)=alpha*pPSO(it-1); %Based on exponencial decreasing 
    pPSO(it)=0.9^(it/8);
    pOther=1-pPSO(it);
    Prob=[pPSO(it),pOther];
    r=rand; % 0-1
    Prob_vec=cumsum(Prob); 
    Method=find(r<=Prob_vec,1,'first');
    %Method=1;
    if Method==2 && it>0.90*MaxIt
        Position = OtherMethodSearch(GlobalBest(it-1).Position,nPop,MaxIt,it,nVar,minPositionsMatrix,maxPositionsMatrix);        
    else
        Position = Position + Velocity; % Update Velocity with velocity PSO function
        % Velocity Mirror Effect - Optinal
        IsOutside=(Position<minPositionsMatrix | Position>maxPositionsMatrix);
        Velocity(IsOutside)=-Velocity(IsOutside);
        Method=1;
    end 
   
    % Apply Position Limits Boundary Control
    Position=update(Position,minPositionsMatrix,maxPositionsMatrix,BRM);
    
    % Evaluation
    FM_pop=Position;
    evals=evals+nPop;
    for i=1:nPop
        [xOpt(i)] = decoding_sol(id,FM_pop(i,:),Data,xOpt(i));
        [S_val_temp(i,1),xOpt(i)] = Fit_Evaluation(xOpt(i),Data,id);
        %[FM_pop(i,:)] = encoding_repaired_sol(xOpt(i),FM_pop(i,:),id,Data);
        %[S_val_temp, ~]=feval(fnc,FM_ui,caseStudyData, otherParameters);
    end
    
    Position=FM_pop;
    Cost=S_val_temp;
    
    % Update Personal Best
    ind=Cost<LocalBest.Cost; 
    LocalBest.Cost(ind)=Cost(ind);
    LocalBest.Position(ind,:)=Position(ind,:);
    % Update Global Best
    if min(LocalBest.Cost)<GlobalBest(it-1).Cost
        [~,idk]=min(LocalBest.Cost);
        GlobalBest(it).Position=LocalBest.Position(idk,:);
        GlobalBest(it).Cost=LocalBest.Cost(idk,1);
        BestSol.Result(it,1)=GlobalBest(it).Cost;
        Best_xOpt=xOpt(idk);
    else
        GlobalBest(it).Position=GlobalBest(it-1).Position;
        GlobalBest(it).Cost=GlobalBest(it-1).Cost;
        BestSol.Result(it,1)=BestSol.Result(it-1,1);
        
    end
    BestCost(it)=GlobalBest(it).Cost;
    fitMaxVector(1,it) = BestCost(it);
%     if BestCost(it)<BestCost(it-1)
%         fprintf('Fitness value: %f\n',BestCost(it))
%         fprintf('Iteration: %d\n',it)
%     end
end
BestSol.Obf=GlobalBest(end).Cost; 
BestSol.Position=GlobalBest(end).Position;
%fprintf('Fitness value: %f\n',BestCost(it))
Fit_and_p=[fitMaxVector(1,it) 0]; %;p2;p3;p4]
%figure;
%plot(BestCost);
%figure;
%plot(pPSO);
end

%% Adition function
% VECTORIZED THE CODE INSTEAD OF USING FOR
function pop=genpop(a,b,lowMatrix,upMatrix)
newlowMatrix=lowMatrix;
newupMatrix=upMatrix;
newlowMatrix(lowMatrix>upMatrix)=upMatrix(lowMatrix>upMatrix);
newupMatrix(upMatrix<lowMatrix)=lowMatrix(upMatrix<lowMatrix);

pop=unifrnd(newlowMatrix,newupMatrix,a,b);
end

function Position=update(Position,minPositionsMatrix,maxPositionsMatrix,BRM)
switch BRM
    case 1 % max and min replace
        Position = max(Position,minPositionsMatrix);
        Position = min(Position,maxPositionsMatrix);
    case 2 %Random reinitialization
        IsOutside=[find(Position<minPositionsMatrix);find(Position>maxPositionsMatrix)];
        Position(IsOutside)=unifrnd(minPositionsMatrix(IsOutside),maxPositionsMatrix(IsOutside),length(IsOutside),1);
    case 3 %Bounce Back 
        [IsOutsidemin] = find(Position<minPositionsMatrix);
        Position(IsOutsidemin)=unifrnd(minPositionsMatrix(IsOutsidemin),minPositionsMatrix(IsOutsidemin)-(Position(IsOutsidemin)+minPositionsMatrix(IsOutsidemin).*-1),length(IsOutsidemin),1);
        [IsOutsidemax] = find(Position>maxPositionsMatrix);
        Position(IsOutsidemax)=unifrnd(maxPositionsMatrix(IsOutsidemax)-(Position(IsOutsidemax)-maxPositionsMatrix(IsOutsidemax)),maxPositionsMatrix(IsOutsidemax),length(IsOutsidemax),1);
end
end

function Position = OtherMethodSearch(GBestPosition,nPop,MaxIt,it,nVar,minPositionsMatrix,maxPositionsMatrix)
    %% based on Vortex Search Algorithm 
    %radius decrement process
    x = 0.1; %x = 0.1 for gammaincinv(x,a) function
    a = (MaxIt-it) / MaxIt; % a value for gammaincinv function
    ginv = (1/x)*gammaincinv(x,a); % compute the new ginv value
    r = ginv * ((maxPositionsMatrix - minPositionsMatrix) / 2); %decrease the radius
    Mu=GBestPosition; %center is always shifted to the best solution found so far
    %create candidate solutions within the circle by using gaussian
    %distribution with standard deviation r and center Mu
    C = r.*randn(nPop,nVar);
    Position = bsxfun(@plus, C, Mu); %candidate solutions
end







