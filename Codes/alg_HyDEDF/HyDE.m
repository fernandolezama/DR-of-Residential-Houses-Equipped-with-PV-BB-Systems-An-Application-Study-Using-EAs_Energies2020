% Function: function [Fit_and_p,FVr_bestmemit, fitMaxVector] =HyDE(deParameters,otherParameters,low_habitat_limit,up_habitat_limit,initialSolution)
% Author: Fernando Lezama, GECAD/ISEP 2019 (Contact ing.flezama@gmail.com)
% Description:	Minimization of a user-supplied function with respect to x(1:I_D), using Hybrid-adaptive differential evolution with decay function (HyDE-DF). 
% For this algorithm we download DE source code from Rainer Storn, Ken Price, Arnold Neumaier, Jim Van Zandt (http://www1.icsi.berkeley.edu/~storn/code.html#matl) and modified it.
% Due to the vectorized expressions deopt executes  fairly fast in MATLAB's interpreter environment.
% 
% Please cite the following work when using HyDE-DF
% * Lezama et. al: HyDE-DF: A novel self-adaptive version of differential evolution for numerical optimization. Swarm and evolutionary computation. 2019
% * Lezama et. al: Hybrid-adaptive differential evolution with decay function (HyDE-DF) applied to the 100-digit challenge competition on single objective numerical optimization. In Proceedings of the Genetic and Evolutionary Computation Conference Companion (GECCO '19). 2019 DOI: https://doi.org/10.1145/3319619.3326747
% * Lezama et. al: A New Hybrid-Adaptive Differential Evolution for a Smart Grid Application Under Uncertainty. In IEEE Congress on Evolutionary Computation (CEC '19) (pp. 1-8). IEEE. 2018
% 
% ---------------------------------------------------------------------------------------------------
% Parameters:	deParameters 	(I)    	Struct with DE required parameters.
% 		otherParameters (I) 	Struct with Problem data information.
% 		low_habitat_limit (I)	Vector of variables lower bounds
% 		up_habitat_limit (I)	Vector of variables upper bounds
% 		initialSolution (I)	Taylored Initial Solutions (this is optional)
% 
% %-------------------------------members of deParameters structure----------------------------------------------------
% deParameters.I_bnd_constr:	Boundary constraint method	1--> Variable limit
%                                                           2--> Random
%                                                           3--> Bounce-back
% deParameters.I_strategy: 	To select other DE variants	1--> DE/rand/1
%                                                       2--> DE/local-to-best/1 
%                                                       3--> Vortex, HyDE, or HyDE-DF
% deParameters.I_strategyVersion: if deParameters.I_strategy==3     1--> Vortex
%                                                                   2--> HyDE-DF ******Our algorithm**********
%                                                                   3--> HyDE
% (**Select deParameters.I_strategy==3 and deParameters.I_strategyVersion==2 to use HyDE-DF**)
% 
% deParameters.I_itermax:	Max iteration limit
% deParameters.I_NP:	Size of population
% deParameters.F_weight: 	initial mutation factor
% deParameters.F_CR:	initial crossover factor
% 
% %-------------------------------members of otherParameters structure----------------------------------------------------
% otherParameters.objfun:		objective function name (e.g., 'stepint')
% otherParameters.dim:		Dimension of the function (e.g., 5)
% otherParameters.lowerlimit: 	variable lower bound 
% otherParameters.upperlimit: 	variable upper bound 
% 
% -------------------------------------------------------------------------------------                                       
% Return value:    
% Fit_and_p (O)	 	Best objective function value
% FVr_bestmemit (O)    	Best solution vector.
% fitMaxVector (O)	Vector with objective function value at each iteration
% 
% 
% NOTE:
% This program is free software; Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the “Software”),
% to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software,
% and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
% 
% Distribution of the original version of HyDE-DF can be done as long as the original contributors’ names for are endorsed. Derived works from this software are possible but endorsements
% of the original contributors's name without specific permissions are not.
% 
% The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

%function [Fit_and_p,FVr_bestmemit, fitMaxVector] = ...
%    HyDE(deParameters,otherParameters,low_habitat_limit,up_habitat_limit,initialSolution)

function [Fit_and_p,FVr_bestmemit, fitMaxVector,Best_xOpt] = ...
     HyDE(deParameters,Data,low_habitat_limit,up_habitat_limit,id,iRuns)
%biobj=otherParameters.biobj;
%idf=otherParameters.idf;


%-----This is just for notational convenience and to keep the code uncluttered.--------
I_NP         = deParameters.I_NP;
F_weight     = deParameters.F_weight;
F_CR         = deParameters.F_CR;
I_D          = numel(up_habitat_limit); %Number of variables or dimension
deParameters.nVariables=I_D;
FVr_minbound = low_habitat_limit;
FVr_maxbound = up_habitat_limit;
I_itermax    = deParameters.I_itermax;

%Repair boundary method employed
BRM=deParameters.I_bnd_constr; %1: bring the value to bound violated
                               %2: repair in the allowed range

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
I_strategy   = deParameters.I_strategy; %important variable
%fnc= otherParameters.objfun;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%-----Check input variables---------------------------------------------
if (I_NP < 5)
   I_NP=5;
   fprintf(1,' I_NP increased to minimal value 5\n');
end
if ((F_CR < 0) || (F_CR > 1))
   F_CR=0.5;
   fprintf(1,'F_CR should be from interval [0,1]; set to default value 0.5\n');
end
if (I_itermax <= 0)
   I_itermax = 200;
   fprintf(1,'I_itermax should be > 0; set to default value 200\n');
end

%-----Initialize population and some arrays-------------------------------
%FM_pop = zeros(I_NP,I_D); %initialize FM_pop to gain speed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pre-allocation of loop variables
fitMaxVector = nan(1,I_itermax+1);
% limit iterations by threshold
gen = 1; %iterations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%----FM_pop is a matrix of size I_NPx(I_D+1). It will be initialized------
%----with random values between the min and max values of the-------------
%----parameters-----------------------------------------------------------
% FLC modification - vectorization
minPositionsMatrix=repmat(FVr_minbound,I_NP,1);
maxPositionsMatrix=repmat(FVr_maxbound,I_NP,1);
deParameters.minPositionsMatrix=minPositionsMatrix;
deParameters.maxPositionsMatrix=maxPositionsMatrix;

% generate initial population.
 rand('state',    iRuns)
 
FM_pop=genpop(I_NP,I_D,minPositionsMatrix,maxPositionsMatrix);
% if nargin>5
%     noInitialSolutions = size(initialSolution,1);
%     FM_pop(1:noInitialSolutions,:)=initialSolution;
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%------Evaluate the best member after initialization----------------------
% Modified by FLC

% for i=1:I_NP
% [S_val_R(i,1),g_R,h_R]=feval(fnc,FM_pop(i,:),biobj,idf);
% [S_val(i,1),g_error,h_error]=contraint_handling(S_val_R(i,1),g_R,h_R);
% end
for i=1:I_NP
    [xOpt(i)] = Create_Results_str(Data);
    [xOpt(i)] = decoding_sol(id,FM_pop(i,:),Data,xOpt(i));
    [S_val(i,1),xOpt(i)] = Fit_Evaluation(xOpt(i),Data,id);
%    [FM_pop(i,:)] = encoding_repaired_sol(xOpt(i),FM_pop(i,:),id,Data);
end

[S_bestval,I_best_index] = min(S_val); % This mean that the best individual correspond to the best worst performance
Best_xOpt=xOpt(I_best_index);
FVr_bestmemit = FM_pop(I_best_index,:); % best member of current iteration
fitMaxVector(1,gen) = S_bestval;

% The user can decide to save the mean, best, or any other value here

%------DE-Minimization---------------------------------------------
%------FM_popold is the population which has to compete. It is--------
%------static through one iteration. FM_pop is the newly--------------
%------emerging population.----------------------------------------
FVr_rot  = (0:1:I_NP-1);               % rotating index array (size I_NP)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%HYDE
if deParameters.I_strategy==3
        F_weight_old=repmat(F_weight,I_NP,3);
        F_weight= F_weight_old;
        F_CR_old=repmat(F_CR,I_NP,1);
        F_CR=F_CR_old;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
I_strategyVersion=deParameters.I_strategyVersion;

while gen<=I_itermax  %%&&  fitIterationGap >= threshold
    %a = itr / MaxItr; % a value for gammaincinv function
    other.a=(I_itermax-gen)/I_itermax;
    other.lowerlimit=FVr_minbound; %lower limit of the problem
    other.upperlimit=FVr_maxbound; %upper limit of the problem
    
     if deParameters.I_strategy==3
                value_R=rand(I_NP,3);
                ind1=value_R<0.1;
                ind2=rand(I_NP,1)<0.1;
                F_weight(ind1)=0.1+rand(sum(sum(ind1)),1)*0.9;
                F_weight(~ind1)=F_weight_old(~ind1);
                F_CR(ind2)=rand(sum(ind2),1);
                F_CR(~ind2)=F_CR_old(~ind2);
     end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [FM_ui,FM_base,~]=generate_trial(I_strategy,F_weight, F_CR, FM_pop, FVr_bestmemit,I_NP, I_D, FVr_rot,I_strategyVersion,other);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 

    %% Boundary Control
    FM_ui=update(FM_ui,minPositionsMatrix,maxPositionsMatrix,BRM,FM_base);

      %Evaluation of new Pop
    for i=1:I_NP
        [xOpt(i)] = decoding_sol(id,FM_ui(i,:),Data,xOpt(i));
        [S_val_temp(i,1),xOpt(i)] = Fit_Evaluation(xOpt(i),Data,id);
        %[FM_ui(i,:)] = encoding_repaired_sol(xOpt(i),FM_ui(i,:),id,Data);
        %[S_val_temp, ~]=feval(fnc,FM_ui,caseStudyData, otherParameters);
    end
    
    %%%%%%%%%%%%%%%%
    
     %% Elitist Selection
    ind=find(S_val_temp<S_val);
    S_val(ind)=S_val_temp(ind);
    FM_pop(ind,:)=FM_ui(ind,:);
    
    
  
   %% update best results
    [S_bestval,I_best_index] = min(S_val);
    Best_xOpt=xOpt(I_best_index);
    FVr_bestmemit = FM_pop(I_best_index,:); % best member of current iteration
    % store fitness evolution and obj fun evolution as well
%     if gen>1 & fitMaxVector(1,gen)<fitMaxVector(1,gen-1)
%         fprintf('Fitness value: %f\n',fitMaxVector(1,gen) )
%         fprintf('Generation: %d\n',gen)
%     end
 
    gen=gen+1;
    fitMaxVector(1,gen) = S_bestval;
    
    
    
     if deParameters.I_strategy==3 %jDE
        F_weight_old(ind,:)=F_weight(ind,:);
        F_CR_old(ind)=F_CR(ind);
     end

end %---end while ((I_iter < I_itermax) ...
%p1=sum(Best_otherInfo.penSlackBusFinal);
Fit_and_p=fitMaxVector(1,gen); %;p2;p3;p4]


 
% VECTORIZED THE CODE INSTEAD OF USING FOR
function pop=genpop(a,b,lowMatrix,upMatrix)
pop=unifrnd(lowMatrix,upMatrix,a,b);

% VECTORIZED THE CODE INSTEAD OF USING FOR
function p=update(p,lowMatrix,upMatrix,BRM,FM_base)
switch BRM
    case 1 %Our method
        %[popsize,dim]=size(p);
        [idx] = find(p<lowMatrix);
        p(idx)=lowMatrix(idx);
        [idx] = find(p>upMatrix);
        p(idx)=upMatrix(idx);
    case 2 %Random reinitialization
        [idx] = [find(p<lowMatrix);find(p>upMatrix)];
        replace=unifrnd(lowMatrix(idx),upMatrix(idx),length(idx),1);
        p(idx)=replace;
    case 3 %Bounce Back
      [idx] = find(p<lowMatrix);
      p(idx)=unifrnd(lowMatrix(idx),FM_base(idx),length(idx),1);
        [idx] = find(p>upMatrix);
      p(idx)=unifrnd(FM_base(idx), upMatrix(idx),length(idx),1);
end

function [FM_ui,FM_base,msg]=generate_trial(method,F_weight, F_CR, FM_pop, FVr_bestmemit,I_NP,I_D,FVr_rot,I_strategyVersion,other)
    FM_popold = FM_pop;                  % save the old population
    FVr_ind = randperm(4);               % index pointer array
    FVr_a1  = randperm(I_NP);                   % shuffle locations of vectors
    FVr_rt  = rem(FVr_rot+FVr_ind(1),I_NP);     % rotate indices by ind(1) positions
    FVr_a2  = FVr_a1(FVr_rt+1);                 % rotate vector locations
    FVr_rt  = rem(FVr_rot+FVr_ind(2),I_NP);
    FVr_a3  = FVr_a2(FVr_rt+1);                
    FM_pm1 = FM_popold(FVr_a1,:);             % shuffled population 1
    FM_pm2 = FM_popold(FVr_a2,:);             % shuffled population 2
    FM_pm3 = FM_popold(FVr_a3,:);             % shuffled population 3
  
    if length(F_CR)==1  %Meaning the same F_CR for all individuals
        FM_mui = rand(I_NP,I_D) < F_CR;  % all random numbers < F_CR are 1, 0 otherwise
        FM_mpo = FM_mui < 0.5;    % inverse mask to FM_mui
    else %Meaning a different F_CR for each individual
        FM_mui = rand(I_NP,I_D) < repmat(F_CR,1,I_D);  % all random numbers < F_CR are 1, 0 otherwise
        FM_mpo = FM_mui < 0.5;    % inverse mask to FM_mui
    end

    switch method %different implementations available
        case 1 %DE/rand1
            FM_ui = FM_pm3 + F_weight*(FM_pm1 - FM_pm2);   % differential variation
            FM_ui = FM_popold.*FM_mpo + FM_ui.*FM_mui;     % crossover
            FM_base = FM_pm3;
            msg=' DE/rand/bin';
        case 2
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
            %VEC by FLC
            FM_bm=repmat(FVr_bestmemit,I_NP,1);
            FM_ui = FM_popold + F_weight*(FM_bm-FM_popold) + F_weight*(FM_pm1 - FM_pm2);
            FM_ui = FM_popold.*FM_mpo + FM_ui.*FM_mui;
            FM_base = FM_bm;
            msg=' DE/current-to-best/1';
        case 3 %jDEPerturbated_v3 v4... v7
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
            FM_bm=repmat(FVr_bestmemit,I_NP,1);
            if length(F_weight)==1  %Meaning the same F_weight for all individuals
                FM_ui = FM_popold + F_weight*(FM_bm-FM_popold) + F_weight*(FM_pm1 - FM_pm2);
            else
                if  I_strategyVersion==1 %Emulate Vortex Algorithm
                        a=other.a;
                        ginv = (1/0.1)*gammaincinv(0.1,a); % compute the new ginv value
                        r = ginv * ((other.upperlimit - other.lowerlimit) / 2); %decrease the radius
                        C = r.*randn(I_NP,I_D);
                        FM_ui = bsxfun(@plus, C, FM_bm(1,:));
                end
                
                if  I_strategyVersion==2 %HyDE-DF
                    a=other.a; %Linear decrease
                   % ginv=exp((1-(1/a^2))); %Exponential decreasing funtion
                    
                     ginv = (1/0.1)*gammaincinv(0.1,a); %Testing vortex effect
                    FM_ui = FM_popold + repmat(F_weight(:,3),1,I_D).*(FM_pm1 - FM_pm2)  + ginv*(repmat(F_weight(:,1),1,I_D).*(FM_bm.*(repmat(F_weight(:,2),1,I_D)+randn(I_NP,I_D))-FM_popold));   % differential variation
                end
                
                if  I_strategyVersion==3 %HyDE
                    FM_ui = FM_popold + repmat(F_weight(:,1),1,I_D).*(FM_bm.*(repmat(F_weight(:,2),1,I_D)+randn(I_NP,I_D))-FM_popold) + repmat(F_weight(:,3),1,I_D).*(FM_pm1 - FM_pm2);
                end
            end
            
            FM_ui = FM_popold.*FM_mpo + FM_ui.*FM_mui;
            FM_base = FM_bm;
            msg=' HyDE/current-to-best/1';   
    end
return

