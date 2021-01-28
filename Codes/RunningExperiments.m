clear variables
close all
clc
addpath('functions')

Nplayers=2; %Change this to 20 for 20 players
filename = ['Data_Paper_2020_Journal_' num2str(Nplayers) '_players'];
%filename = 'Data_Paper_2020_Journal_20_players';

load(filename);

%% function for create the vector of limits
[upperB,lowerB,id] = Create_limits_vector(Data);

for exp=1:6
    
    Select_Algorithm=exp;
    %1: Deterministic Cplex
    %2: DE algorithm (test algorithm)
    %3: PSO-LVS
    ResDB=struc([]);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Load MH parameters (e.g., get MH parameters from DEparameters.m file)
    switch Select_Algorithm
        case 1
            addpath('DEalg')
            algorithm='DE_rand'; %'The participants should include their algorithm here'
            DEparameters %Function defined by the participant
            No_solutions=deParameters.I_NP; %Notice that some algorithms are limited to one individual
            
        case 2
            addpath('PSOalg')
            algorithm='PSO_LVS'; %'The participants should include their algorithm here'
            psoParameters %Function defined by the participant
            No_solutions=PSOparameters.nPop; %Notice that some algorithms are limited to one individual
        case 3
            addpath('alg_HyDEDF')
            algorithm='HyDE_DF'; %'The participants should include their algorithm here'
            HyDEparameters %Function defined by the participant
            No_solutions=deParameters.I_NP; %Notice that some algorithms are limited to one individual
            deParameters.I_strategy=3;
            deParameters.I_strategyVersion=2;
        case 4
            addpath('alg_HyDEDF')
            algorithm='HyDE'; %'The participants should include their algorithm here'
            HyDEparameters %Function defined by the participant
            No_solutions=deParameters.I_NP; %Notice that some algorithms are limited to one individual
            deParameters.I_strategy=3;
            deParameters.I_strategyVersion=3;
        case 5
            addpath('alg_HyDEDF')
            algorithm='VS'; %'The participants should include their algorithm here'
            HyDEparameters %Function defined by the participant
            No_solutions=deParameters.I_NP; %Notice that some algorithms are limited to one individual
            deParameters.I_strategy=3;
            deParameters.I_strategyVersion=1;
        case 6
            addpath('alg_HyDEDF')
            algorithm='DE_best'; %'The participants should include their algorithm here'
            HyDEparameters %Function defined by the participant
            No_solutions=deParameters.I_NP; %Notice that some algorithms are limited to one individual
            deParameters.I_strategy=2;
            deParameters.I_strategyVersion=1;
            
            
        otherwise
            fprintf(1,'No algorithm selected\n');
    end
    
    fileResultsname=['Results\' algorithm 'NP_20_Nplayers' num2str(Nplayers) '.mat'];
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Load Data base
    noRuns=30; %this can be changed but final results should be based on 20 trials
    
    %% Label of the algorithm and the case study
    Tag.algorithm=algorithm;
    Select_Algorithm
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Call the MH for optimizationclear
    for iRuns=1:noRuns %Number of trails
        fprintf('Run: %d\n',iRuns)
        tOpt=tic;
        %rand('state',sum(noRuns*100*clock))% ensure stochastic indpt trials
        rand('state',    iRuns)
        switch Select_Algorithm
            case 1
                [ResDB(iRuns).Fit_and_p, ...
                    ResDB(iRuns).sol, ...
                    ResDB(iRuns).fitVector...
                    ResDB(iRuns).Best_xOpt]= ...
                    deopt_simple(deParameters,Data,lowerB,upperB,id,iRuns);
            case 2
                [ResDB(iRuns).Fit_and_p, ...
                    ResDB(iRuns).sol, ...
                    ResDB(iRuns).fitVector ...
                    ResDB(iRuns).Best_xOpt]= ...
                    PSO_LVS(PSOparameters,Data,lowerB,upperB,id,iRuns);
            case {3,4,5,6}
                [ResDB(iRuns).Fit_and_p, ...
                    ResDB(iRuns).sol, ...
                    ResDB(iRuns).fitVector ...
                    ResDB(iRuns).Best_xOpt]= ...
                    HyDE(deParameters,Data,lowerB,upperB,id,iRuns);
        end
        
      
        ResDB(iRuns).tOpt=toc(tOpt); % time of each trial
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Save the results and stats
        fprintf('Fitness value: %f\n',ResDB(end).Fit_and_p(1))
        %Save_results
    end
    %     Results(iRuns).Time=ResDB.tOpt;
    %    Results(iRuns).f=ResDB(iRuns).Fit_and_p;
    %    Results(iRuns).x=
    
    save(fileResultsname)
    
end

%tTotalTime=toc(tTotalTime); %Total time
%% End of MH Optimization

