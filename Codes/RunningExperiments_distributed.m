clear variables
close all
clc
addpath('functions')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load Data base
Nplayers=2; %Change this to 20 for 20 players
filename = ['Data_Paper_2020_Journal_' num2str(Nplayers) '_players'];
load(filename);
%filename = 'Data_Paper_2020_Journal_20_players';

for exp=1:6
    for player=1:Nplayers
        DataP.Ppv_max=Data.Ppv_max(:,player);
        DataP.Load_Total=Data.Load_Total(:,player);
        DataP.Load_DR_appli=Data.Load_DR_appli(:,3*player-2:3*player);
        DataP.Prices=Data.Prices;
        DataP.Bat=Data.Bat(:,player);
        DataP.grid=Data.grid(:,player);
        DataP.n_periods=Data.n_periods;
        DataP.n_pv=1;
        DataP.n_stor=1;
        DataP.n_prosumers=1;
        DataP.n_DR_loads=3;
        DataP.n_DR_loads_per_prosumer=3;
        %%
        player
        %% function for create the vector of limits
        [upperB,lowerB,id] = Create_limits_vector(DataP);
        
        Select_Algorithm=exp;
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
        
        fileResultsname=['Results_distributed\' algorithm 'NP_20_Nplayers' num2str(Nplayers) '.mat'];
        noRuns=20; %this can be changed but final results should be based on 20 trials
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
                        deopt_simple(deParameters,DataP,lowerB,upperB,id,iRuns);
                case 2
                    [ResDB(iRuns).Fit_and_p, ...
                        ResDB(iRuns).sol, ...
                        ResDB(iRuns).fitVector ...
                        ResDB(iRuns).Best_xOpt]= ...
                        PSO_LVS(PSOparameters,DataP,lowerB,upperB,id,iRuns);
                case {3,4,5,6}
                    [ResDB(iRuns).Fit_and_p, ...
                        ResDB(iRuns).sol, ...
                        ResDB(iRuns).fitVector ...
                        ResDB(iRuns).Best_xOpt]= ...
                        HyDE(deParameters,DataP,lowerB,upperB,id,iRuns);
            end
            ResDB(iRuns).tOpt=toc(tOpt); % time of each trial
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% Save the results and stats
            %fprintf('Fitness value: %f\n',ResDB(end).Fit_and_p(1))
            %Save_results
        end
        %Save results players
        Player(player).ResDB=ResDB;
    end
    
    %% Put all the results together again
    for iRuns=1:noRuns %Number of trails
        temp1=0;
        temp2=[];
        temp2b=[];
        temp3=0;
        temp4=0;
        temp5=0;
        
        for player=1:Nplayers
            temp1=temp1+Player(player).ResDB(iRuns).Fit_and_p(1,1);
            temp2=[temp2 Player(player).ResDB(iRuns).sol(1,1:96)];
            temp2b=[temp2b Player(player).ResDB(iRuns).sol(1,97:end)];
            temp3=temp3+Player(player).ResDB(iRuns).fitVector(1,:);
            %temp4=temp4+Player(player).ResDB(iRuns).fitVector(1,:);
            temp5=temp5+Player(player).ResDB(iRuns).tOpt(1,1);
        end
        
        ResDB(iRuns).Fit_and_p(1,1)=temp1;
        ResDB(iRuns).sol=[temp2 temp2b];
        ResDB(iRuns).fitVector=temp3;
        %load(filename);
        [upperB,lowerB,id] = Create_limits_vector(Data);
        [xOpt] = Create_Results_str(Data);
        [xOpt] = decoding_sol(id,ResDB(iRuns).sol,Data,xOpt);
        [fit,xOpt] = Fit_Evaluation(xOpt,Data,id);
        ResDB(iRuns).Best_xOpt=xOpt;
        ResDB(iRuns).tOpt(1,1)=temp5;
    end
    
    save(fileResultsname)
end

%tTotalTime=toc(tTotalTime); %Total time
%% End of MH Optimization

