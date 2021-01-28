% Author: Fernando Lezama, GECAD/ISEP 2019
% Description:	Initial parameters as used in SWEVO paper
% 
% Please cite the following work when using HyDE-DF
% * Lezama et. al: HyDE-DF: A novel self-adaptive version of differential evolution for numerical optimization. Swarm and evolutionary computation. 2019
% * Lezama et. al: Hybrid-adaptive differential evolution with decay function (HyDE-DF) applied to the 100-digit challenge competition on single objective numerical optimization. In Proceedings of the Genetic and Evolutionary Computation Conference Companion (GECCO '19). 2019 DOI: https://doi.org/10.1145/3319619.3326747
% * Lezama et. al: A New Hybrid-Adaptive Differential Evolution for a Smart Grid Application Under Uncertainty. In IEEE Congress on Evolutionary Computation (CEC '19) (pp. 1-8). IEEE. 2018


deParameters.I_bnd_constr = 1; %Using bound constraints /is possible to change direct in DE
% 1 repair to the lower or upper violated bound (why)
% 2 rand value in the allowed range
% 3 bounce back

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Settings used in WCCI paper
deParameters.I_strategy=3;
deParameters.I_strategyVersion=3;
deParameters.I_itermax= 5e3;
            
deParameters.I_NP=20;
deParameters.F_weight=0.5;
deParameters.F_CR=0.9;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%