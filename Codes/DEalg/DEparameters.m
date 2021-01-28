% Author:           Rainer Storn, Ken Price, Arnold Neumaier, Jim Van Zandt
% Modified by FLC \GECAD 04/winter/2017

deParameters.I_NP= 20; % population in DE
deParameters.F_weight= 0.5; %Mutation factor
deParameters.F_CR= 0.9; %Recombination constant
deParameters.I_itermax= 4000; % number of max iterations/gen
deParameters.I_strategy   = 1; %DE strategy

deParameters.I_bnd_constr = 3; %Using bound constraints 
% 1 repair to the lower or upper violated bound 
% 2 rand value in the allowed range
% 3 bounce back

%
