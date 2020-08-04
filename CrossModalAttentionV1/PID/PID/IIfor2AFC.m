function  [I_II, S_R_info, C_R_info, S_C_info, non_readout_sensory_info, internal_choice_info, ...
    S_C_info_from_unobserved_R] = IIfor2AFC( S,R,C,preprocess )

% Expects vectors of binars S and C that indicate stimulus and choice and
% matrix of neural responses R.

% intersection_information: Computes intersection information II of a perceptual discrimination
% dataset where the experimenter recorded, in each of n_trial trials, the stimulus S, some neural feature R, and the choice C.
% 
% The information-theoretic intersection information II is defined and described in Pica et al (2017), Advances in Neural
% Information Processing, 3687-3697, "Quantifying how much sensory information in
% a neural code is relevant for behavior".

% The definition of II also defines immediately the leftover information components:
% I(S:R)-II is sensory information that is not readout for behavior, "non_readout_sensory_info"
% I(R:C)-II is choice information that is not related to the stimulus, "internal_choice_info"
% I(S:C)-II is correspondence between stimulus and choice that is due to other neural responses 
% than the observed R, "S_C_info_from_unobserverd_R" 

% The presumptive statistical model for the joint distribution of S,R,C
% is parameterized by log p(S,R,C) ~ S*w_s.R + C*w_c.R + a_s*S + a_c*C + b + J*S*C + log p_emp(R)
% where p_emp(R) is the empirical marginal distribution of neural responses


end

