%% run MCMC pRF for BrainLife.io
close all; clear all;

%% ----- Reading Nifti stimulus ----- %%
stim = read_stimulus_nii('stimulus.nii');
tSeries = read_tSeries_nii('tSeries.nii');

%% ----- Loading variables ----- %%
% most likely to be changed once in BL.io 
variables = v_definition(tSeries);

%% ----- Computing pRF----- %%
[bayes_pRF_fit] = mcmc_pRF_main(stim,tSeries,variables);

%% ----- Storing pRF output----- %%
storing_mcmc_pRF(bayes_pRF_fit,'sub1');