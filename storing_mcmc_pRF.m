function storing_mcmc_pRF(bayes_pRF,subject)

matfile=['bayes_pRF_fit_' num2str(subject) '.mat'];    
save (matfile,'bayes_pRF');

fprintf('storing %s ...\n',subject);

end