function [data_bayes_pRF] = mcmc_pRF_main(stimulus,tSeries_data,v)

v.n_samples=50;
for u=1%:size(tSeries_data,2)
	 
    Y=tSeries_data(:,u);
    Y=Y./max(Y);
    varY = var(Y);% compute variance of Y
    err = varY;% set error matrix (U) to variance of Y
    varBase = ones(size(Y,1),1);% make vector of all ones for variable baseline
    
    % make hrf for the convolution ( normal hrf)
    hrf = makeHRF(0:1.5:24);
    
    % reshape stimulus
    stimulus_f=stimulus(:,:,1:size(Y,1));
    [a b c]=size(stimulus_f);
    stimulus2d = reshape(stimulus_f,[a*b c])';
    stim_time = 1:1.5:size(stimulus2d,1)*1.5;
    
    
    if exist('posterior','var');
        [v.nparam,col]=size(posterior); % allow to run the script twice to just double the number of samples.
        posterior=cat(2,posterior,zeros(v.nparam,v.n_samples));
    else
        posterior=zeros(v.nparam,v.n_samples);%each row will contain one parameter each collumn one (accepted) simulation result
        col=1;
    end
    
    if exist('posterior_latent','var');
        [v.nparam,col]=size(posterior_latent); % allow to run the script twice to just double the number of samples.
        posterior_latent=cat(2,posterior_latent,zeros(v.nparam,v.n_samples));
    else
        posterior_latent=zeros(v.nparam,v.n_samples);%each row will contain one parameter each collumn one (accepted) simulation result
        col=1;
    end
    
    
    for kk=1:v.n_samples
        
        if kk<=v.n_samples/10
            proposal_width=6; % coarse fit
        elseif kk<=v.n_samples/5 &&  kk>v.n_samples/10
            proposal_width=2;
        elseif kk<=v.n_samples/2.5 &&  kk>v.n_samples/5
            proposal_width=1;
        elseif kk<=v.n_samples/2 &&  kk>v.n_samples/2.5
            proposal_width=0.8;
        else
            proposal_width=0.5; % fine fit
        end
        
        
        if v.how_beta==0
            
            l_rho_proposal=normrnd(v.l_rho, proposal_width);
            l_theta_proposal=normrnd(v.l_theta, proposal_width);
            l_sigma_proposal=normrnd(v.l_sigma, proposal_width);
            
            [p_current,xi_current,var_current,loglike_current,prior_current] = computing_mcmc (Y,v.x,v.y,v.radius,hrf,stimulus2d,v.l_rho,v.l_theta,v.l_sigma,v.r_min,v.l_beta,v.how_beta);
            [p_proposal,xi_proposal,var_proposal,loglike_proposal,prior_proposal] = computing_mcmc (Y,v.x,v.y,v.radius,hrf,stimulus2d,l_rho_proposal,l_theta_proposal,l_sigma_proposal,v.r_min,v.l_beta,v.how_beta);
            
            pvalue_current(kk)=p_current;
            pvalue_proposal(kk)=p_proposal;
            var_curr(kk)=var_current;
            var_pro(kk)=var_proposal;
            
            p_accept(kk)=exp(p_proposal-p_current);
            accept = rand(1);
            accepted(kk)=abs(accept)<p_accept(kk);
            
            if accepted(kk)==1
                
                v.l_rho = l_rho_proposal;
                v.l_theta=l_theta_proposal;
                v.l_sigma=l_sigma_proposal;
                
                pstore(kk)=p_proposal;
                var_u(kk)=var_proposal;
                loglikelihood(kk)=loglike_proposal;
                prior(kk)=prior_proposal;
                posterior_latent(:,col+kk-1)=[v.l_rho;v.l_theta;v.l_sigma;v.l_beta];
                posterior(:,col+kk-1)=xi_proposal;
                
            else
                pstore(kk)=p_current;
                var_u(kk)=var_current;
                loglikelihood(kk)=loglike_current;
                prior(kk)=prior_current;
                posterior_latent(:,col+kk-1)=[v.l_rho;v.l_theta;v.l_sigma;v.l_beta];
                posterior(:,col+kk-1)=xi_current;
                
            end
            
        elseif how_beta==1
            
            l_rho_proposal=normrnd(v.l_rho, proposal_width);
            l_theta_proposal=normrnd(v.l_theta, proposal_width);
            l_sigma_proposal=normrnd(v.l_sigma, proposal_width);
            l_beta_proposal=normrnd(v.l_beta, proposal_width);
            
            [p_current,xi_current,var_current,loglike_current,prior_current] = computing_mcmc (Y,v.x,v.y,v.radius,hrf,stimulus2d,v.l_rho,v.l_theta,v.l_sigma,v.r_min,v.l_beta,v.how_beta);
            [p_proposal,xi_proposal,var_proposal,loglike_proposal,prior_proposal] = computing_mcmc (Y,v.x,v.y,v.radius,hrf,stimulus2d,l_rho_proposal,l_theta_proposal,l_sigma_proposal,v.r_min,l_beta_proposal,v.how_beta);
            
            pvalue_current(kk)=p_current;
            pvalue_proposal(kk)=p_proposal;
            var_curr(kk)=var_current;
            var_pro(kk)=var_proposal;
            
            p_accept(kk)=exp(p_proposal-p_current);
            accept = rand(1);
            accepted(kk)=abs(accept)<p_accept(kk);
            
            if accepted(kk)==1
                
                l_rho = l_rho_proposal;
                l_theta=l_theta_proposal;
                l_sigma=l_sigma_proposal;
                l_beta=l_beta_proposal;
                
                pstore(kk)=p_proposal;
                var_u(kk)=var_proposal;
                loglikelihood(kk)=loglike_proposal;
                prior(kk)=prior_proposal;
                posterior_latent(:,col+kk-1)=[v.l_rho;v.l_theta;v.l_sigma;v.l_beta];
                posterior(:,col+kk-1)=xi_proposal;
                
            else
                pstore(kk)=p_current;
                var_u(kk)=var_current;
                loglikelihood(kk)=loglike_current;
                prior(kk)=prior_current;
                posterior_latent(:,col+kk-1)=[v.l_rho;v.l_theta;v.l_sigma;v.l_beta];
                posterior(:,col+kk-1)=xi_current;
                
            end
        end
        
    end
                                                                       
    [p_max{u},p_avg{u},varExpl{u},pstore_b{u},posterior_latent_b{u},posterior_b{u},loglikelihood_b{u},prior_b{u}] = mcmc_burn_in(Y,var_u,posterior_latent,posterior,pstore,loglikelihood,prior,v.radius,v.r_min,v.n,v.how_beta,v.burn_in);

    data_bayes_pRF = struct('p_max', p_max,'varExpl', varExpl, 'pstore_b', pstore_b,'posterior_latent_b',posterior_latent_b,'loglikelihood_b',loglikelihood_b,'prior_b',prior_b);
    
    clear var_u; clear posterior_latent; clear posterior; clear pstore;
    clear prior; clear loglikelihood;
    
end

end
