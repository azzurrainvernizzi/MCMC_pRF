function [p_max,p_avg,varExpl,pstore,posterior_latent,posterior,loglikelihood,prior] = mcmc_burn_in (Y,var_u,posterior_latent,posterior,pstore,loglikelihood,prior,radius,r_min,n,how_beta,burn_in)
%-------------------------------------------------------------------------

% Variables definition
p_max = zeros(1,size(posterior_latent,1));
p_avg = zeros(1,size(posterior_latent,1));

% Functions definition
rho=@(l_rho,radius) radius*normcdf(l_rho,0,1);

theta=@(l_theta) 2*pi.*normcdf(l_theta,0,1)-pi;

beta=@(l_beta) exp(l_beta);

alpha=@(l_alpha)  2*pi.*normcdf(l_alpha,0,1)-pi; 

sigma=@(l_sigma,radius, r_min)(radius-r_min).*normcdf(l_sigma,0,1)+r_min;

mu_x0=@(l_rho,radius,l_theta) rho(l_rho,radius).*cos(theta(l_theta));

mu_y0=@(l_rho,radius,l_theta) rho(l_rho,radius).*sin(theta(l_theta));

%-------------------------------------------------------------------------
if burn_in==0
    % Burn-in - remove the first 1000 of elements on 10000 runs
    if size(posterior_latent,2)<n
        burn_in=1;
    else
        
        n_burn=(size(posterior_latent,2)/n);
        pstore=pstore(n_burn+1:end);
        posterior_latent= posterior_latent(:,n_burn+1:end);
        var_u=var_u(:,n_burn+1:end);
        posterior= posterior(:,n_burn+1:end);
        prior=prior(:,n_burn+1:end);
        loglikelihood=loglikelihood(:,n_burn+1:end);
        
        A=find(pstore==max(pstore));
        max_A=posterior_latent(:,A(1));
        b=mean(posterior_latent(:,A),2);
        max_position=A(1);
        varExpl = 1 - var_u(max_position)./var(Y);
        
        p_max(1)=mu_x0(max_A(1),radius,max_A(2));
        p_max(2)=mu_y0(max_A(1),radius,max_A(2));
        p_max(3)=sigma(max_A(3),radius,r_min);
        p_avg(1)=mu_x0(b(1),radius,b(2));
        p_avg(2)=mu_y0(b(1),radius,b(2));
        p_avg(3)=sigma(b(3),radius,r_min);
       
        if how_beta==0
            p_max(end)=max_A(end);
            p_avg(end)=b(end);
        elseif how_beta==1
            p_max(end)=beta(max_A(end));
            p_avg(end)=beta(b(end));
        end
        
    end
    
elseif burn_in==1
    
    A=find(pstore==max(pstore));
    max_A=posterior_latent(:,A(1));
    max_position=A(1);
    b=mean(posterior_latent(:,A),2);
    max_position=A(1);
    varExpl = 1 - var_u(max_position)./var(Y);
    
    p_max(1)=mu_x0(max_A(1),radius,max_A(2));
    p_max(2)=mu_y0(max_A(1),radius,max_A(2));
    p_max(3)=sigma(max_A(3),radius,r_min);
    p_avg(1)=mu_x0(b(1),radius,b(2));
    p_avg(2)=mu_y0(b(1),radius,b(2));
    p_avg(3)=sigma(b(3),radius,r_min);
    
    if how_beta==0
        p_max(end)=max_A(end);
        p_avg(end)=b(end);
    elseif how_beta==1
        p_max(end)=beta(max_A(end));
        p_max(end)=beta(b(end));
    end
    
end

end
