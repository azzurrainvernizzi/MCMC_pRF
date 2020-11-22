function out= v_definition(tSeries);

%% Initialization of variables for MCMC pRF

f=rng;
out.r_min=0.5; % minimum position to be estimated
out.radius=10.5;%%7;  % stimulus radius
out.n_volumes=size(tSeries,1);%128;
out.TR=1.5;
out.nparam=4;
out.n_samples=17500;
out.how_beta=0; % 0=estimate beta using classical glm or 1=get a posterior distribution of beta
out.burn_in=0; % 0= apply burn-in, 1=no burn-in
%how_pRF=0; % 0=classical pRF; 1=elliptical; 2=DoG.
out.n=10;
out.X=linspace(-out.radius,out.radius,1000);
[x, y] = meshgrid(-14:0.28:14);%(-radius:0.14:radius); % define visual space grid
out.x=x;
out.y=y;
out.accepted= false(1,out.n_samples);

% init latent parameters
out.l_rho=0.5;
out.l_theta=1;
out.l_sigma=1;
out.l_beta=1;

end

