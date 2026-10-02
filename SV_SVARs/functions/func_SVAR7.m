function [store_h,store_Sigh,store_B0, store_alpha, Sig_mean, Sig_quant, store_phi, store_tau, store_lambda, store_h0] = func_SVAR7(data, HC, nsim, burnin, thin, p, var_id, varB0)
addpath('functions/utility')

% var_id should be a vector of variable indices, e.g., [1 2 3] for the first three variables


effsamp=(nsim)/thin;

y0limit = 4*p; % These observations are used to get the variances for the Minnesota prior
Y0 = data(1:y0limit,var_id);  % save the first y0limit obs as the initial conditions
Y = data((y0limit+1):end,var_id);

[T,n] = size(Y);
tmpY = [Y0(end-p+1:end,:); Y];
X = zeros(T,n*p);
for ii=1:p
    X(:,(ii-1)*n+1:ii*n) = tmpY(p-ii+1:end-ii,:);
end
X = [ones(T,1) X];


k = 1+n*p;
k_alpha = n*k;
sig2 = get_resid_var(Y0,Y);
[C,idx_kappa1,idx_kappa2] = get_C(n,p,sig2);




    % priors
Hyper.nuh = 5*ones(n,1); 
Hyper.Sh = .01*ones(n,1).*(Hyper.nuh-1);
Hyper.phi0 = .95*ones(n,1); 
Hyper.Vphi = .000000001^2*ones(n,1);
Hyper.B0 = eye(n); 
Hyper.VB0 = varB0*ones(n);


    % initialization for storeage
store_kappa = zeros(effsamp,2);
store_alpha = zeros(effsamp,k_alpha);
store_B0 = zeros(effsamp,n^2);
store_h = zeros(effsamp,T,n);
count_phi = zeros(n,1);
store_phi = zeros(effsamp,n);

store_Sigh = zeros(effsamp,T,n);
store_h0 = zeros(effsamp,n);
store_tau = zeros(effsamp, n);
store_lambda = zeros(effsamp, n, T);


start_time = clock;

%% MCMC starts here
randn('seed',sum(clock*100)); rand('seed',sum(clock*1000));
disp('Starting MCMC for BSVAR.... ');

    % initialize the Markov chain
kappa = [.1,.1,NaN,100]; % kappa(1): own lag; kappa(2): other lag; kappa(4): intercepts
A = (X'*X + .01*speye(k))\(X'*Y);
U = Y-X*A;
alpha = A(:);
Sig_hat = U'*U/T;
B0 = diag(1./sqrt(diag(Sig_hat)));
h = repmat(log(diag(Sig_hat))',T,1);
phi = Hyper.phi0 + sqrt(Hyper.Vphi).*randn(n,1);
E = zeros(T,n);
z_psi1 = 1./gamrnd(.5,1,n*p,1);
z_psi2 = 1./gamrnd(.5,1,(n-1)*n*p,1);
z_kappa = 1./gamrnd(.5,1,2,1);
psi_kappa1 = 1./gamrnd(.5,z_psi1);
psi_kappa2 = 1./gamrnd(.5,z_psi2);
Psi = ones(k*n,1);
Psi(idx_kappa1) = psi_kappa1; Psi(idx_kappa2) = psi_kappa2;

tauh=ones(n,1)*0.1;
lambdah=ones(n,T)*10;
Sigh=(lambdah.*tauh)';


for isim = 1:nsim + burnin    
        % sample B0
    U = Y-X*A;
    for ii=1:n
        b0i = Hyper.B0(ii,:); % does not matter row or col, prior is identity
        UOhiU = (U'*sparse(1:T,1:T,exp(-h(:,ii))))*U;
        iVbi = sparse(1:n,1:n,1./Hyper.VB0(ii,:)); % diagonal sparse mat      
        Kbi = iVbi + UOhiU;
        bi_hat = Kbi\(iVbi*b0i'); % transpose b0i so it's col
        Ci = chol((Kbi),'lower')/sqrt(T);
        B0_minusi = B0([1:ii-1 ii+1:end],:)'; % selecting rows and transposing
        ortB0_minusi = null(B0_minusi');
        
        V = zeros(n,n); zeta = zeros(n,1);        
        for jj=1:n
            if jj==1                
                v1 = Ci\ortB0_minusi; 
                v1 = v1/norm(v1); %divide by norm
                V = [v1 null(v1')];
                zeta1_hat = bi_hat'*(Ci*v1);
                zeta(1) = anormrnd(zeta1_hat,1/T); % absolute normal
            else
                zetaj_hat = bi_hat'*(Ci*V(:,jj));
                zeta(jj) = zetaj_hat + 1/sqrt(T)*randn; %uncentred param
            end
        end
        bi = (Ci')\sum(V.*repmat(zeta',n,1),2); %all rows of repmat are the same (zetas), so each col v_j times zeta_j
        B0(ii,:) = bi*sign(bi(ii)); % fix the sign of the i-th element to be positive
    end

    
    
        % sample alpha    
    np = length(idx_kappa1);
    CPsi = C.*Psi;
    k_beta = length(CPsi);
    Vbeta = zeros(k_beta,1);
    Vbeta(1:np+1:end) = kappa(4)*sig2;          % intercepts
    Vbeta(idx_kappa1) = kappa(1)*CPsi(idx_kappa1); % own lags
    Vbeta(idx_kappa2) = kappa(2)*CPsi(idx_kappa2); % other lags

    Lambda = vec(exp(h/2));
    for ii=1:n
        A(:,ii) = 0;
        yi = vec((Y-X*A)*B0')./Lambda;
        Wi = kron(B0(:,ii),X)./Lambda;     
        iValphai = sparse(1:k,1:k,1./Vbeta((ii-1)*k+1:ii*k));
        Kalphai = iValphai + Wi'*Wi;
        try
            CKalphai = chol(Kalphai, 'lower');
        catch
                    
            % Add a small ridge until it succeeds
            jitter = 1e-10 ;
            if jitter == 0, jitter = 1e-10; end
        
            success = false;
            while ~success
                try
                    CKalphai = chol(Kalphai + jitter*eye(size(Kalphai)), 'lower');
                    success = true;
                catch
                    jitter = jitter * 10;   % increase if still failing
                end
            end
         end

        alphai_hat = (CKalphai')\(CKalphai\(Wi'*yi));        
        alphai = alphai_hat + CKalphai'\randn(k,1); %non-centered param
        A(:,ii) = alphai;         
    end  
    alpha = A(:); %stack columns

    % set h0 to zero
    h0= zeros(1,n)';

           % sample h
    E = (Y - X*A)*B0';
    for ii=1:n
        ystar = log(E(:,ii).^2 + 0.000000000001);
        h(:,ii) = sample_SV2(ystar,h(:,ii),h0(ii),phi(ii),Sigh(:,ii));
    end 
        
    
        % sample Sigh
    if HC==1
        e = [(h(1,:)-h0').*sqrt(1-phi.^2)'; h(2:end,:)-phi'.*h(1:end-1,:)]';

        v=1./gamrnd(1,1./(1+1./lambdah));
        lambdah=1./gamrnd(1, 1./(1./v + (e.^2)./(2.*tauh)));

        xi=1./gamrnd(1,1./(1+1./tauh));
        tauh=1./gamrnd((T+1)/2, 1./(1./xi +0.5*sum(e.^2./lambdah,2)));

        Sigh=(lambdah.*tauh)'+1e-10;
    else    
        e = [(h(1,:)-h0').*sqrt(1-phi.^2)'; h(2:end,:)-phi'.*h(1:end-1,:)];
        Sigh = 1./gamrnd(Hyper.nuh+T/2, 1./(Hyper.Sh + sum(e.^2)'/2));
        Sigh=repmat(Sigh,1,T)';
    end 


        % sample phi
    Xphi = h(1:end-1,:);
    zphi = h(2:end,:);
    Kphi = 1./Hyper.Vphi + sum( Xphi.^2 ./ Sigh(1:T-1,:) )';
    phi_hat = (Hyper.phi0./Hyper.Vphi + sum( (Xphi .* zphi) ./ sqrt(Sigh(2:T,:)) )' )./Kphi;
    phic = phi_hat + 1./sqrt(Kphi).*randn(n,1);
    flag_phi = zeros(n,1);
    for ii = 1:n
        g_phi = @(x) .5*log(1-x^2) -.5*(1-x^2)/Sigh(1,ii)*h(1,ii)^2;
       if abs(phic(ii))<.99
            alpMH = exp(g_phi(phic(ii))-g_phi(phi(ii)));
            if alpMH>rand
                phi(ii) = phic(ii);
                flag_phi(ii) = 1;
            end
       end 
    end
    count_phi = count_phi+flag_phi;


        % sample psi
    tmpv1 = 1./z_psi1 + alpha(idx_kappa1).^2./(2*C(idx_kappa1)*kappa(1));
    tmpv2 = 1./z_psi2 + alpha(idx_kappa2).^2./(2*C(idx_kappa2)*kappa(2));
    psi_kappa1 = 1./gamrnd(1,1./tmpv1);
    psi_kappa2 = 1./gamrnd(1,1./tmpv2);
    Psi(idx_kappa1) = psi_kappa1; 
    Psi(idx_kappa2) = psi_kappa2;
    
        % sample z_psi
    z_psi1 = 1./gamrnd(1,1./(1+1./psi_kappa1));
    z_psi2 = 1./gamrnd(1,1./(1+1./psi_kappa2));
        
        % sample kappa1 and kappa2
    tmpc1 = 1/z_kappa(1) + sum(alpha(idx_kappa1).^2./(2*psi_kappa1.*C(idx_kappa1)));
    tmpc2 = 1/z_kappa(2) + sum(alpha(idx_kappa2).^2./(2*psi_kappa2.*C(idx_kappa2)));    
    kappa(1) = 1./gamrnd((n*p+1)/2,1./tmpc1);   
    kappa(2) = 1./gamrnd(((n-1)*n*p+1)/2,1./tmpc2);
    
        % sample z_kappa
    z_kappa = 1./gamrnd(1,1./(1+1./kappa(1:2)));    

        % store
    if isim > burnin && mod(isim, thin)== 0
        isave = (isim-burnin)/thin;  
        store_kappa(isave,:) = kappa(1:2);
        store_alpha(isave,:) =  alpha';
        store_B0(isave,:) =  reshape(B0',n^2,1); % stacked by rows
        store_h(isave,:,:) = h; 
        store_Sigh(isave,:,:) = Sigh;
        store_h0(isave,:) = h0';  
        store_phi(isave,:) = phi';
        store_tau(isave,:) = tauh';
        store_lambda(isave,:,:) = lambdah;
   
        
    end
    
    if (mod(isim, 10000) ==0)
        disp([num2str(isim) ' loops... ' ])
    end     
end

disp( ['MCMC takes '  num2str( etime( clock, start_time) ) ' seconds' ] );
disp(' ' );

%% Sig

store_Sig = zeros(effsamp,T,n,n);
for isim = 1:effsamp
    B0 = reshape(store_B0(isim,:),n,n)'; % stacked by rows
    h_squeezed = squeeze(store_h(isim,:,:));
    store_Sig(isim,:,:,:) = construct_Sigt(h_squeezed,B0);
end
Sig_mean = mean(store_Sig);
Sig_quant = quantile(store_Sig,[0.05, 0.16, 0.50, 0.84, 0.95]);

