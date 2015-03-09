clear all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   December 5, 2014
%   In this file we intend to evaluate hom PGAS performs on factorial
%   finite state machines with multiple parallel chains, bounded memory and
%   complex Gaussian observations. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%    Defining the model parameters 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%  We use the following transition model
%       Pr[ x_t(m) = i | x_{t-1}(m) = 0 } = a_m if i=0 and (1-am)/Q if
%       i=1,2,...,Q.
%       Pr[ x_t(m) = i | x_{t-1}(m) > 0 } = b_m if i=0 and (1-bm)/Q if
%       i=1,2,...,Q.
%       It is further assumed that x_0(m)=0 for all m, i.e., for m = 1, 2,
%       ..., Nt. 
%       Note: in the current implementation, we use x_t to denote state
%       variables that take values in the range 0, 1, ..., Q. The actual
%       values of the state is z_{tm} = C(x_t(m)), where C is a vector that
%       describes the current constellation. 
%
%     The measurement model can be described as follows
%       yt   =   H_0 z_t + H_1 z_{t-1} + ... + H_{L-1} z_{t-L+1} + n_t
%       where n_t ~ CN( 0, s_y^2 I). That is, the elements of n_t are
%       i.i.d. complex Gaussian random variables with variance s_y^2.
%
%       Nt denotes the number of transmitters, Nr the number of
%       receivers, L the maximum memory and T the length of the measurement
%       sequence.


Scenario    =   2;
switch Scenario
    case 1
        Nt  =   1;      Nr  =   1;
        L   =   4;      sy2 =   0.3;
        H   =   zeros(Nr,Nt,L); H(1,1,1)    = 0.3; H(1,1,2) = 1.3;
        H(1,1,3) = 0.2; H(1,1,4) = -1.2;
        C   =   [0,   1+i,   1-i,    -1-i,   -1+i];  Q   =   4;
        a   =   zeros(Nt,1);    b   =   zeros(Nt,1);
        a   =   0.2;    b   =   0.1;
        T   =   20;
    case 2
        Nt  =   3;      Nr  =   10;
        L   =   4;      sy2 =   0.5;
        H = sqrt(1/2)*randn(Nr,Nt,L)+1i*sqrt(1/2)*randn(Nr,Nt,L);
        %H   =   zeros(Nr,Nt,L); H(:,:,1)    = 0.6*eye(Nt); H(:,:,2) = 0.3*eye(Nt);
        %H(:,:,3) = 0.2*eye(Nt); H(:,:,4) = -0.2*eye(Nt);
        C   =   [0,   1+i,   1-i,    -1-i,   -1+i]/sqrt(2);  Q   =   4;
        a   =   0.4*(0.2*rand(Nt,1)+ones(Nt,1));    
        b   =   0.1*(0.5*rand(Nt,1)+ones(Nt,1));
        T   =   1000;
end






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%    Generating the true sequences 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% First we initialise the matrices where we intend to store the sequences:
X       =   zeros(Nt,T);  % Sequences from t=1 to t=T for Nt transmitters
Y       =   zeros(Nr,T);  % Nr receivers observe data at T time instances.


% Here we list the active transmitters
Act     =   zeros(Nt,1);


% Now we can generate the state and observation sequences
for     t   =   1   :   T
    
    %   We sample from different distributions for active and passive
    %   transmitters 
    X(:,t)  =   (Act.*binornd(ones(Nt,1),1-b)+...
           (1-Act).*binornd(ones(Nt,1),1-a)).*unidrnd(Q,Nt,1);
    %   It may have been better to use randsample but I don't see how to
    %   use that to generate vectors of random numbers. 
 
    Act     =   X(:,t)>0;
    
    %   To generate measurements we first compute the noise free
    %   measurements by summing over the length of the memory of the
    %   channel
    Ypred   =   zeros(Nr,1);
    for r = 0 : min(L-1,t-1)
         Ypred   =   Ypred   +   H(:,:,r+1)*C(X(:,t-r)+1).';
    end        
    Y(:,t)  =  Ypred + sqrt(sy2/2)*randn(Nr,2)*[1;i];
end

% As a reality check for the X samples we can check that
% sum(X==0)/T*(1-a+b)/b is close to 1 for each transmitter.







%profile on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%    Generating samples using the PGAS kernel 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N_PF    =   3000;             %   Number of particles in the PF
N_PG    =   500;              %   Number of particles in the PGAS kernel
M       =   10;             %   Number of MCMC samples
X_PG    =   zeros(Nt,M,T);    %   Stores MCMC samples of trajectories
xc      =   zeros(Nt,T);      %   The particle that we condition on

% Note that X_PG(:,m,t) is a vector that describes the transmitted symbols
% at time t, according to the m'th sample of the sequence of symbols. 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TO THE READER:
% In what follows we run a Sequential Monte Carlo algorithm for m=1, which
% corresponds to Alg 1 in the PGAS paper, and then we run the PGAS kernel
% for m = 2, 3, ..., M. This kernel corresponds to Alg. 2 in the PGAS paper.
%
% To clarify the different steps, we include e.g. [Line 5] to point out
% that the code corresponds to line 5 in Alg. 2 in the PGAS paper. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%

for m = 1 : M
    %%
    % To know how much longer we have to wait:
    if mod(m,100)==0
        M-m
    end
    
    if m == 1
        % The particle filter used to generate the first trajectory makes
        % use of more particles
        N   =   N_PF; 
    else
        % than the PGAS kernel
        N   =   N_PG;
    end
    % We now initialize three important variables:
    Xt      =   zeros(Nt,N,T);  %   Stores particles within PGAS
    a_ind   =   zeros(N,T);     %   Stores the ancestor indices
    W       =   zeros(N,T);     %   Stores the particle weights
    
    % In our first version we also store the last L elements 
    % in the particle trajectories. 
    X0_hist =   zeros(Nt,N,L); % Contains x_{t-L+1:t}^i, i=1,...,N
    X1_hist =   zeros(Nt,N,L); % Contains x_{t-L:t-1}^i, i=1,...,N
    
    % NOTE 1: Xt also contains particles but to find out how these are
    % connected across time we need to consider the ancestor indices a_ind.
    % NOTE 2: it is not necessary to store X0_hist and X1_hist since the same
    % information can be obtained from a_ind and Xt. We can probably speed
    % up the current implementation by not introducing the two variables. 
    % NOTE 3: within the PGAS algorithm, we use X0_hist to compute the
    % weights for particles 1 to N-1 and X1_hist to compute weights that we
    % use to sample the ancestors for particle N. 
    
    % In order to use built-in functions in Matlab it is convenient for us to 
    % produce matrices containing the am and bm coefficients
    A   =   repmat(a,1,N);
    B   =   repmat(b,1,N);
    
    
for t = 1 : T 
    % To story away x_{t-L:t-1}^i for i = 1, 2, ..., N:
    X1_hist     =   X0_hist;
    if t == 1
        % At time t = 1 we sample the states from the prior at time 1.
        % We know that all transmitters were passive at time 0     [Line 1]   
        Xt(:,:,t)   =   binornd(ones(Nt,N),1-A).*unidrnd(Q,Nt,N);
        
        % Note that the particles do not yet have any ancestors that 
        % we can sample. However, if m > 1 the N'th particle should be the
        % particle that we condition on                            [Line 2]
        if m ~= 1
            Xt(:,N,t) = xc(:,t);
        end
        
        % For convenience, we store the particle trajectories here:
        X0_hist      =   cat(3,zeros(Nt,N,L-1),Xt(:,:,t));
        
    else
        % At later times we first sample the ancestors, using multinomial 
        % resampling, to obtain a(:,t) = ind                       [Line 5]
        [~,ind] = histc(rand(N,1), [0 cumsum(W(:,t-1)')]);
        
        % Then we propagate the selected particles from time t-1 to t, 
        % to obtain Xt(:,:,t)                                      [Line 5]
        Act         =   Xt(:,ind,t-1)>0;
        Xt(:,:,t)   =   (Act.*binornd(ones(Nt,N),1-B)+...
           (1-Act).*binornd(ones(Nt,N),1-A)).*unidrnd(Q,Nt,N);
       % Note: if we wish to improve the update rate for complex scenarios
       % we should probably try to improve how we propagate particles from
       % time t-1 to t. Specifically, it seems highly inefficient to
       % sample using "unidrnd(Q,Nt,N)" when Q and Nt are large. 
        
        % Within the PGAS algorithm, the ancestor index for the N'th
        % particle is sampled differently:
        if m ~= 1
            % Set the N'th particle to the particle that we condition on
            Xt(:,N,t) =   xc(:,t);                             %   [Line 6]
            
            % The ancestor probabilities (the weights) are obtained by
            % multiplying three factors. The objective here is to calculate
            % these factors in order to obtain what we denote w_a. [Line 7]
            
            
            % To compute the factor that depends on y, we currently use a 
            % double for-loop (certainly, this can be improved upon):
            logWY    =   zeros(N,L-1);
            for tau     =   t : min(t+L-2,T)
                Ypred   =   zeros(Nr,N);
                for r = 0 : min(L-1,tau-1)
                    % Symbols at different times contribute to the same measurement
                    % and r defines the delay.
                    % For small delays r, the receivers observe "symbols" 
                    % transmitted after time t, thus described by the sequence 
                    % that we condition on. For larger delays, the symbols 
                    % are instead described by particles at time t-1.
                    % The smallest delay for which the particle sequence in the 
                    % interval 1 to t-1 has an impact is when tau-r=t-1.
                    if tau-r>t-1
                        Ypred   =  Ypred + H(:,:,r+1)*repmat(C(xc(:,tau-r)+1).',1,N);
                    else 
                        Ypred   =  Ypred + H(:,:,r+1)*C(X1_hist(:,:,end+tau-r-t+1)+1);
                    end
                end
                Ydiff                =   Ypred - repmat(Y(:,tau),1,N);
                logWY(:,tau-t+1)     =   sum(-abs(Ydiff).^2/sy2,1);
            end
            logWY   =   sum(logWY,2);
            WY      =   exp(logWY-max(logWY));
            
            % Another factor represents the transition probabilities between
            % particle i at time t-1 and the particle that we condition on at
            % time t. To compute these probabilities we need to distinguish
            % between four events for each symbol and particle i:
            % 0    ->  0,   passive remains passive, Prob = am
            % 0    -> (>0), passive becomes active,  Prob = (1-am)/Q
            % (>0) -> (>0), active remains active,   Prob = (1-bm)/Q
            % (>0) -> 0,    active becomes passive,  Prob = bm
            WZ_mat  =   zeros(Nt,N);    % Transition probability for each element
            X1      =   Xt(:,:,t-1);    % The particles at time t-1
            Act1    =   X1>0;
            Act0    =   repmat(xc(:,t)>0,1,N);
            % We can now go through the four cases and compute transition
            % probabilities for each symbol and particle:
            WZ_mat  =   (1-Act1).*(1-Act0).*A + (1-Act1).*Act0.*(1-A)/Q ...
                + Act1.*Act0.*(1-B)/Q +Act1.*(1-Act0).*B;
            logWZ   =   sum(log(WZ_mat),1); % Log-transition probabilities for each particle
            WZ      =   exp(logWZ-max(logWZ))';
            
            % Finally, we can compute the weights of interest
            w_a     =   W(:,t-1).*WY.*WZ;
            w_a     =   w_a/sum(w_a);
            % from which we generate the N'th ancestor             [Line 8]
            ind(N) = find(rand(1) < cumsum(w_a),1); 
        end
        % We have now computed all the ancestor indices
        a_ind(:,t)  =   ind;        % Stores results of     [Lines 5 and 8]
        
        % We can also store away the particles recent histories,  ~[Line 9] 
        X0_hist     =   cat(3,X1_hist(:,ind,2:end),Xt(:,:,t));
    end
    
    
    % Compute the importance weights for all the particles, W(:,t).
    % First we compute the expected measurements for different particles
    % by considering the particles and their histories     [Lines 3 and 10]
    Ypred   =   zeros(Nr,N);
    for r = 0 : min(L-1,t-1)
         Ypred   =   Ypred   +   H(:,:,r+1)*C(X0_hist(:,:,end-r)+1);
    end        
    Ydiff       =   Ypred - repmat(Y(:,t),1,N);
    logW        =   sum(-abs(Ydiff).^2/sy2,1)';
    W(:,t)      =   exp(logW-max(logW));
    W(:,t)      =   W(:,t)/sum(W(:,t));
end % This marks the end of the for-loop over time, t.
    
    
% We can now compute the generated trajectories from the ancestor indices
%                                                                  [Line 9]
ind = a_ind(:,T);
for(t = T-1:-1:1)
    Xt(:,:,t) = Xt(:,ind,t);
    ind = a_ind(ind,t);
end
% Two comments:
%   --  Previously we only stored the most recent parts of the particle
%   trajectories, which is more efficient as long L<<T. Here we construct
%   the entire trajectories from Xt and a_ind.
%   --  It may be slightly faster to first select one trajectory at random
%   (see below) and then only construct that tractory, but this loop is
%   anyway very fast so it probably does not make much difference. 

% Finally we select one trajectory at random                      [Line 12]
J = find(rand(1) < cumsum(W(:,T)),1);
% that we can store away 
X_PG(:,m,:)     =   Xt(:,J,:);
xc              =   reshape(Xt(:,J,:),Nt,T);

%%
end



% profile viewer
% Conclusion: the command binornd is slow. Any faster alternatives?




%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%    Illustration of the obtained results 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% We can check the error probabilities which should be very small when the
% measurement noise is small.

E_prob  =   zeros(Nt,T);
for t = 1 : T 
    Xt_samples  =   reshape(X_PG(:,:,t),Nt,M);
    E_prob(:,t)   =   mean(abs(Xt_samples - repmat(X(:,t),1,M))>0, 2);
end
E_prob



% We should also study the update rate as a function of time: 
Update_rate     =   zeros(1,T);
for t = 1:T
    Xvec              =   reshape(X_PG(:,:,t),Nt,M);
    Update_rate(t)  =   mean(sum(abs(diff(Xvec')'),1)>0);
end
figure(3)
clf
plot(Update_rate,'LineWidth',1.5)
title('Update rate','FontSize',14)
xlabel('Time','FontSize',14)












