%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   December 1, 2014
%   In this file we intend to evaluate hom PGAS performs on linear and
%   Gaussian state space models of different dimensionality. We will study
%   mixing properties in terms of update rates as well as by comparing the
%   obtained posterior moments with the exact moments given by an RTS
%   smoother.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%    Defining the model parameters 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%  We use the following notation:
% The state space model (SSM)
%       xt   =   F xt1  +   qt,  qt ~ N(0,Q),  Motion model
%       yt   =   H xt   +   rt,  rt ~ N(0,R),  Measurement model
%       x1   ~   N (x10, P10),                 Prior at time 1
%  We use nx and ny to denote the dimension of xt and yt,
%  respectively, and T to denote the length of the sequence.


Scenario    =   2;
switch Scenario
    case 1
        nx  =   1;      ny  =   1;
        F   =   1;      Q   =   1.5;  
        H   =   1;      R   =   2;
        x10 =   4;      P10 =   5;  T   =   30;
    case 2
        nx  =   10;          ny  =   nx;
        F   =   eye(nx);    Q   =   eye(nx);
        H   =   eye(nx);    R   =   eye(nx);
        x10 =   ones(nx,1); P10 =   5*eye(nx);
        T   =   30;
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%    Generating the true sequences 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% First we initialise the matrices where we intend to store the sequences:
X       =   zeros(nx,T); % sequences go from t=1 to t=T
Y       =   zeros(ny,T);   

% Now we can generate the state and observation sequences
for     t   =   1   :   T
    if  t   ==  1
        X(:,t)  =   mvnrnd(x10',P10)';
        Y(:,t)  =   mvnrnd(H*X(:,t),R)';
    else
        X(:,t)  =   mvnrnd(F*X(:,t-1),Q);
        Y(:,t)  =   mvnrnd(H*X(:,t),R)';
    end
end






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%    Computing posterior moments using RTS smoothing 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% To compute the exact posterior moments we first run a Kalman filter and
% then an RTS smoother

% Let us store the posterior moments in the following matrices:
Xtt     =   zeros(nx,T);    % Contains E[xt | y1,..,yt],
Ptt     =   zeros(nx,nx,T); % Cov[ xt | y1,..,yt], for t=1,...,T
Xt1     =   zeros(nx,T);    % Predicted mean
Pt1     =   zeros(nx,nx,T); % and covariance
XtT     =   zeros(nx,T);    % Smoothed posterior mean
PtT     =   zeros(nx,nx,T); % and covariance


% As a first step we run a Kalman filter
for     t   =   1   :   T
    %   The prediction step
    if  t   ==  1
        Xt1(:,t)        =   x10;
        Pt1(:,:,t)      =   P10;
    else
        Xt1(:,t)        =   F*Xtt(:,t-1);
        Pt1(:,:,t)      =   F*Ptt(:,:,t-1)*F' + Q;
    end
    
    %   The update step
    Stt         =   H * Pt1(:,:,t) * H' + R;    % Innovation covariance
    Kt          =   Pt1(:,:,t) * H' * Stt^(-1); % Kalman gain
    Xtt(:,t)    =   Xt1(:,t) + Kt * (Y(:,t) - H*Xt1(:,t)); % Posterior mean
    Ptt(:,:,t)  =   Pt1(:,:,t) - Kt*Stt*Kt';
end

% Given the Kalman filter output, we now run an RTS smoothing algorithm
for     t   =   T:-1:1
    if  t   ==  T
        XtT(:,t)    =   Xtt(:,T); % Filtering solution is identical to 
        PtT(:,:,t)  =   Ptt(:,:,T); %smoothing at time T.
    else
    Gt          =   Ptt(:,:,t) * F' * Pt1(:,:,t+1)^(-1); % RTS gain
    XtT(:,t)    =   Xtt(:,t) + Gt*(XtT(:,t+1)-Xt1(:,t+1)); 
    PtT(:,:,t)  =   Ptt(:,:,t)+Gt*(PtT(:,:,t+1)-Pt1(:,:,t+1))*Gt';
    end
end







%%

% profile on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%    Generating samples using the PGAS kernel 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N_PF    =   1000;           %   Number of particles in the PF
N_PG    =   5;              %   Number of particles in the PGAS kernel
M       =   3000;           %   Number of MCMC samples
X_PG    =   zeros(nx,M,T);  %   Stores MCMC samples of trajectories
xc      =   zeros(nx,T);    %   The particle that we condition on

% Note that X_PG(:,m,t) is the state vector for the m'th trajectory
% sample at time t. 


for m = 1 : M
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
    Xt      =   zeros(nx,N,T);  %   Stores particles within PGAS
    a       =   zeros(N,T);     %   Stores the ancestor indices
    W       =   zeros(N,T);     %   Store the particle weights
    
for t = 1 : T   
    if t == 1
        % At time t = 1 we sample the states from the prior at time 1. 
        Xt(:,:,t)   =   mvnrnd(x10',P10,N)';
        
        % Note that the particles do not yet have any ancestors that 
        % we can sample. However, if m > 1 the N'th particle should be the
        % particle that we condition on
        if m ~= 1
            Xt(:,N,t) = xc(:,t);
        end
    else
        % At later times we first sample the ancestors, using multinomial 
        % resampling, to obtain a(:,t) = ind
        [~,ind] = histc(rand(N,1), [0 cumsum(W(:,t-1)')]);
        
        % Then we propagate the particles from time t-1 to t, to obtain
        % Xt(:,:,t)
        Xpred       =   F*Xt(:,:,t-1); % We store this to save computations below.
        Xt(:,:,t)   =   Xpred(:,ind) + mvnrnd(zeros(1,nx),Q,N)';
        
        % Within the PGAS algorithm, the ancestor index for the N'th
        % particle is sampled differently:
        if m ~= 1
            % Set the N'th particle to the particle that we condition on
            Xt(:,N,t) =   xc(:,t);
            
            % The ancestor probabilities are
            w_a     =   W(:,t-1).*mvnpdf(repmat(xc(:,t)',N,1),Xpred',Q);
            w_a     =   w_a/sum(w_a);
            % from which we generate the N'th ancestor
            ind(N) = find(rand(1) < cumsum(w_a),1);
        end
        % We have now computed all the ancestor indices
        a(:,t)  =   ind;
    end
    % Compute the importance weights for all the particles, W(:,t).
    W(:,t)      =   mvnpdf(repmat(Y(:,t)',N,1),(H*Xt(:,:,t))',R);
    W(:,t)      =   W(:,t)/sum(W(:,t));
end

% We can now compute the generated trajectories from the ancestor indices
ind = a(:,T);
for(t = T-1:-1:1)
    Xt(:,:,t) = Xt(:,ind,t);
    ind = a(ind,t);
end

% Finally we select one trajectory at random
J = find(rand(1) < cumsum(W(:,T)),1);
% that we can store away 
X_PG(:,m,:)     =   Xt(:,J,:);
xc              =   reshape(Xt(:,J,:),nx,T);


end

% profile viewer
% Conclusion: mvnpdf and repmat slows down this implementation.
%             We can perform repmat(Y..) one time instead of M times.
%             We can probably also speed up mvnpdf. 





%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%    Illustration of the obtained results 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%   We should evaluate the solution in order to study its performance.
%   Suitable things to check:


%   1)  Posterior mean of PGAS versus RTS
X_PG_mean = zeros(nx,T);
for t = 1:T
    X_PG_mean(:,t)  =  mean(X_PG(:,100:end,t),2);
end

figure(1)
clf
plot(X_PG_mean(1,:),'LineWidth',1.5)
hold on
plot(XtT(1,:),'r-.','LineWidth',1.5)
legend('PGAS mean','RTS mean')
xlabel('Time','FontSize',14)
title('Posterior mean','FontSize',14)






%   2)  Posterior covariance of PGAS versus RTS
%       When nx>1 we simply consider the variance of the first element in
%       the state vector. 

P_PG    =   zeros(1,T);
for t = 1:T
    P_PG(t)  =  cov(X_PG(1,100:end,t));
end

figure(2)
clf
plot(P_PG,'LineWidth',1.5)
hold on
plot(reshape(PtT(1,1,:),T,1),'r-.','LineWidth',1.5)
axis([0 T 0 max(P_PG)*1.1])
legend('PGAS cov','RTS cov')
xlabel('Time','FontSize',14)
title('Posterior covariance','FontSize',14)







%   3)  The update rate as a function of time. 
Update_rate     =   zeros(1,T);
for t = 1:T
    Xvec            =   reshape(X_PG(1,:,t),M,1);
    Update_rate(t)  =   sum(abs(diff(Xvec))>0)/M;
end
figure(3)
clf
plot(Update_rate,'LineWidth',1.5)
title('Update rate','FontSize',14)
xlabel('Time','FontSize',14)





