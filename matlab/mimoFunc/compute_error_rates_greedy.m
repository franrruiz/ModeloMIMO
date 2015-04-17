function [ADER SER_ALL SER_ACT MMSE vec_ord rot ADER_indiv SER_ALL_indiv SER_ACT_indiv MMSE_indiv] = compute_error_rates_greedy(data,samples,hyper,param)
% Returns:
% -ADER: Activity detection error rate (i.e., prob that a trasmitter is
%        idle when it is said to be trasnmitting or vice-versa)
% -SER_ALL: Symbol error rate. It is computed as if 0 were another symbol of
%           the constellation
% -SER_ACT: Symbol error rate, conditioned on the actual trasmitter being
%           active
% -MMSE: Mean Square Error (w.r.t. the channel coefficients)
% -vec_ord: Vector containing the order of the estimated transmitters
%           needed to match the true ones
% -rot: Vector containing the coefficients to rotate the constellation in
%       order to match the true symbols
%
% The variables ending with "_indiv" correspond to the individual ADER, SER
% or MSE (for each transmitter)
% 

[Mest T] = size(samples.Z);
Nt = size(data.symbols,1);
despMax = max(param.L,10);
thr = min(abs(param.constellation))/10;
thrSER = 0.1;

if(Mest==Nt)
    vec_ord = zeros(1,Nt);
    rot = zeros(1,Nt)+1i*zeros(1,Nt);
    alreadyChosen = [];
    notChosen = 1:Nt;
    
    vec_ord_inf = zeros(1,Mest);
    rot_inf = zeros(1,Mest);
    
    % We match the inferred transmitters with the true ones
    for m=1:Mest
        flagChosen = 0;
        min_ser_m = inf;
        ll=-despMax;
        while(ll<=despMax && ~flagChosen)
            % Shifted replica of Z(m,:)
            Zm = [zeros(1,max(0,ll)) samples.Z(m,max(1-ll,1):min(T-ll,T)) zeros(1,max(0,-ll))];
            r = 0;
            while(r<=3 && ~flagChosen)
                % Rotate the constellation by a factor of pi*r/2
                Zm_rot = Zm*exp(1i*r*pi/2);
                % Compute the SER vs all transmitters
                ser_vs_all = sum(abs(repmat(Zm_rot,Nt,1)-data.symbols)>thr,2)/T;
                % If the SER of the m-th Tx has improved
                if(min(ser_vs_all(notChosen))<min_ser_m)
                    % Find those tx's for which the SER is below threshold
                    idx_all = find(ser_vs_all(notChosen)<thrSER);
                    if(isempty(idx_all))
                        [min_ser_m idxSelected] = min(ser_vs_all(notChosen));
                        vec_ord_inf(m) = notChosen(idxSelected);
                        rot_inf(m) = r;
                    else
                        [min_ser_m idxSelected] = min(ser_vs_all(notChosen));
                        vec_ord_inf(m) = notChosen(idxSelected);
                        rot_inf(m) = r;
                        alreadyChosen = [alreadyChosen notChosen(idxSelected)];
                        notChosen(notChosen==notChosen(idxSelected)) = [];
                        flagChosen = 1;
                    end
                end
                r = r+1;
            end
            ll = ll+1;
        end
    end
end






return;
if(Mest<Nt)
    samples.H = cat(2,samples.H,zeros(param.Nr,Nt-Mest,param.L));
	samples.Z = [samples.Z; zeros(Nt-Mest,size(samples.Z,2))];
    samples.seq = [samples.seq; zeros(Nt-Mest,size(samples.Z,2))];
	Mest = Nt;
end

if(param.L>size(data.channel,3))
    data.channel = cat(3,data.channel,zeros(param.Nr,Nt,param.L-size(data.channel,3)));
elseif(param.L<size(data.channel,3))
    warning('Cannot compute performance due to wrong choice of param.L');
    ADER = NaN;
    SER_ALL = NaN;
    SER_ACT = NaN;
    MMSE = NaN;
    vec_ord = NaN;
    rot = NaN;
    ADER_indiv = NaN;
    SER_ALL_indiv = NaN;
    SER_ACT_indiv = NaN;
    MMSE_indiv = NaN;
    return;
end

idxNZ = find(data.symbols~=0);
idxNZ_indiv = cell(1,Nt);
for m=1:Nt
    idxNZ_indiv{m} = find(data.symbols(m,:)~=0);
end

if(Mest==Nt)
    % Hay Mest! permutaciones
    permutac = perms(1:Mest);
    ADER = Inf;
    maxI = factorial(Nt);
    if(~flagExhaustiveOrder)
        maxI = 1;
    end
    for i=1:maxI
        ADER_aux = sum(sum((samples.seq(permutac(i,:),:)~=0)~=(data.seq~=0)))/(Nt*T);
        if(ADER_aux<ADER)
            vec_ord = permutac(i,:);
            ADER = ADER_aux;
            SER_ALL = Inf;
            maxM = 4^Nt-1;
            if(~flagExhaustiveRot)
                maxM = 0;
            end
            for m=0:maxM
                rot_aux = exp(1i*pi*de2bi(m,Nt,4)/2);
                Zaux = samples.Z(permutac(i,:),:).*repmat(rot_aux.',1,T);
                Haux = samples.H(:,permutac(i,:),:).*repmat(1./rot_aux,[param.Nr,1,param.L]);
                SER_ALL_aux = sum(abs(Zaux(:)-data.symbols(:))>thr)/(Nt*T);
                if(SER_ALL_aux<SER_ALL)
                    SER_ALL = SER_ALL_aux;
                    SER_ACT = sum(abs(Zaux(idxNZ)-data.symbols(idxNZ))>thr)/length(idxNZ);
                    rot = rot_aux;
                    MMSE = sum(abs(Haux(:)-data.channel(:)).^2)/numel(data.channel);
                    
                    ADER_indiv = sum((Zaux~=0)~=(data.seq~=0),2)/T;
                    SER_ACT_indiv = zeros(Nt,1);
                    SER_ALL_indiv = zeros(Nt,1);
                    MMSE_indiv = zeros(Nt,1);
                    for nt=1:Nt
                        SER_ACT_indiv(nt) = sum(abs(Zaux(nt,idxNZ_indiv{nt})-data.symbols(nt,idxNZ_indiv{nt}))>thr)/length(idxNZ_indiv{nt});
                        SER_ALL_indiv(nt) = sum(abs(Zaux(nt,:)-data.symbols(nt,:))>thr)/T;
                        MMSE_indiv(nt) = sum(sum(abs(Haux(:,nt,:)-data.channel(:,nt,:)).^2))/numel(data.channel(:,nt,:));
                    end
                end
            end
        end
    end
elseif(Mest>Nt)
    combinac = combnk(1:Mest,Nt);
    ADER = Inf;
    for j=1:size(combinac,1)
        permutac = perms(combinac(j,:));
        idxNotUsed = setdiff(1:Mest,combinac(j,:));
        for i=1:factorial(Nt)
            ADER_aux = (sum(sum((samples.seq(permutac(i,:),:)~=0)~=(data.seq~=0)))+sum(sum(samples.seq(idxNotUsed,:)~=0)))/(Nt*T);
            if(ADER_aux<ADER)
                vec_ord = permutac(i,:);
                ADER = ADER_aux;
                SER_ALL = Inf;
                maxM = 4^Nt-1;
                if(~flagExhaustiveRot)
                    maxM = 0;
                end
                for m=0:maxM
                    rot_aux = exp(1i*pi*de2bi(m,Nt,4)/2);
                    Zaux = samples.Z(permutac(i,:),:).*repmat(rot_aux.',1,T);
                    Haux = samples.H(:,permutac(i,:),:).*repmat(1./rot_aux,[param.Nr,1,param.L]);
                    SER_ALL_aux = (sum(abs(Zaux(:)-data.symbols(:))>thr)+sum(sum(samples.seq(idxNotUsed,:)~=0)))/(Nt*T);
                    if(SER_ALL_aux<SER_ALL)
                        SER_ALL = SER_ALL_aux;
                        SER_ACT = sum(abs(Zaux(idxNZ)-data.symbols(idxNZ))>thr)/length(idxNZ);
                        rot = rot_aux;
                        MMSE = sum(abs(Haux(:)-data.channel(:)).^2)/numel(data.channel);
                        
                        ADER_indiv = sum((Zaux~=0)~=(data.seq~=0),2)/T;
                        SER_ACT_indiv = zeros(Nt,1);
                        SER_ALL_indiv = zeros(Nt,1);
                        MMSE_indiv = zeros(Nt,1);
                        for nt=1:Nt
                            SER_ACT_indiv(nt) = sum(abs(Zaux(nt,idxNZ_indiv{nt})-data.symbols(nt,idxNZ_indiv{nt}))>thr)/length(idxNZ_indiv{nt});
                            SER_ALL_indiv(nt) = sum(abs(Zaux(nt,:)-data.symbols(nt,:))>thr)/T;
                            MMSE_indiv(nt) = sum(sum(abs(Haux(:,nt,:)-data.channel(:,nt,:)).^2))/numel(data.channel(:,nt,:));
                        end
                    end
                end
            end
        end
    end
else
    error(['Mest y Nt no son comparables: Mest=' num2str(Mest) ' y Nt=' num2str(Nt)]);
end

