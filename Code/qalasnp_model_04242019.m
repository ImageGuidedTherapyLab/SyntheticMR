function Mmeas = qalasnp_model_04242019(x, parnames, parvals)

% y = line_model(x, parnames, parvals)
%
% This function will return a line given by the equation y = mx + b, where
% m is the line's gradient and b is its y-intersept. The
% input parameters are:
%   x - the x values at which y will be calculated
%   parnames - a cell array containing the parameters names. These can be
%       in any order, but must include the following parameters:
%           {'m', 'b'}
%   parvals - a cell array containing the values of the parameters given in
%       parnames. These must be in the same order as in parnames.
%
%--------------------------------------------------------------------------
%           This is the format required by nested_sampler.m.
%--------------------------------------------------------------------------

% check that parnames and parvals have the same length
lpn = length(parnames);
lpv = length(parvals);
if lpn ~= lpv
    error('Error: parnames and parvals are not the same length!');
end

nparams = lpn;

% extract parameter values
% for ii=1:((lpn-5)/3)
%     M0(ii)=parvals{3*ii-2};
%     T1(ii)=parvals{3*ii-1};
%     T2(ii)=parvals{3*ii};
% %     eval(sprintf('M0_%i=parvals{ii};',ii));
% %     eval(sprintf('T1_%i=parvals{ii};',ii));
% %     eval(sprintf('T2_%i=parvals{ii};',ii));
% end

for ii=1:nparams
    switch parnames{ii}
        case 'M0_GM'
            M0_GM = parvals{ii};
        case 'T1_GM'
            T1_GM = parvals{ii};
        case 'T2_GM'
            T2_GM = parvals{ii};
        case 'M0_WM'
            M0_WM = parvals{ii};
        case 'T1_WM'
            T1_WM = parvals{ii};
        case 'T2_WM'
            T2_WM = parvals{ii};
        case 'M0_CSF'
            M0_CSF = parvals{ii};
        case 'T1_CSF'
            T1_CSF = parvals{ii};
        case 'T2_CSF'
            T2_CSF = parvals{ii};
        case 'TR'
            TR = parvals{ii};
        case 'TE_T2prep'
            TE_T2prep = parvals{ii};
        case 'flipAngle'
            flipAngle = parvals{ii};
        case 'nacq'
            nacq = parvals{ii};
        case 'dt'
            dt = parvals{ii};
        case 'materialID'
            materialID = parvals{ii};
    end
end

    switch materialID
        case 1
            M0=M0_GM; T1=T1_GM; T2=T2_GM;
        case 2
            M0=M0_WM; T1=T1_WM; T2=T2_WM;
        case 3
            M0=M0_CSF; T1=T1_CSF; T2=T2_CSF;
    end
    % szm=size(M0);
    star=(1-exp(-TR./T1(:)'))./(1-cosd(flipAngle).*exp(-TR./T1(:)'));
    % T2 sensitization
    M(1,:)=M0(:)';        % assume initialization?
    M(2,:)=M(1,:).*exp(-TE_T2prep./T2(:)');
    M(3,:)=M0(:)'.*star-(M0(:)'.*star-M(2,:)).*exp(-dt(3)./(T1(:)'.*star));
    M(4,:)=M0(:)'-(M0(:)'-M(3,:)).*exp(-dt(4)./T1(:)');
    % T1 sensitization
    M(5,:)=-M(4,:);             % Assume perfect inversion? 100 ms inversion time
    M(6,:)=M0(:)'-(M0(:)'-M(5,:)).*exp(-dt(6)./T1(:)');
    % Post T1 sens acquisitions
    for iii=1:nacq-1
        M(5+2*iii,:)=M0(:)'.*star-(M0(:)'.*star-M(4+2*iii,:)).*exp(-dt(5+2*iii)./(T1(:)'.*star));
        M(6+2*iii,:)=M0(:)'-(M0(:)'-M(5+2*iii,:)).*exp(-dt(6+2*iii)./T1(:)');
    end
    Mmeas=sind(flipAngle).*[M(2,:);M(6:2:end-1,:)];
    
    for jjj=1:100
        %     [M,Mmeas]=qalas1time(M(end),M0,T1,T2,TR,TE_T2prep,flipAngle,nacq,dt);
        star=(1-exp(-TR./T1(:)'))./(1-cosd(flipAngle).*exp(-TR./T1(:)'));
        % T2 sensitization
        M(1,:)=M(end,:);        % assume initialization?
        M(2,:)=M(1,:).*exp(-TE_T2prep./T2(:)');
        M(3,:)=M0(:)'.*star-(M0(:)'.*star-M(2,:)).*exp(-dt(3)./(T1(:)'.*star));
        M(4,:)=M0(:)'-(M0(:)'-M(3,:)).*exp(-dt(4)./T1(:)');
        % T1 sensitization
        M(5,:)=-M(4,:);             % Assume perfect inversion? 100 ms inversion time
        M(6,:)=M0(:)'-(M0(:)'-M(5,:)).*exp(-dt(6)./T1(:)');
        % Post T1 sens acquisitions
        for iii=1:nacq-1
            M(5+2*iii,:)=M0(:)'.*star-(M0(:)'.*star-M(4+2*iii,:)).*exp(-dt(5+2*iii)./(T1(:)'.*star));
            M(6+2*iii,:)=M0(:)'-(M0(:)'-M(5+2*iii,:)).*exp(-dt(6+2*iii)./T1(:)');
        end
        Mmeas=sind(flipAngle).*[M(2,:);M(6:2:end-1,:)];
        if norm(M(1,:)-M(end,:))<=0.0001; break; end;
    end
    % Mmeas=reshape(Mmeas,[5,szm]);
%     switch tisind
%         case 1
%             Mmeas_GM=Mmeas(:);
%         case 2
%             Mmeas_WM=Mmeas(:);
%         case 3
%             Mmeas_CSF=Mmeas(:);
%     end    


% Mmeas=[Mmeas_GM,Mmeas_WM,Mmeas_CSF];

return
