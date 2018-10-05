% DG_IAR_constrainedIAR.m
%
% -------------------------------------------------------------------------
% DESCRIPTION: Resolves carrier-phase ambiguities to integer values using
%              a constrained tree search 
% -------------------------------------------------------------------------
% INPUTS:   x       - vector of unresolved ambiguities
%           Q       - subset of satellite covariance matrix containing
%                     abiguity entries
%           fSat    - Structure containing all information about each
%                   satellite within constellation (including absolute 
%                    position/vel as estimated by the filter
%           obs     -current measurements
%           P       -covariance matrix of the filter
%                   
%     
% -------------------------------------------------------------------------
% OUTPUTS:  N_dd_fix       - best candidate of fixed double differenced ambiguities
%           N_dd_fix_2     - second best candidate of fixed dd amb
%           
% -------------------------------------------------------------------------
% AUTHOR: Malik Boudiaf
%         April 10, 2018
% -------------------------------------------------------------------------

function [N_dd_fix,idx_dd_fix,N_dd_fix_2]=DG_IAR_contrainedIAR(x,fSat,N_dd,Q_dd,obs,y,idx_obs,P)

dim_x=length(x)/2;
nb_amb_unfixed=length(y);
idx_dd_fix=1:nb_amb_unfixed;
z=[];
R=fSat.sig_ddcp^2*eye(nb_amb_unfixed);
H=DG_IAR_computeH(fSat,x,obs,idx_obs,xtrue_c,xtrue_d);
A=DG_IAR_computeA(obs,idx_obs);
ksi=sqrt((round(N_dd)-N_dd)'/Q_dd*(round(N_dd)-N_dd)); %size of ellipsoidal search space
idx=dim_x+(1:3);
cur_baseline_est= x(idx)-x(1:3); %Baseline estimated by the filter

cur_sigma= diag(diag(P(1:3,1:3))+diag(P(idx:idx))); % Uncertainty on estimated baseline
mul=norm(cur_baseline_est); % Our initial guess on the baseline length
prefactors=linspace(0,5,12); %defining tolerance on the tree pruning
nb_candidates=[];
j=0;
for pref=prefactors
    Z=hmatrix;
    DG_IAR_treeSearch(1,z,N_dd,Q_dd,Z,ksi,A,H,y,R,cur_baseline_est,cur_sigma,pref,mul,fSat,x,obs,idx_obs);
    if isempty(Z.array)
       [N_dd_fix,~,~,~,~]= DG_IAR_mLambda(N_dd,Q_dd);
    else
        N_dd_fix=DG_IAR_bestCandidate(Z,N_dd,Q_dd);
        tic;
    end
    nb_candidates=[nb_candidates size(Z.array,2)];
    j=j+1;
end
    
if Z.nonEmpty()
    N_dd_fix_2=DG_IAR_bestCandidate(Z,N_dd,Q_dd);
else
    [N_dd_fix_2,~,~,~] = DG_IAR_IntBootstrap(N_dd,Q_dd);
end
end