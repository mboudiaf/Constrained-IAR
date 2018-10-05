% DG_IAR_treeSearch.m
%
% -------------------------------------------------------------------------
% DESCRIPTION: Searches the huge tree for integers vector candidates. Uses
% the constraint on the baseline to prune the tree efficiently. The
% algorithm is written recursively
% -------------------------------------------------------------------------
% INPUTS:   k
%           z
%           N_dd        Real valued double differenced ambiguity estimates
%           Q_dd        Covariance matrix associated to real valued
%           ambiguities estimates
%           Z           Structure that will contain vector candidates
%           ksi         Size of search space
%           A
% -------------------------------------------------------------------------
% OUTPUTS: None. Only Z structure is updated with all integers vectors
% candidates
% -------------------------------------------------------------------------
% AUTHOR: Malik Boudiaf
%         May 2018
% -------------------------------------------------------------------------

function DG_IAR_treeSearch(k,z,N_dd,Q_dd,Z,ksi,A,H,y,R,cur_baseline_est,cur_sigma,prefl,mul,fSat,cur_state,obs,idx_obs)
n=max(size(N_dd));
[~,~,D,~] = DG_IAR_IntBootstrap(N_dd(1:k),Q_dd(1:k,1:k));
v=ksi^2-sum((z(1:k-1)-N_dd(1:k-1)).^2./diag(D(1:k-1,1:k-1)));
if v>=0
    delta=sqrt(D(end,end))*sqrt(v);
    lbound= N_dd(k)-delta;
    ubound= N_dd(k)+delta;
    validIntegers=(ceil(lbound):floor(ubound));
    for validInteger=validIntegers%% Potential integers that meet requirements
        potentialInt=validInteger;
%         [tree,new_parent]=tree.addnode(parent,potentialInt);
        [baseline_ok,new_baseline_est,new_sigma]=DG_IAR_BaselineReq([z;potentialInt],H,A,Q_dd,cur_sigma,cur_baseline_est,k,prefl,mul,R,fSat,cur_state,obs,idx_obs,y);
        if baseline_ok
            z=[z;potentialInt];
            if k~=n
%                 tree=DG_IAR_treeSearch(new_parent,tree,k+1,z,N_dd,Q_dd,Z,ksi,A,H,y,R,new_baseline_est,new_sigma,prefl,mul,fSat,cur_state,obs,idx_obs);
                DG_IAR_treeSearch(k+1,z,N_dd,Q_dd,Z,ksi,A,H,y,R,new_baseline_est,new_sigma,prefl,mul,fSat,cur_state,obs,idx_obs);
                z=z(1:end-1);
            else
                Z.append(z);
                z=z(1:end-1);
            end
        end
    end
end
end