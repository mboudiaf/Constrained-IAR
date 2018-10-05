% DG_IAR_BaselineReq.m
%
% -------------------------------------------------------------------------
% DESCRIPTION: Verifies whether the current candidate's associated baseline
% is not too far from filter's estimated baseline. Used to prune the search
% tree
% -------------------------------------------------------------------------
% INPUTS:   k
%           [baseline_est;z]           current candidate to be tested
%           N_dd        Real valued double differenced ambiguity estimates
%           Q_dd        Covariance matrix associated to real valued
%           ambiguities estimates
%           Z           Structure that will contain vector candidates
%           ksi         Size of search space
% -------------------------------------------------------------------------
% OUTPUTS:  bool                        Is the test passed or not
%           baseline_est                Current estimated baseline in the
%           tree
% -------------------------------------------------------------------------
% AUTHOR: Malik Boudiaf
%         May 2018
% -------------------------------------------------------------------------

function [bool,baseline_est,sigma]=DG_IAR_BaselineReq(z,H,A,Q_dd,sigma,baseline_est,k,prefl,mul,R,fSat,x,obs,idx_obs,y)

%% Augment the model with the new potential ambiguity and derive new matrices
Ck=[H(1:k,:) A(1:k,1:k)];
Rk=R(1:k,1:k);
X=[baseline_est;z];
sigma=blkdiag(sigma,Q_dd(k,k));
exp_meas=Ck*X;
true_meas=A(1:k,1:k)*y(1:k);

%% Get modeled measurements
% Earth Orientation Parameters
% MJD  = s3_epoch_caltomjd(fSat.date_meas);
% eopv = s3_eop_geteop(fSat.lsdata,fSat.eop,MJD);
% dUTC = eopv(8);
% dUT1 = eopv(4);
% xp   = eopv(2);
% yp   = eopv(3);
% dim_x   = fSat.dim_x;
% C_LIGHT = fSat.C_LIGHT;
% 
% R_eci2ecef = s3_refsys_rotation_ecitoecef(MJD,dUTC,dUT1,xp,yp,fSat.eopmdl);
% GPSt  = s3_epoch_caltogps(fSat.date_meas);
% GNSSt = GNSStime(GPSt(1),GPSt(2));
% prn_ref=obs(idx_obs(end,1),1).prn;
% sys_ref=obs(idx_obs(end,1),1).system;
% band_ref=obs(idx_obs(end,1),1).band;
% lam_ref= s3_getGNSSwavelength(band_ref);
% eph_ref=s3_GNSSgetephemeris(fSat.eph,sys_ref,prn_ref,GNSSt);
% idx_c  = 10    + dim_x;
% cbiasa = x(10)/C_LIGHT;
% cbiasb = x(idx_c)/C_LIGHT;
% 
% r_rcva=R_eci2ecef*x(1:3);
% r_rcvb=r_rcva+R_eci2ecef*baseline_est;
% cpAref = s3_GNSScarrierphase(band_ref,eph_ref,GNSSt,r_rcva,cbiasa,x(10+idx_obs(end,1)))*lam_ref;
% cpBref = s3_GNSScarrierphase(band_ref,eph_ref,GNSSt,r_rcvb,cbiasb,x(idx_c+idx_obs(end,2)))*lam_ref;
% sdcp_ref_exp=cpAref-cpBref;
% sdcp_ref_true=(obs(idx_obs(end,1),1).carrierphase-obs(idx_obs(end,2),2).carrierphase)*lam_ref;
% exp_meas=[];
% true_meas=[];
% for j=1:k
%     prn=obs(idx_obs(j,1),1).prn;
%     sys=obs(idx_obs(j,1),1).system;
%     eph=s3_GNSSgetephemeris(fSat.eph,sys,prn,GNSSt);
%     band=obs(idx_obs(j,1),1).band;
%     lam=s3_getGNSSwavelength(band);
%     
%     cpA = s3_GNSScarrierphase(band,eph,GNSSt,r_rcva,cbiasa,x(10+idx_obs(j,1)))*lam;
%     cpB = s3_GNSScarrierphase(band,eph,GNSSt,r_rcvb,cbiasb,x(idx_c+idx_obs(j,2)))*lam;
%     
%     sdcp=cpA-cpB;
%     ddcp_exp=sdcp-sdcp_ref_exp;
%     
%     cpA=obs(idx_obs(j,1),1).carrierphase*lam;
%     cpB=obs(idx_obs(j,2),2).carrierphase*lam;
%     
%     sdcp=cpA-cpB;
%     ddcp_true=sdcp-sdcp_ref_true;
%     
%     exp_meas=[exp_meas;ddcp_exp];
%     true_meas=[true_meas;ddcp_true];
% end

%% Updating our baseline estimate (EKF meas update)

K=Ck'/(Ck*sigma*Ck'+Rk);
X=X+sigma*K*(true_meas-exp_meas);
sigma=sigma-sigma*K*Ck*sigma;
baseline_est=X(1:3);
sigmal=sqrt(trace(sigma(1:3,1:3)));



%% Testing estimated baseline
l_est=norm(baseline_est);
% nu1=atan2d(baseline_est(3),(sqrt(baseline_est(2)^2+baseline_est(1))^2));
% nu2=atan2d(baseline_est(2),baseline_est(1));
% if abs(l-mul)<=mu*sl && abs(nu1-munu1)<=munu1*snu1 && abs(nu2-munu2)<=munu2*snu2
if abs(l_est-mul)<=prefl*sigmal
    bool=true;
else
    bool=false;
end
end