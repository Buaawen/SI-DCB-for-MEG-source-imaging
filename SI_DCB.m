function [s]=SI_DCB(B,L,G)
%==========================================================================
% Input: B :             MEG measurements
%        L :             Lead field matrix
%        G :             DCB


%==========================================================================
tic

[nSensor,nSource]=size(L);                      % number of sensors and sources
nSnap=size(B,2);                                % number of time points  


% Default Control Parameters
Cov_n                   = ones(nSensor,1);     % noise's convariance
epsilon                 = 1e-6;             % stop criterion
MAX_ITER                = 50;               % maximum iteration numbers
print                   = 1;                % whether print the iteration process
Y=B;



l0          = ones(nSensor,1);  % initial values of system noise's convariance
q          = ones(nSource,1)*nSensor/trace(L*L');  % initial values of system noise's convariance
cost_old   = 0;
s          = zeros(nSource,nSnap);                 % initial state estimates
evidence   = zeros(MAX_ITER,1);                    % evidence values
%##########################################################################
%                   Leadfield Matrix normalization
%##########################################################################

%%
%ADD BY WEN
suoyin=[];
for i=1:nSource
ind=find(G(:,i)>0);
suoyin{i}=ind;
end

gamma          = ones(nSource,1);%*nSensor/trace(L*L');  % initial values of system noise's convariance
gamma_g=G*gamma;
q=gamma_g;

%%===================================================================
%                      iteration
%==========================================================================
    fprintf('\nRunning SI-DCB for MEG source localization...\n');


for iter = 1:MAX_ITER
    tic;
    gamma_old=gamma;
    gammag_old=q;
    %-----------estimate s0 s1----------------------%

    sigma_y0=diag(l0)+L*(q.*L');
    invsigma_y0=inv(sigma_y0);
    DTinvSig=  L'*invsigma_y0;
    s=speye(nSource)\(q.*DTinvSig*(B));

    %--------------------------------update q0 q1---------------------------------%
    
    
T=nSnap;
norm_sfs=norms(s');


norm_sfs2=norm_sfs'.*norm_sfs';

invga=(1./(gammag_old)).*(1./(gammag_old));
norm_invga=norm_sfs2.*invga;

digdt=[];
for i=1:nSource
    Di=L(:,i);
    digdt(i,1)=Di'*invsigma_y0*Di;
end
    for i=1:nSource
        gammaA=sum(G(suoyin{i},i).*norm_invga(suoyin{i}));%
        gammaB=T*sum(G(suoyin{i},i).*digdt(suoyin{i}));%
        gamma(i)=gamma_old(i)*(abs(gammaA/gammaB));%

    end
    gamma_g2=G*gamma;%
    gamma_g=gamma_g2;
q=gamma_g;
%--------------------------update l0---------------------------------%
     d_y=Y-L*s;
    for i = 1:nSensor
        %----------Gradient descent----------%
        l0(i) = norm(d_y(i,:),'fro')/sqrt(nSnap*invsigma_y0(i,i));
    end

%--------------------------check stop condition---------------------------%

%% identical {a_i}
cost = -0.5*(trace(B*B'/sigma_y0)+ nSnap*log(det(sigma_y0)) + nSensor*nSnap*log(2*pi));% log p(B,a,q)

MSE = ((cost-cost_old)/cost);
cost_old = cost;
if 1
    disp(['SI-DCB Iteration: ',num2str(iter),'       MSE: ',num2str(MSE),' ])%,'  ])
end
if abs(MSE)<epsilon
    break;
end
evidence(iter) = cost;






toc;
 end

 evidence=evidence(1:iter-1);
 fprintf('\n SI-DCB Finished!\n');
 disp(' ');
 toc