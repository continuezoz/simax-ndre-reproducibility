function prob=multistate_problem(n,r,gamma,beta)
%MULTISTATE_PROBLEM One-dimensional multistate transport model of Li, Chu,
% Juang and Lin (IMA J. Numer. Anal. 31, 2011), eq. (1.8):
%   R' = B^- - R F^- - F^+ R + R B^+ R,  R(0) = 0,
% with F^+- = (I-F)D^+-, B^+- = B D^+-, D^+- = diag(k^+-) > 0 and
% low-rank nonnegative F = F1*F2', B = B1*B2' satisfying (1.7):
% sum_i (f_ij + b_ij) = gamma < 1 for every incoming state j.
% States are energy groups u_j in [0,1] with cross-sections sigma_j =
% 10^(3u_j); k^+- = sigma*(1 -+ beta). Emission profiles are r Gaussian
% bumps (backward ones shifted to lower energies), membership weights are
% the same bumps normalized per incoming state. In the notation
% X' = -A X - X D + X S X + Q of this paper: A = F^+, D = F^-, S = B^+,
% Q = B^- = B1*(D^- B2)'. With ndre_scaled the scaled problem is returned.
if nargin<2 || isempty(r), r=10; end
if nargin<3 || isempty(gamma), gamma=0.9999; end
if nargin<4 || isempty(beta), beta=0.1; end
u=linspace(0,1,n)'; sigma=10.^(3*u); kp=sigma*(1-beta); km=sigma*(1+beta);
c=((1:r)-0.5)/r; w=1/r;
bump=@(shift) exp(-(u-(c+shift)).^2/(2*w^2));
F1=bump(0); F1=F1./sum(F1,1); B1=bump(-w/2); B1=B1./sum(B1,1);
member=bump(0); member=member./sum(member,2);
phi=0.6; F2=gamma*phi*member; B2=gamma*(1-phi)*member;
G=km.*B2; DF2=kp.*F2; EF2=km.*F2; DB2=kp.*B2;
prob=struct('n',n,'p',n,'F',B1,'G',G,'normQ',sqrt(sum(sum((B1'*B1).*(G'*G)))), ...
    'blockWidth',2*r, ...
    'Aapply',@(X)kp.*X-F1*(DF2'*X),'Dapply',@(X)km.*X-F1*(EF2'*X), ...
    'DTapply',@(X)km.*X-EF2*(F1'*X),'Sapply',@(X)B1*(DB2'*X), ...
    'Asolve',@(X)smw(X,kp,F1,DF2),'Dsolve',@(X)smw(X,km,F1,EF2), ...
    'DTsolve',@(X)smw(X,km,EF2,F1), ...
    'Ashift',@(X,s)smw(X,kp-s,F1,DF2),'DTshift',@(X,s)smw(X,km-s,EF2,F1), ...
    'shiftRange',[min(kp),max(km)]);
prob.rightBasisApply=prob.DTapply; prob.rightBasisSolve=prob.DTsolve;
prob.family="multistate"; prob.seed=NaN; prob.parameter=[gamma,beta,r];
prob.factors=struct('kp',kp,'km',km,'F1',F1,'F2',F2,'B1',B1,'B2',B2);
if ndre_scaled, prob=ndre_setup(prob); end
end

function X=smw(X,d,U,V)
% (diag(d) - U*V')\X by the Sherman--Morrison--Woodbury formula.
Z=X./d; DU=U./d; X=Z+DU*((eye(size(U,2))-V'*DU)\(V'*Z));
end
