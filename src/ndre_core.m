function api=ndre_core
%NDRE_CORE Shared kernels: the single-state transport problem, Galerkin
% projection, the stationary solution of the reduced equation, the centered
% flow, scaled Radon composition (mobius), error curves, and the EBA--BDF1
% baseline of Angelova, Hached and Jbilou (2020): extended block Arnoldi
% projection with backward Euler integration of the projected equation.
api=struct('problem',@make_problem,'options',@default_options, ...
    'adaptive',@adaptive_solve,'basis',@build_basis, ...
    'project',@project_problem,'prepare',@centered_prepare, ...
    'evaluate',@centered_evaluate,'errors',@trajectory_errors, ...
    'mobius',@mobius_trajectory);
end
function opts=default_options(prob)
r=size(prob.F,2);
opts=struct('tf',1,'h',0.01,'tol',1e-8/sqrt(2),'gate',1e-8, ...
    'newtonTol',1e-11,'newtonMax',50,'maxDim',320,'maxEBASteps',160/r, ...
    'checkInterval',5,'verbose',false,'saveTrajectory',true, ...
    'outputTimes',0:0.01:1,'checkDimension',10*r);
if prob.family~="transport", opts.maxDim=160; opts.maxEBASteps=40; end
end

function prob=make_problem(family,n,p,seed,parameter)
% Single-state transport model of section 2; parameter = [c alpha].
% p and seed are kept for the interface; the model has p = n.
if nargin<3 || isempty(p), p=n; end
if nargin<4 || isempty(seed), seed=20261007+n+p; end
if nargin<5, parameter=3; end
family=string(family);
assert(family=="transport",'ndre:family','Only the transport family is implemented.');
c=0.9999; alpha=1e-8;
if numel(parameter)==2, c=parameter(1); alpha=parameter(2); end
[omega,w]=load_gauss_legendre_matlab(n);
e=ones(n,1); q=w./(2*omega);
a=1./(c*omega*(1+alpha)); d=1./(c*omega*(1-alpha));
prob=struct('n',n,'p',n,'F',e,'G',e,'normQ',n,'blockWidth',2, ...
    'Aapply',@(X)a.*X-e*(q'*X),'Dapply',@(X)d.*X-q*(e'*X), ...
    'DTapply',@(X)d.*X-e*(q'*X),'Sapply',@(X)q*(q'*X), ...
    'Asolve',@(X)rankone_inverse(X,a,e,q),'Dsolve',@(X)rankone_inverse(X,d,q,e), ...
    'DTsolve',@(X)rankone_inverse(X,d,e,q), ...
    'Ashift',@(X,s)rankone_inverse(X,a-s,e,q), ...
    'DTshift',@(X,s)rankone_inverse(X,d-s,e,q),'shiftRange',[min(a),max(a)]);
prob.rightBasisApply=prob.DTapply; prob.rightBasisSolve=prob.DTsolve;
prob.family=family; prob.seed=seed; prob.parameter=parameter;
if family=="transport" && ndre_scaled, prob=ndre_setup(prob); end
end

function X=rankone_inverse(B,d,u,v)
Z=B./d; q=u./d; X=Z+q*((v'*Z)/(1-v'*q));
end

function [V,W]=build_basis(prob,basisName,dimension)
% Extended block Arnoldi bases of the given dimension (references).
assert(string(basisName)=="EBA",'ndre:basis','Only the EBA basis is implemented.');
[V,W]=extended_pair(prob,dimension);
end

function result=adaptive_solve(prob,method,opts)
% Published EBA--BDF1 solver with its own endpoint stopping rule.
assert(string(method)=="EBA-BDF1",'ndre:method','Only EBA-BDF1 is implemented.');
result=solve_original_eba_bdf1(prob,opts); result.accepted=result.converged;
result.indicator=result.boundaryFrobeniusResidual;
result.branchMargin=NaN; result.middleRcond=NaN;
end

function curve=trajectory_errors(a,b)
assert(numel(a.time)==numel(b.time) && max(abs(a.time-b.time))<1e-12,'Different time grids.');
[~,RV]=qr([a.V,b.V],0); [~,RW]=qr([a.W,b.W],0);
count=numel(b.time); absoluteError=zeros(count,1); relativeError=NaN(count,1);
referenceNorm=zeros(count,1);
for j=1:count
    absoluteError(j)=norm(RV*blkdiag(a.Yhistory(:,:,j),-b.Yhistory(:,:,j))*RW','fro');
    referenceNorm(j)=norm(b.Yhistory(:,:,j),'fro');
    if referenceNorm(j)>0, relativeError(j)=absoluteError(j)/referenceNorm(j); end
end
curve=table(b.time(:),referenceNorm,absoluteError,relativeError, ...
    'VariableNames',{'time','referenceNorm','absoluteError','relativeError'});
end

function sol=mobius_trajectory(red,times,substeps,extraScale)
% Independent scaled composition of the Riccati linear-fractional flow.
% F(Y)=Qh+Lh*Y/(I-Sh*Y)*Rh. No stationary Riccati solve is used.
% Store Lh-I and Rh-I during doubling to avoid losing slow-mode increments.
if nargin<3, substeps=1; end
if nargin<4, extraScale=0; end
k=size(red.A,1); ell=size(red.D,1);
H=[red.D,-red.S;red.Q,-red.A]; h=(times(2)-times(1))/substeps;
scalings=max(0,ceil(log2(max(h*norm(H,1),realmin)/0.25)))+extraScale;
B=(h/2^scalings)*H; term=B; EminusI=B;
for q=2:40
    term=term*B/q; EminusI=EminusI+term;
    if norm(term,1)<=eps*norm(EminusI,1), break; end
end
U=eye(ell)+EminusI(1:ell,1:ell);
RminusI=-(U\EminusI(1:ell,1:ell));
Sh=-(U\EminusI(1:ell,ell+1:end));
Qh=EminusI(ell+1:end,1:ell)/U;
LminusI=EminusI(ell+1:end,ell+1:end)-Qh*EminusI(1:ell,ell+1:end);
minRcond=1;
for q=1:scalings
    L=eye(k)+LminusI; R=eye(ell)+RminusI;
    U=eye(k)-Qh*Sh; V=eye(ell)-Sh*Qh;
    minRcond=min([minRcond,rcond(U),rcond(V)]);
    assert(minRcond>1e-14,'ndre:mobiusDoubling','Doubling denominator is ill-conditioned.');
    nextQ=Qh+L*(U\Qh)*R;
    nextS=Sh+R*(V\Sh)*L;
    nextL=2*LminusI+LminusI*LminusI+L*(U\(Qh*Sh))*L;
    nextR=2*RminusI+RminusI*RminusI+R*(V\(Sh*Qh))*R;
    Qh=nextQ; Sh=nextS; LminusI=nextL; RminusI=nextR;
end
L=eye(k)+LminusI; R=eye(ell)+RminusI;
Y=zeros(k,ell); if isfield(red,'Y0'), Y=red.Y0; end
Yhistory=zeros(k,ell,numel(times)); Yhistory(:,:,1)=Y;
for j=2:numel(times)
    for s=1:substeps
        M=eye(ell)-Sh*Y; rc=rcond(M); minRcond=min(minRcond,rc);
        assert(rc>1e-14,'ndre:mobiusStep','Flow denominator is ill-conditioned.');
        Y=Qh+(L*Y/M)*R;
    end
    Yhistory(:,:,j)=Y;
end
sol=struct('Y',Y,'Yhistory',Yhistory,'time',times,'middleRcond',minRcond, ...
    'scalings',scalings,'substeps',substeps);
end

function result = solve_original_eba_bdf1(prob,opts)
% Build one EBA block at a time and stop with the residual of Theorem 3 of
% Angelova et al.
totalTimer = tic;
rankF = size(prob.F,2);
blockWidth = 2*rankF;

basisTimer = tic;
[V, lastV, nextV, boundaryA] = initialize_eba_side(prob.Aapply,prob.Asolve,prob.F,blockWidth);
if isfield(prob,'rightBasisApply')
    rightApply = prob.rightBasisApply;
    rightSolve = prob.rightBasisSolve;
else
    rightApply = prob.Dapply;
    rightSolve = prob.Dsolve;
end
[W, lastW, nextW, boundaryD] = initialize_eba_side(rightApply,rightSolve,prob.G,blockWidth);
basisSeconds = toc(basisTimer);

status = "MAX_EBA_STEPS";
converged = false;
boundaryResidual = NaN;
normalizedResidual = NaN;
liftedResidual = NaN;
boundaryFrobeniusResidual = NaN;
solverSeconds = 0;
maxNewtonResidual = 0;
maxNewtonBackwardError = 0;
ebaSteps = 0;
dimension = size(V,2);
curve = table([],[],[],[],[],[],[], ...
    'VariableNames',{'ebaSteps','dimension','leftResidual','rightResidual', ...
    'maxResidual','boundaryFrobeniusResidual','directFrobeniusResidual'});

for m = 1:opts.maxEBASteps
    % Generate the relation block T_{m+1,m} used by the stopping residual.
    basisTimer = tic;
    if m > 1
        [nextV, boundaryA] = expand_eba_side(prob.Aapply,prob.Asolve, ...
            V,lastV,prob.blockWidth);
        [nextW, boundaryD] = expand_eba_side(rightApply,rightSolve, ...
            W,lastW,prob.blockWidth);
        basisSeconds = basisSeconds + toc(basisTimer);
    end

    % As in the published method, the projected equation is solved every
    % opts.checkInterval EBA iterations.
    if mod(m,opts.checkInterval) ~= 0 && m < opts.maxEBASteps && ...
            size(nextV,2) == blockWidth && size(nextW,2) == blockWidth
        ebaSteps = m;
        dimension = size(V,2);
        V = [V,nextV]; %#ok<AGROW>
        W = [W,nextW]; %#ok<AGROW>
        lastV = nextV;
        lastW = nextW;
        continue
    end

    red = project_original(prob,V,W);
    solverTimer = tic;
    sol = integrate_projected_bdf1(red,opts);
    solverSeconds = solverSeconds + toc(solverTimer);
    maxNewtonResidual = max(maxNewtonResidual,sol.maxNewtonResidual);
    maxNewtonBackwardError = max(maxNewtonBackwardError,sol.maxNewtonBackwardError);
    Y = sol.Y;

    % Residual of the published method, formed from the boundary blocks.
    last = size(Y,1)-blockWidth+1:size(Y,1);
    leftBoundary = boundaryA*Y(last,:);
    rightBoundary = Y(:,last)*boundaryD';
    boundaryResidual = max(norm(leftBoundary,'fro'),norm(rightBoundary,'fro'));
    normalizedResidual = boundaryResidual / max(prob.normQ,eps);
    boundaryFrobeniusResidual = hypot(norm(leftBoundary,'fro'), ...
        norm(rightBoundary,'fro')) / max(prob.normQ,eps);
    % Full-operator residual, formed only on request.
    if isfield(opts,'diagnostics') && opts.diagnostics
        RA = prob.Aapply(V)-V*red.A;
        RD = prob.DTapply(W)-W*red.D';
        liftedResidual = hypot(norm(RA*Y,'fro'),norm(Y*RD','fro')) / max(prob.normQ,eps);
    else
        liftedResidual = NaN;
    end
    curve(end+1,:) = {m,size(V,2),norm(leftBoundary,'fro')/prob.normQ, ...
        norm(rightBoundary,'fro')/prob.normQ,normalizedResidual, ...
        boundaryFrobeniusResidual,liftedResidual}; %#ok<AGROW>
    ebaSteps = m;
    dimension = size(V,2);
    if normalizedResidual < opts.tol
        converged = true;
        status = "OK";
        break
    end

    if size(nextV,2) == 0
        status = "BREAKDOWN";
        break
    elseif size(nextV,2) ~= prob.blockWidth
        status = "PARTIAL_DEFLATION";
        break
    end

    if m == opts.maxEBASteps
        break
    end
    % The current next block becomes the last block for the next EBA step.
    V = [V,nextV]; %#ok<AGROW>
    W = [W,nextW]; %#ok<AGROW>
    lastV = nextV;
    lastW = nextW;
end

result = struct('ebaSteps',ebaSteps,'dimension',dimension, ...
    'boundaryResidual',boundaryResidual,'normalizedResidual',normalizedResidual, ...
    'boundaryFrobeniusResidual',boundaryFrobeniusResidual, ...
    'liftedContinuousResidual',liftedResidual,'basisSeconds',basisSeconds, ...
    'solverSeconds',solverSeconds,'totalSeconds',toc(totalTimer), ...
    'maxNewtonResidual',maxNewtonResidual, ...
    'maxNewtonBackwardError',maxNewtonBackwardError, ...
    'converged',converged,'status',status,'curve',curve);
result.V=V; result.W=W; result.Y=Y; % Retain solution for untimed error checks.
if isfield(sol,'Yhistory')
    result.Yhistory=sol.Yhistory; result.time=sol.time;
end
end


function [V,lastBlock,nextBlock,boundary] = initialize_eba_side(applyM,solveM,F,blockWidth)
candidate = [F,solveM(F)];
V = append_orthogonal_block(candidate,zeros(size(F,1),0));
if size(V,2) ~= blockWidth
    error('ndre:rankBreakdown','Initial EBA block is rank deficient.');
end
lastBlock = V;
[nextBlock,boundary] = next_eba_block(applyM,solveM,V,blockWidth,V);
end


function [nextBlock,boundary] = expand_eba_side(applyM,solveM,V,lastBlock,blockWidth)
[nextBlock,boundary] = next_eba_block(applyM,solveM,lastBlock,blockWidth,V);
end


function [nextBlock,boundary] = next_eba_block(applyM,solveM,lastBlock,blockWidth,V)
if nargin < 5, V = []; end
s = blockWidth/2;
candidate = [applyM(lastBlock(:,1:s)), solveM(lastBlock(:,s+1:end))];
nextBlock = append_orthogonal_block(candidate,V);
if isempty(nextBlock)
    boundary = zeros(0,size(lastBlock,2));
else
    boundary = nextBlock' * applyM(lastBlock);
end
end


function Q = append_orthogonal_block(candidate,V)
% Two-pass MGS followed by an unpivoted thin QR, as in a block EBA.
Z = candidate;
if ~isempty(V)
    Z = Z - V*(V'*Z);
    Z = Z - V*(V'*Z);
end
[Q,R] = qr(Z,0);
tol = max(size(candidate))*eps*max(norm(candidate,2),1);
if any(abs(diag(R)) <= tol)
    Q = Q(:,abs(diag(R)) > tol);
end
end


function red = project_original(prob,V,W)
red = struct('A',V'*prob.Aapply(V),'D',W'*prob.Dapply(W), ...
    'S',W'*prob.Sapply(V),'F',V'*prob.F,'G',W'*prob.G, ...
    'normQ',prob.normQ);
end


function sol = integrate_projected_bdf1(red,opts)
k = size(red.A,1); ell = size(red.D,1);
steps = round(opts.tf/opts.h);
Y = zeros(k,ell);
if isfield(red,'Y0'), Y=red.Y0; end
if isfield(red,'Q'), Q=red.Q; else, Q=red.F*red.G'; end
storeHistory=isfield(opts,'saveTrajectory') && opts.saveTrajectory;
if storeHistory
    if isfield(opts,'outputTimes'), times=opts.outputTimes; else, times=(0:steps)*opts.h; end
    sampleSteps=round(times/opts.h);
    assert(max(abs(sampleSteps*opts.h-times))<1e-12,'Output nodes must be BDF steps.');
    Yhistory=zeros(k,ell,numel(times)); Yhistory(:,:,1)=Y; nextSample=2;
end
maxNewtonResidual = 0;
maxNewtonBackwardError = 0;
for j = 1:steps
    base = Y;
    h = opts.h;
    Yn = Y;
    for it = 1:opts.newtonMax
        f = -red.A*Yn - Yn*red.D + Yn*red.S*Yn + Q;
        G = Yn - base - h*f;
        scale = norm(Yn,'fro')+norm(base,'fro')+h*( ...
            norm(red.A*Yn,'fro')+norm(Yn*red.D,'fro')+ ...
            norm(Yn*red.S*Yn,'fro')+norm(Q,'fro'));
        if norm(G,'fro') <= opts.newtonTol*max(scale,eps)
            break
        end
        % Newton step: split the identity evenly between the two
        % Sylvester coefficients.
        J = 0.5*eye(k) + h*(red.A-Yn*red.S);
        K = 0.5*eye(ell) + h*(red.D-red.S*Yn);
        d = sylvester(J,K,-G);
        Yn = Yn + d;
    end
    f = -red.A*Yn - Yn*red.D + Yn*red.S*Yn + Q;
    relativeNewtonResidual = norm(Yn-base-h*f,'fro') / ...
        max(1,norm(Yn,'fro'));
    scale = norm(Yn,'fro')+norm(base,'fro')+h*( ...
        norm(red.A*Yn,'fro')+norm(Yn*red.D,'fro')+ ...
        norm(Yn*red.S*Yn,'fro')+norm(Q,'fro'));
    backwardError = norm(Yn-base-h*f,'fro') / max(scale,eps);
    maxNewtonResidual = max(maxNewtonResidual,relativeNewtonResidual);
    maxNewtonBackwardError = max(maxNewtonBackwardError,backwardError);
    % Newton aims at newtonTol; a step whose iteration stagnates below
    % 100*newtonTol is accepted, and the largest backward error is returned.
    if backwardError > 100*opts.newtonTol
        error('ndre:newtonNonconvergence', ...
            'Projected BDF Newton failed at step %d, dimension %d, backward error %.3e after %d iterations.', ...
            j,k,backwardError,it);
    end
    Y = Yn;
    if storeHistory && nextSample<=numel(times) && j==sampleSteps(nextSample)
        Yhistory(:,:,nextSample)=Y; nextSample=nextSample+1;
    end
end
sol = struct('Y',Y,'maxNewtonResidual',maxNewtonResidual, ...
    'maxNewtonBackwardError',maxNewtonBackwardError);
if storeHistory
    assert(nextSample==numel(times)+1,'Incomplete BDF trajectory.');
    sol.Yhistory=Yhistory; sol.time=times;
end
end


function red = project_problem(prob,V,W)
red = struct('V',V,'W',W,'A',V'*prob.Aapply(V), ...
    'D',W'*prob.Dapply(W),'S',W'*prob.Sapply(V), ...
    'Q',(V'*prob.F)*(W'*prob.G)', ...
    'normQ',prob.normQ,'prob',prob);
end


function prepared=centered_prepare(red,wantKernel)
% wantKernel=false skips the Sylvester kernel L, which only the direct
% centered evaluator needs; increment evaluation uses Ystar, Ac and Dc.
if nargin<2, wantKernel=true; end
ell = size(red.D,1);
H = [red.D,-red.S;red.Q,-red.A];
[U,T] = schur(H,'real');
select = real(eig(T)) > 0;
if sum(select) ~= ell
    error('ndre:branchCount','No stable graph branch of required size.');
end
[U,~] = ordschur(U,T,select);
Z = U(:,1:ell);
if rcond(Z(1:ell,:)) < 1e-14
    error('ndre:graph','The stationary graph is ill-conditioned.');
end
Ystar = Z(ell+1:end,:)/Z(1:ell,:);
qScale = max(1,norm(red.Q,'fro'));
for it = 1:30
    N = -red.A*Ystar-Ystar*red.D+Ystar*red.S*Ystar+red.Q;
    if norm(N,'fro') <= 1e-12*qScale, break; end
    Ac = red.A-Ystar*red.S; Dc = red.D-red.S*Ystar;
    Ystar = Ystar+sylvester(Ac,Dc,N);
end
N = -red.A*Ystar-Ystar*red.D+Ystar*red.S*Ystar+red.Q;
scale = norm(red.A*Ystar,'fro')+norm(Ystar*red.D,'fro')+ ...
    norm(Ystar*red.S*Ystar,'fro')+norm(red.Q,'fro');
nareBackwardError = norm(N,'fro')/max(scale,eps);
nareScaledResidual = norm(N,'fro')/qScale;
if nareScaledResidual >= 1e-8
    error('ndre:nare','Stationary NARE scaled residual %.3e.',nareScaledResidual);
end
Ac = red.A-Ystar*red.S; Dc = red.D-red.S*Ystar;
minClosedLoopRealPart = min(real([eig(Ac);eig(Dc)]));
branchMarginThreshold = 1e-12;
if minClosedLoopRealPart < branchMarginThreshold
    error('ndre:branchMargin','The branch separation margin %.3e is below %.3e.', ...
        minClosedLoopRealPart,branchMarginThreshold);
end
L = []; sylvesterScaledResidual = NaN;
if wantKernel
    L = sylvester(Dc,Ac,red.S);
    sylvesterScaledResidual = norm(Dc*L+L*Ac-red.S,'fro') / ...
        max(1,norm(red.S,'fro'));
    if sylvesterScaledResidual >= 1e-8
        error('ndre:sylvester','Centered Sylvester residual %.3e.', ...
            sylvesterScaledResidual);
    end
end
prepared=struct('Ystar',Ystar,'Ac',Ac,'Dc',Dc,'L',L, ...
    'nareBackwardError',nareBackwardError,'nareScaledResidual',nareScaledResidual, ...
    'sylvesterScaledResidual',sylvesterScaledResidual, ...
    'minClosedLoopRealPart',minClosedLoopRealPart, ...
    'branchMarginThreshold',branchMarginThreshold);
prepared.E0=-Ystar;
if isfield(red,'Y0'), prepared.E0=red.Y0-Ystar; end
end

function sol=centered_evaluate(prepared,t,wantDerivative)
Ystar=prepared.Ystar; Ac=prepared.Ac; Dc=prepared.Dc; L=prepared.L;
ell=size(Dc,1); E0=prepared.E0;
P = expm(-t*Ac); R = expm(-t*Dc);
K = L-R*L*P;
M = eye(ell)-K*E0;
middleRcond = rcond(M);
if middleRcond < 1e-14
    error('ndre:middle','Centered middle matrix is ill-conditioned.');
end
E = (P*E0/M)*R;
Y = Ystar+E;
% Differentiate the computed quotient, retaining Sylvester roundoff.
sol=prepared; sol.Y=Y; sol.middleRcond=middleRcond;
if wantDerivative
    Kdot = Dc*R*L*P + R*L*Ac*P;
    U = P*E0/M;
    sol.dY = -Ac*E-E*Dc+U*Kdot*(E0/M)*R;
end
end


function [V,W] = extended_pair(prob,targetDim)
r = size(prob.F,2);
blockWidth = 2*r;
if mod(targetDim,blockWidth) ~= 0
    error('ndre:dimension','EBA target dimension is not a full block.');
end
V = extended_side(prob.Aapply,prob.Asolve,prob.F,targetDim,r);
W = extended_side(prob.DTapply,prob.DTsolve,prob.G,targetDim,r);
end


function V = extended_side(applyM,solveM,F,targetDim,r)
blockWidth = 2*r;
V = eba_append([F,solveM(F)],zeros(size(F,1),0));
if size(V,2) ~= blockWidth
    error('ndre:initialRank','Initial EBA block rank deficient.');
end
last = V;
for j = 2:(targetDim/blockWidth)
    Z = [applyM(last(:,1:r)),solveM(last(:,r+1:end))];
    Q = eba_append(Z,V);
    if size(Q,2) ~= blockWidth
        error('ndre:blockBreakdown','EBA block deflated.');
    end
    V = [V,Q]; %#ok<AGROW>
    last = Q;
end
end


function Q = eba_append(candidate,V)
Z = candidate;
if ~isempty(V)
    Z = Z-V*(V'*Z);
    Z = Z-V*(V'*Z);
end
[Q,R] = qr(Z,0);
tol = max(size(candidate))*eps*max(norm(candidate,2),1);
if any(abs(diag(R)) <= tol)
    Q = Q(:,abs(diag(R)) > tol);
end
end

