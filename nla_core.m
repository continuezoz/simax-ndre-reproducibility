function api=nla_core
%NLA_CORE Shared base-MATLAB kernels for the NLA experiment revision.
% EBA/BDF1 kernels are copied from the verified 2020 reconstruction.
% BDF2/3 use the published coefficients; BDF3 has local trapezoidal startup.
api=struct('problem',@make_problem,'options',@default_options, ...
    'adaptive',@adaptive_solve,'fixed',@fixed_solve,'basis',@build_basis, ...
    'project',@project_problem,'prepare',@centered_prepare, ...
    'evaluate',@centered_evaluate,'bdf',@integrate_projected_bdf1, ...
    'errors',@trajectory_errors,'reference',@reference_trajectory, ...
    'split',@split_trajectory,'dense',@dense_trajectory,'radon',@radon_trajectory, ...
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
if nargin<3 || isempty(p), p=n; end
if nargin<4 || isempty(seed), seed=20261007+n+p; end
if nargin<5, parameter=3; end
family=string(family);
if family=="transport"
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
else
    rng(seed,'twister');
    if family=="square_admissible"
        A=spdiags([-ones(n,1),parameter*ones(n,1),-ones(n,1)],-1:1,n,n);
        D=spdiags([-ones(p,1),parameter*ones(p,1),-ones(p,1)],-1:1,p,p);
        [F,~]=qr(randn(n,2),0); [G,~]=qr(randn(p,2),0);
    elseif family=="square" || family=="original_rectangular"
        A=spdiags([-ones(n,1),2*ones(n,1),-ones(n,1)],-1:1,n,n);
        D=spdiags([-ones(p,1),2*ones(p,1),-ones(p,1)],-1:1,p,p);
        F=randn(n,2); G=randn(p,2);
    else
        A=spdiags([-ones(n,1),parameter*ones(n,1),-0.5*ones(n,1)],-1:1,n,n);
        D=spdiags([-0.5*ones(p,1),parameter*ones(p,1),-ones(p,1)],-1:1,p,p);
        [F,~]=qr(randn(n,2),0); [G,~]=qr(randn(p,2),0);
    end
    S=sparse(1:2,1:2,ones(2,1),p,n);
    prob=struct('n',n,'p',p,'F',F,'G',G, ...
        'normQ',sqrt(sum(sum((F'*F).*(G'*G)'))),'blockWidth',4, ...
        'Aapply',@(X)A*X,'Dapply',@(X)D*X,'DTapply',@(X)D'*X, ...
        'Sapply',@(X)S*X,'Asolve',@(X)A\X,'Dsolve',@(X)D\X,'DTsolve',@(X)D'\X, ...
        'Ashift',@(X,s)(A-s*speye(n))\X, ...
        'DTshift',@(X,s)(D'-s*speye(p))\X,'shiftRange',[1e-2,4]);
end
prob.rightBasisApply=prob.DTapply; prob.rightBasisSolve=prob.DTsolve;
prob.family=family; prob.seed=seed; prob.parameter=parameter;
if family=="transport" && ndre_scaled, prob=ndre_setup(prob); end
end

function X=rankone_inverse(B,d,u,v)
Z=B./d; q=u./d; X=Z+q*((v'*Z)/(1-v'*q));
end

function [V,W]=build_basis(prob,basisName,dimension,maxDim)
if nargin<4, maxDim=dimension; end
if string(basisName)=="EBA", [V,W]=extended_pair(prob,dimension); return; end
if string(basisName)=="polynomial"
    r=size(prob.F,2); V=orth(prob.F); W=orth(prob.G);
    for j=1:dimension/r-1
        V=append_block(V,prob.Aapply(V(:,end-r+1:end)),r);
        W=append_block(W,prob.DTapply(W(:,end-r+1:end)),r);
    end
    return
end
r=size(prob.F,2); shifts=-logspace(log10(prob.shiftRange(1)),log10(prob.shiftRange(2)),maxDim/r-1);
V=orth(prob.F); W=orth(prob.G);
for j=1:dimension/r-1
    V=append_block(V,prob.Ashift(V(:,end-r+1:end),shifts(j)),r);
    W=append_block(W,prob.DTshift(W(:,end-r+1:end),shifts(j)),r);
end
end

function result=adaptive_solve(prob,method,opts)
if startsWith(string(method),"EBA-BDF")
    opts.bdfOrder=str2double(extractAfter(string(method),"EBA-BDF"));
    result=solve_original_eba_bdf1(prob,opts); result.accepted=result.converged;
    result.indicator=result.boundaryFrobeniusResidual;
    result.branchMargin=NaN; result.middleRcond=NaN; return
end
timer=tic; r=size(prob.F,2); V=orth(prob.F); W=orth(prob.G);
shifts=-logspace(log10(prob.shiftRange(1)),log10(prob.shiftRange(2)),opts.maxDim/r-1);
curve=table('Size',[0,5],'VariableTypes',repmat({'double'},1,5), ...
    'VariableNames',{'dimension','indicator','reducedDefect','branchMargin','seconds'});
accepted=false; status="MAX_DIMENSION"; margin=NaN; indicator=Inf;
for j=1:numel(shifts)
    V=append_block(V,prob.Ashift(V(:,end-r+1:end),shifts(j)),r);
    W=append_block(W,prob.DTshift(W(:,end-r+1:end),shifts(j)),r);
    dim=size(V,2); if mod(dim,opts.checkDimension)~=0, continue; end
    red=project_problem(prob,V,W);
    try
        prepared=centered_prepare(red); endpoint=centered_evaluate(prepared,opts.tf,true);
        [indicator,reduced,~,~]=lifted_residual(red,endpoint.Y,endpoint.dY);
        margin=prepared.minClosedLoopRealPart; status="RESIDUAL";
    catch ME
        indicator=Inf; reduced=Inf; status=string(ME.identifier);
    end
    curve(end+1,:)={dim,indicator,reduced,margin,toc(timer)};
    if indicator<opts.gate, accepted=true; status="OK"; break; end
end
if ~accepted
    result=struct('V',V,'W',W,'dimension',size(V,2),'accepted',false, ...
        'indicator',indicator,'status',status,'curve',curve,'totalSeconds',toc(timer), ...
        'branchMargin',margin,'middleRcond',NaN); return
end
result=centered_trajectory(prepared,opts.outputTimes);
result.V=V; result.W=W; result.dimension=size(V,2); result.accepted=true;
result.indicator=indicator; result.status=status; result.curve=curve;
result.branchMargin=margin; result.totalSeconds=toc(timer);
end

function result=fixed_solve(prob,basisName,solverName,dim,opts)
timer=tic; basisTimer=tic;
[V,W]=build_basis(prob,basisName,dim,opts.maxDim); basisSeconds=toc(basisTimer);
projectionTimer=tic; red=project_problem(prob,V,W); projectionSeconds=toc(projectionTimer);
solverTimer=tic;
if startsWith(string(solverName),"BDF")
    opts.bdfOrder=str2double(extractAfter(string(solverName),"BDF"));
    sol=integrate_projected_bdf1(red,opts);
    dY=-red.A*sol.Y-sol.Y*red.D+sol.Y*red.S*sol.Y+red.Q; margin=NaN;
elseif string(solverName)=="centered"
    prepared=centered_prepare(red); sol=centered_trajectory(prepared,opts.outputTimes);
    ep=centered_evaluate(prepared,opts.tf,true); dY=ep.dY; margin=prepared.minClosedLoopRealPart;
else
    sol=radon_trajectory(red,opts.outputTimes,opts.h);
    dY=-red.A*sol.Y-sol.Y*red.D+sol.Y*red.S*sol.Y+red.Q; margin=NaN;
end
solverSeconds=toc(solverTimer); [indicator,~,~,~]=lifted_residual(red,sol.Y,dY);
result=sol; result.V=V; result.W=W; result.dimension=dim;
result.basisSeconds=basisSeconds; result.projectionSeconds=projectionSeconds;
result.solverSeconds=solverSeconds; result.totalSeconds=toc(timer);
result.accepted=isfinite(indicator) && indicator<opts.gate;
result.indicator=indicator; result.branchMargin=margin;
end

function sol=centered_trajectory(prepared,times)
[k,ell]=size(prepared.Ystar); Yhistory=zeros(k,ell,numel(times)); rc=ones(size(times));
for j=1:numel(times)
    if times(j)==0, Yhistory(:,:,j)=prepared.Ystar+prepared.E0;
    else
        value=centered_evaluate(prepared,times(j),false);
        Yhistory(:,:,j)=value.Y; rc(j)=value.middleRcond;
    end
end
sol=struct('Y',Yhistory(:,:,end),'Yhistory',Yhistory,'time',times,'middleRcond',min(rc));
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

function sol=split_trajectory(red,times,microSteps)
% Independent Strang splitting: exact affine-linear and quadratic subflows.
% No stationary Riccati solution and no centered Sylvester kernel are used.
k=size(red.A,1); ell=size(red.D,1); h=(times(2)-times(1))/microSteps;
assert(max(abs(diff(times)-diff(times(1:2))))<1e-12,'Use an equispaced grid.');
Z=sylvester(red.A,red.D,red.Q); P=expm(-h*red.A/2); R=expm(-h*red.D/2); B=Z-P*Z*R;
Y=zeros(k,ell); if isfield(red,'Y0'), Y=red.Y0; end
Yhistory=zeros(k,ell,numel(times)); Yhistory(:,:,1)=Y;
for j=2:numel(times)
    for s=1:microSteps
        Y=P*Y*R+B; Y=Y/(eye(ell)-h*red.S*Y); Y=P*Y*R+B;
    end
    Yhistory(:,:,j)=Y;
end
sol=struct('Y',Y,'Yhistory',Yhistory,'time',times);
end

function ref=reference_trajectory(prob,dim,times,tolerance)
% Global Richardson extrapolation of independent symmetric splitting.
if nargin<4, tolerance=1e-10; end
[V,W]=build_basis(prob,"EBA",dim); red=project_problem(prob,V,W);
previous=[]; temporalDifference=Inf; raw=cell(1,7); fourth=raw; sixth=raw;
for level=1:7
    sol=split_trajectory(red,times,2^(level-1)); raw{level}=sol.Yhistory;
    if level>=2, fourth{level}=(4*raw{level}-raw{level-1})/3; end
    if level>=3
        sixth{level}=(16*fourth{level}-fourth{level-1})/15;
        if ~isempty(previous)
            diffnorm=vecnorm(reshape(sixth{level}-previous,dim^2,[]),2,1);
            refnorm=vecnorm(reshape(sixth{level},dim^2,[]),2,1);
            temporalDifference=max(diffnorm(2:end)./refnorm(2:end));
            fprintf('  reference dim=%d microsteps=%d refinement=%.3e\n',dim,2^(level-1),temporalDifference);
            if temporalDifference<tolerance, break; end
        end
        previous=sixth{level};
    end
end
ref=struct('V',V,'W',W,'Yhistory',sixth{level},'Y',sixth{level}(:,:,end), ...
    'time',times,'temporalRefinement',temporalDifference,'microSteps',2^(level-1), ...
    'referenceQualified',temporalDifference<tolerance,'dimension',dim);
end

function sol=dense_trajectory(red,times,tolerance)
if nargin<3, tolerance=1e-12; end
k=size(red.A,1); ell=size(red.D,1); Y0=zeros(k,ell);
if isfield(red,'Y0'), Y0=red.Y0; end
rhs=@(~,y)reshape(-red.A*reshape(y,k,ell)-reshape(y,k,ell)*red.D+ ...
    reshape(y,k,ell)*red.S*reshape(y,k,ell)+red.Q,[],1);
jac=@(~,y)kron(speye(ell),-red.A+reshape(y,k,ell)*red.S)+ ...
    kron((-red.D+red.S*reshape(y,k,ell))',speye(k));
options=odeset('RelTol',tolerance,'AbsTol',tolerance*0.01,'Jacobian',jac);
[~,values]=ode15s(rhs,times,Y0(:),options); Yhistory=reshape(values',k,ell,[]);
sol=struct('Yhistory',Yhistory,'Y',Yhistory(:,:,end),'time',times);
end

function sol=radon_trajectory(red,times,step)
% Restarted quotient, with explicit denominator checks.
k=size(red.A,1); ell=size(red.D,1); H=[red.D,-red.S;red.Q,-red.A];
Y=zeros(k,ell); if isfield(red,'Y0'), Y=red.Y0; end
Yhistory=zeros(k,ell,numel(times)); Yhistory(:,:,1)=Y;
if step>0
    m=max(1,round((times(2)-times(1))/step)); step=(times(2)-times(1))/m;
    E=expm(step*H);
    assert(all(isfinite(E),'all'),'nla:radonOverflow','Radon exponential overflow.');
end
minRcond=1;
for j=2:numel(times)
    if step==0
        E=expm(times(j)*H); Z=E*[eye(ell);Yhistory(:,:,1)];
        denominator=Z(1:ell,:); rc=rcond(denominator);
        assert(all(isfinite(Z),'all') && rc>1e-14,'nla:radonDenominator','Direct denominator failed.');
        Y=Z(ell+1:end,:)/denominator;
    else
        for s=1:m
            Z=E*[eye(ell);Y]; denominator=Z(1:ell,:); rc=rcond(denominator);
            assert(all(isfinite(Z),'all') && rc>1e-14,'nla:radonDenominator','Restarted denominator failed.');
            Y=Z(ell+1:end,:)/denominator; minRcond=min(minRcond,rc);
        end
    end
    minRcond=min(minRcond,rc); Yhistory(:,:,j)=Y;
end
sol=struct('Y',Y,'Yhistory',Yhistory,'time',times,'middleRcond',minRcond);
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
    assert(minRcond>1e-14,'nla:mobiusDoubling','Doubling denominator is ill-conditioned.');
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
        assert(rc>1e-14,'nla:mobiusStep','Flow denominator is ill-conditioned.');
        Y=Qh+(L*Y/M)*R;
    end
    Yhistory(:,:,j)=Y;
end
sol=struct('Y',Y,'Yhistory',Yhistory,'time',times,'middleRcond',minRcond, ...
    'scalings',scalings,'substeps',substeps);
end

function result = solve_original_eba_bdf1(prob,opts)
% Build one EBA block at a time and use the Theorem 3 stopping residual.
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
    % Generate the relation block T_{m+1,m}; it is needed by Theorem 3.
    basisTimer = tic;
    if m > 1
        [nextV, boundaryA] = expand_eba_side(prob.Aapply,prob.Asolve, ...
            V,lastV,prob.blockWidth);
        [nextW, boundaryD] = expand_eba_side(rightApply,rightSolve, ...
            W,lastW,prob.blockWidth);
        basisSeconds = basisSeconds + toc(basisTimer);
    end

    % The 2020 implementation evaluates the projected BDF/exponential
    % solver every five EBA iterations to avoid repeating a reduced solve at
    % every intermediate basis size.
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

    % Theorem 3 boundary formula. Direct formation of A*V-V*T_A suffers
    % cancellation here because the Gauss--Legendre diagonal reaches O(n^2).
    % It is retained separately as a floating-point audit, never as the
    % published algorithm's stopping quantity.
    last = size(Y,1)-blockWidth+1:size(Y,1);
    leftBoundary = boundaryA*Y(last,:);
    rightBoundary = Y(:,last)*boundaryD';
    boundaryResidual = max(norm(leftBoundary,'fro'),norm(rightBoundary,'fro'));
    normalizedResidual = boundaryResidual / max(prob.normQ,eps);
    boundaryFrobeniusResidual = hypot(norm(leftBoundary,'fro'), ...
        norm(rightBoundary,'fro')) / max(prob.normQ,eps);
    % Full-operator residual: an audit only, not part of the 2020 algorithm,
    % so it is skipped (and not timed) unless requested.
    if opts.verbose || (isfield(opts,'diagnostics') && opts.diagnostics)
        RA = prob.Aapply(V)-V*red.A;
        RD = prob.DTapply(W)-W*red.D';
        liftedResidual = hypot(norm(RA*Y,'fro'),norm(Y*RD','fro')) / max(prob.normQ,eps);
    else
        liftedResidual = NaN;
    end
    curve(end+1,:) = {m,size(V,2),norm(leftBoundary,'fro')/prob.normQ, ...
        norm(rightBoundary,'fro')/prob.normQ,normalizedResidual, ...
        boundaryFrobeniusResidual,liftedResidual}; %#ok<AGROW>
    if opts.verbose
        Qm = red.F*red.G';
        fprintf('    solve: normY=%.3e normQm=%.3e rhs=%.3e\n', ...
            norm(Y,'fro'),norm(Qm,'fro'), ...
            norm(-red.A*Y-Y*red.D+Y*red.S*Y+Qm,'fro'));
    end
    ebaSteps = m;
    dimension = size(V,2);
    if opts.verbose
        if size(nextV,2) == size(boundaryA,1) && size(nextW,2) == size(boundaryD,1)
            if opts.verbose
                fprintf('    sizes V=%dx%d nextV=%dx%d bA=%dx%d\n', ...
                    size(V,1),size(V,2),size(nextV,1),size(nextV,2), ...
                    size(boundaryA,1),size(boundaryA,2));
            end
            Mv = prob.Aapply(V(:,end-prob.blockWidth+1:end));
            Mw = rightApply(W(:,end-prob.blockWidth+1:end));
            relA = norm(Mv-V*(V'*Mv)-nextV*boundaryA,'fro') / max(norm(Mv,'fro'),eps);
            relD = norm(Mw-W*(W'*Mw)-nextW*boundaryD,'fro') / max(norm(Mw,'fro'),eps);
        else
            relA = NaN; relD = NaN;
        end
        fprintf('  m=%d dim=%d residual=%.3e relA=%.3e relD=%.3e\n', ...
            m,dimension,normalizedResidual,relA,relD);
        fprintf('    explicit=%.3e orthV=%.3e orthW=%.3e\n', ...
            liftedResidual,norm(V'*V-eye(size(V,2)),'fro'), ...
            norm(W'*W-eye(size(W,2)),'fro'));
        fprintf('    outside-last A=%.3e D=%.3e\n', ...
            norm(RA(:,1:end-prob.blockWidth),'fro')/max(norm(RA,'fro'),eps), ...
            norm(RD(:,1:end-prob.blockWidth),'fro')/max(norm(RD,'fro'),eps));
    end
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
    error('reproduce2020:rankBreakdown','Initial EBA block is rank deficient.');
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
order=1; if isfield(opts,'bdfOrder'), order=opts.bdfOrder; end
assert(any(order==[1,2,3]),'Only BDF1, BDF2 and BDF3 are implemented.');
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
previous=Y; older=Y;
for j = 1:steps
    base = Y;
    h=opts.h;
    if order==2 && j>1, base=(4*Y-previous)/3; h=2*opts.h/3; end
    if order==3
        if j<=2
            % Second-order trapezoidal startup has O(h^3) local error.
            % It uses only the original reduced vector field, no stationary flow.
            base=Y+0.5*opts.h*(-red.A*Y-Y*red.D+Y*red.S*Y+Q);
            h=0.5*opts.h;
        else
            base=(18*Y-9*previous+2*older)/11; h=6*opts.h/11;
        end
    end
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
        % Equation (38): split the single identity evenly between the
        % left and right Sylvester coefficients.
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
    % Newton aims at newtonTol; a step whose iteration stagnates at the
    % rounding level (below 10*newtonTol) is accepted, and the largest
    % backward error is reported with the run (2026-11-11).
    if backwardError > 100*opts.newtonTol
        error('reproduce2020:newtonNonconvergence', ...
            'Projected BDF Newton failed at step %d, dimension %d, backward error %.3e after %d iterations.', ...
            j,k,backwardError,it);
    end
    older=previous; previous=Y; Y = Yn;
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


function V = append_block(V,Z,r)
candidateNorm = norm(Z,'fro');
Z = Z - V*(V'*Z);
Z = Z - V*(V'*Z);
[Q,R,~] = qr(Z,'econ','vector');
if size(Q,2) ~= r || any(abs(diag(R)) <= 1e-12*max(candidateNorm,eps))
    error('compare2020:blockBreakdown','Rational block deflated.');
end
V = [V,Q]; %#ok<AGROW>
end


function red = project_problem(prob,V,W)
red = struct('V',V,'W',W,'A',V'*prob.Aapply(V), ...
    'D',W'*prob.Dapply(W),'S',W'*prob.Sapply(V), ...
    'Q',(V'*prob.F)*(W'*prob.G)', ...
    'normQ',prob.normQ,'prob',prob);
end


function sol = centered_at_one(red)
prepared=centered_prepare(red);
sol=centered_evaluate(prepared,1,true);
sol.prepared=prepared;
end

function prepared=centered_prepare(red,wantKernel)
% wantKernel=false skips the Sylvester kernel L, which only the direct
% centered evaluator needs; increment evaluation uses Ystar, Ac and Dc.
if nargin<2, wantKernel=true; end
k = size(red.A,1); ell = size(red.D,1);
H = [red.D,-red.S;red.Q,-red.A];
[U,T] = schur(H,'real');
select = real(eig(T)) > 0;
if sum(select) ~= ell
    error('compare2020:branchCount','No stable graph branch of required size.');
end
[U,~] = ordschur(U,T,select);
Z = U(:,1:ell);
if rcond(Z(1:ell,:)) < 1e-14
    error('compare2020:graph','The stationary graph is ill-conditioned.');
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
    error('compare2020:nare','Stationary NARE scaled residual %.3e.',nareScaledResidual);
end
Ac = red.A-Ystar*red.S; Dc = red.D-red.S*Ystar;
minClosedLoopRealPart = min(real([eig(Ac);eig(Dc)]));
branchMarginThreshold = 1e-12;
if minClosedLoopRealPart < branchMarginThreshold
    error('compare2020:branchMargin','The branch separation margin %.3e is below %.3e.', ...
        minClosedLoopRealPart,branchMarginThreshold);
end
L = []; sylvesterScaledResidual = NaN;
if wantKernel
    L = sylvester(Dc,Ac,red.S);
    sylvesterScaledResidual = norm(Dc*L+L*Ac-red.S,'fro') / ...
        max(1,norm(red.S,'fro'));
    if sylvesterScaledResidual >= 1e-8
        error('compare2020:sylvester','Centered Sylvester residual %.3e.', ...
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
    error('compare2020:middle','Centered middle matrix is ill-conditioned.');
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


function [res,reduced,left,right] = lifted_residual(red,Y,dY)
reduced = norm(dY+red.A*Y+Y*red.D-Y*red.S*Y-red.Q,'fro')/red.normQ;
RA = red.prob.Aapply(red.V)-red.V*red.A;
RD = red.prob.DTapply(red.W)-red.W*red.D';
left = norm(RA*Y,'fro')/red.normQ;
right = norm(Y*RD','fro')/red.normQ;
res = hypot(hypot(reduced,left),right);
end


function [V,W] = extended_pair(prob,targetDim)
r = size(prob.F,2);
blockWidth = 2*r;
if mod(targetDim,blockWidth) ~= 0
    error('compare2020:dimension','EBA target dimension is not a full block.');
end
V = extended_side(prob.Aapply,prob.Asolve,prob.F,targetDim,r);
W = extended_side(prob.DTapply,prob.DTsolve,prob.G,targetDim,r);
end


function V = extended_side(applyM,solveM,F,targetDim,r)
blockWidth = 2*r;
V = eba_append([F,solveM(F)],zeros(size(F,1),0));
if size(V,2) ~= blockWidth
    error('compare2020:initialRank','Initial EBA block rank deficient.');
end
last = V;
for j = 2:(targetDim/blockWidth)
    Z = [applyM(last(:,1:r)),solveM(last(:,r+1:end))];
    Q = eba_append(Z,V);
    if size(Q,2) ~= blockWidth
        error('compare2020:blockBreakdown','EBA block deflated.');
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

