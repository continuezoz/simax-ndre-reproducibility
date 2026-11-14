function result=eba_deflated_bdf(a,prob,opts)
%EBA_DEFLATED_BDF EBA--BDF1 of Angelova, Hached and Jbilou (2020) with
% deflation of the extended block Arnoldi basis. The published method
% (ndre_core, adaptive 'EBA-BDF1') stops when an extended block loses rank,
% as it does for the multistate model. Here the positive (M*) and negative
% (M^{-1}*) block sequences are orthogonalized separately, dependent
% columns (pivoted QR diagonal below 1e-10 times the block norm, as for the
% rational spaces) are removed, and only the retained columns are expanded.
% Every opts.checkInterval steps the projected equation is integrated by
% backward Euler (a.bdf) and the stopping quantity of the published method,
%   max(||R_A*Y||_F,||Y*R_D'||_F)/||Q||_F  at t_f,
% is formed from the Arnoldi residuals R_A=(I-VV')AV, R_D=(I-WW')D'W.
rightApply=prob.DTapply; rightSolve=prob.DTsolve;
if isfield(prob,'rightBasisApply'), rightApply=prob.rightBasisApply; rightSolve=prob.rightBasisSolve; end
timer=tic;
[V,Pv,Nv]=start_side(prob.Asolve,prob.F);
[W,Pw,Nw]=start_side(rightSolve,prob.G);
converged=false; status="MAX_STEPS"; residual=NaN; sol=[];
for m=1:opts.maxEBASteps
    if mod(m,opts.checkInterval)==0 || m==opts.maxEBASteps
        red=a.project(prob,V,W); sol=a.bdf(red,opts);
        RA=prob.Aapply(V); RA=RA-V*(V'*RA); RA=RA-V*(V'*RA);
        RD=prob.DTapply(W); RD=RD-W*(W'*RD); RD=RD-W*(W'*RD);
        residual=max(norm(RA*sol.Y,'fro'),norm(sol.Y*RD','fro'))/prob.normQ;
        if residual<opts.tol, converged=true; status="OK"; break; end
    end
    if size(V,2)>=opts.maxDim, status="MAX_DIMENSION"; break; end
    [V,Pv,Nv]=expand_side(prob.Aapply,prob.Asolve,V,Pv,Nv);
    [W,Pw,Nw]=expand_side(rightApply,rightSolve,W,Pw,Nw);
    if isempty(Pv) && isempty(Nv) && isempty(Pw) && isempty(Nw), status="INVARIANT"; break; end
end
if ~converged, red=a.project(prob,V,W); sol=a.bdf(red,opts); end
result=struct('V',V,'W',W,'Y',sol.Y,'converged',converged,'status',status, ...
    'residual',residual,'ebaSteps',m,'seconds',toc(timer));
if isfield(sol,'Yhistory'), result.Yhistory=sol.Yhistory; result.time=sol.time; end
end

function [V,P,N]=start_side(solveM,F)
V=zeros(size(F,1),0);
[V,P]=deflate(V,F); [V,N]=deflate(V,solveM(F));
assert(~isempty(V),'ndre:ebaStart','Empty starting block.');
end

function [V,P,N]=expand_side(applyM,solveM,V,P,N)
if ~isempty(P), [V,P]=deflate(V,applyM(P)); end
if ~isempty(N), [V,N]=deflate(V,solveM(N)); end
end

function [V,Q]=deflate(V,Z)
scale=max(norm(Z,'fro'),eps); Z=Z-V*(V'*Z); Z=Z-V*(V'*Z);
[q,R,~]=qr(Z,0); keep=abs(diag(R))>1e-10*scale; Q=q(:,keep); V=[V,Q];
end
