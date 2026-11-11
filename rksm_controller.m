function r=rksm_controller(a,p,strategy,gate,o,guard,shifts)
%RKSM_CONTROLLER Two-sided rational Krylov expansion with sampled stopping.
% Expands both bases with the given ordered poles, one block of width
% size(p.F,2) per pole, and checks the lifted residual every
% o.checkDimension columns, at o.tf only ("endpoint") or additionally on the
% nodes in guard ("transient"). The optional o.initialScale sets X0 = s*Q,
% whose factors lie in the starting spaces.
if isfield(p,'nu')
    % Scaled problem: Algorithm 1 with the a posteriori error estimate; the
    % residual strategy and the guard nodes are not used, gate is tau.
    r=rksm_bound_controller(a,p,1e-12,gate,o,shifts); r.gate=gate;
    r.indicator=Inf; if height(r.trace)>0, r.indicator=r.trace.ratio(end); end; r.failedDimensions=zeros(0,1); r.failureReasons=strings(0,1);
    return
end
timer=tic; V=orth(p.F); W=orth(p.G); b=size(p.F,2);
if ~isfield(o,'initialScale'), o.initialScale=0; end
trace=table('Size',[0 5],'VariableTypes',repmat({'double'},1,5), ...
    'VariableNames',{'dimension','endpointResidual','monitorResidual','peakTime','seconds'});
accepted=false; status="MAX_DIMENSION"; indicator=Inf; failed=zeros(0,1); reasons=strings(0,1);
for shift=shifts
    V=append(V,p.Ashift(V(:,end-b+1:end),shift));
    W=append(W,p.DTshift(W(:,end-b+1:end),shift));
    dim=size(V,2); if mod(dim,o.checkDimension)~=0, continue; end
    red=a.project(p,V,W); red.Y0=o.initialScale*red.Q;
    % Branch preparation without the Sylvester kernel; states at the
    % stopping nodes come from composed increment evaluation.
    try, prep=a.prepare(red,false);
    catch ME, failed(end+1,1)=dim; reasons(end+1,1)=string(ME.identifier); continue; %#ok<AGROW>
    end
    [~,RA]=qr(p.Aapply(V)-V*red.A,0);
    [~,RD]=qr(p.DTapply(W)-W*red.D',0);
    e=increment_defects(red,prep,o.tf,RA,RD);
    indicator=e(1); peakTime=o.tf;
    if indicator<gate && strategy~="endpoint"
        nodes=unique([guard(guard>0),o.tf]);
        z=increment_defects(red,prep,nodes,RA,RD); [peak,j]=max(z(:,1));
        if peak>indicator, indicator=peak; peakTime=nodes(j); end
    end
    trace(end+1,:)={dim,e(1),indicator,peakTime,toc(timer)}; %#ok<AGROW>
    fprintf('CONTROLLER %s gate=%.0e m=%d endpoint=%.3e monitor=%.3e peak=%.3e\n', ...
        strategy,gate,dim,e(1),indicator,peakTime);
    if indicator<gate, accepted=true; status="OK_SAMPLED_RESIDUAL"; break; end
end
r=struct('V',V,'W',W,'accepted',accepted,'status',status,'trace',trace, ...
    'indicator',indicator,'gate',gate,'seconds',toc(timer), ...
    'failedDimensions',failed,'failureReasons',reasons);
end

function Q=append(V,Z)
scale=norm(Z,'fro'); Z=Z-V*(V'*Z); Z=Z-V*(V'*Z);
[q,R]=qr(Z,0);
assert(all(abs(diag(R))>1e-12*max(scale,eps)),'pilot:breakdown','Rational basis deflated.');
Q=[V,q];
end
