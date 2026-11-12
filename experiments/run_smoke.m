function run_smoke
%RUN_SMOKE Algorithm 5.1 on two small generated problems (n=80), checked
% against the full-order solution.
ndre_scaled(true);
a = ndre_core;
problems = {a.problem('transport',80),multistate_problem(80,10,0.9999,0.1)};
for j = 1:numel(problems)
    p = problems{j};
    assert(p.nu>0,'ndre:smokeNu','Scaling did not produce a positive rate.');
    shifts = pole_order(-logspace(log10(p.shiftRange(1)),log10(p.shiftRange(2)),159),"leja");
    r = rksm_bound_controller(a,p,1e-12,1e-6,struct('tf',1,'checkDimension',0),shifts);
    require_accepted(r.status,"smoke problem "+j);
    full = a.project(p,eye(p.n),eye(p.p));
    exact = a.mobius(full,[0 1],1,0);
    approximation = r.V*r.estimate.Yhistory(:,:,end)*r.W';
    error = norm(exact.Y-approximation,'fro');
    assert(error<=r.estimate.epsilon(end)+1e-10*max(1,norm(exact.Y,'fro')), ...
        'ndre:smokeAccuracy','Error at t=1 exceeds the estimate.');
    fprintf('SMOKE problem=%d dimension=%d nu=%.6g error=%.6g estimate=%.6g\n', ...
        j,size(r.V,2),p.nu,error,r.estimate.epsilon(end));
end
fprintf('Both small problems passed.\n');
end
