function summary = check_results
%CHECK_RESULTS Check the recorded results in data/ against the statements
% of the paper. No solver is run.
res = ndre_paths().data;
read = @(name) readtable(fullfile(res,name),'TextType','string','VariableNamingRule','preserve');
D = read('bounds/dense_bounds.csv');
assert(height(D)==24 && all(D.status=="OK_ERROR_ESTIMATE"),'ndre:dense','Expected 24 accepted dense tests.');
assert(all(isfinite(D.projOverEstimateMax)) && max(D.projOverEstimateMax)<1, ...
    'ndre:bound','The projection error exceeds the estimate.');
assert(max(D.projOverTarget)<1,'ndre:target','The dense target was not met.');
W = read('bounds/wellposed.csv');
S = W(W.variant=="scaled",:);
assert(height(W)==16 && all(S.stationaryFailures==0) && all(S.flowFailures==0), ...
    'ndre:wellposed','A scaled projected equation failed.');
assert(all(S.maxNorm2Ystar<1) && all(S.minMargin>0),'ndre:unitBall','Scaled stationary solutions must lie in the unit ball.');
C = read('complete_solves/complete_solves.csv');
assert(height(C)==12 && numel(unique(C.n))==4,'ndre:singleState','Incomplete single-state records.');
I = C(C.method=="incremental-composed",:);
assert(all(isfinite(I.maxRelativeError)) && all(I.maxRelativeError<1e-7), ...
    'ndre:singleStateError','Single-state errors differ from Table 3.');
assert(all(I.tau(I.n<40000)==1e-7) && I.tau(I.n==40000)==1e-6, ...
    'ndre:tolerance','Unexpected tolerances.');
M = read('multistate_transport/multistate_run.csv');
I = M(M.method=="incremental-composed",:);
B = M(M.space=="EBA (published)",:);
assert(height(M)==12 && height(I)==4 && all(I.status=="OK_ERROR_ESTIMATE"), ...
    'ndre:multistate','Incomplete multistate records.');
assert(all(I.maxRelativeError<1e-7) && height(B)==4 && all(B.status=="BREAKDOWN"), ...
    'ndre:baseline','Multistate results differ from Table 4.');
E = read('evaluator_comparison/transport_evaluators.csv');
assert(height(E)==4 && all(isfinite(E.maxRelativeError)),'ndre:evaluators','Missing evaluator results.');
F = read('floor/floor_summary.csv');
assert(height(F)==4 && all(F.plateauMedian>0),'ndre:floor','Missing floor traces.');
for name = ["transport","multistate"]
    H = read("bounds/history_"+name+".csv");
    assert(all(isfinite(H.projectionError)) && all(H.projectionError<=H.estimate+1e-18), ...
        'ndre:history','The projection error exceeds the estimate.');
end
summary = struct('denseCases',height(D),'wellposedRecords',height(W), ...
    'singleStateRows',height(C),'multistateRows',height(M), ...
    'evaluatorRows',height(E),'maxProjectionErrorOverEstimate',max(D.projOverEstimateMax));
disp(summary);
fprintf('All checks of the recorded results passed.\n');
end
