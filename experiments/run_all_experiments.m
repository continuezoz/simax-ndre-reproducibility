function run_all_experiments(n)
%RUN_ALL_EXPERIMENTS All experiments of the paper, in the order in which
% they depend on each other: the dense tests and Figure 1 data (n=600 by
% default), Table 3, Table 4, Table 5 (which reuses the n=4000 reference of
% Table 3) and Figure 2. References already present in data/ are reused.
if nargin<1, n=600; end
stages={ ...
    "Tables 1-2 and stationary bounds", @() run_bound_experiments("all",n); ...
    "Figure 1 data",                    @() run_bound_history(n,1e-7); ...
    "Table 3",                          @() run_complete_solves; ...
    "Table 4",                          @() run_multistate_transport; ...
    "Table 5",                          @() run_evaluator_comparison; ...
    "Figure 2 data",                    @() run_floor_traces};
for j=1:size(stages,1)
    timer=tic; ndre_scaled(true);
    fprintf('STAGE %d: %s\n',j,stages{j,1});
    stages{j,2}();
    fprintf('STAGE %d finished in %.0f s\n',j,toc(timer));
end
end
