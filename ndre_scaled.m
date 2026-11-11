function flag=ndre_scaled(newFlag)
%NDRE_SCALED Switch between the scaled equation of Section 2.1 (default) and
% the original one. With ndre_scaled true, the problem constructors return
% the scaled problem with its constant nu (ndre_setup), and results go to
% results_scaled/. ndre_scaled(false) gives the unscaled constructors, which
% run_bound_experiments and run_bound_history use before calling
% ndre_setup explicitly (the unscaled problem enters only the
% well-posedness comparison).
persistent state
if isempty(state), state=true; end
if nargin>0, state=logical(newFlag); end
flag=state;
end
