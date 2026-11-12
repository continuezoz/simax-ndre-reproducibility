function flag=ndre_scaled(newFlag)
%NDRE_SCALED Switch between the scaled equation of section 2.1 (default) and
% the original one. With ndre_scaled true, the problem constructors return
% the scaled problem with its constant nu (ndre_setup). ndre_scaled(false)
% gives the unscaled constructors; run_bound_experiments and
% run_bound_history use them and call ndre_setup explicitly, since the
% unscaled problem enters the comparison of Table 1.
persistent state
if isempty(state), state=true; end
if nargin>0, state=logical(newFlag); end
flag=state;
end
