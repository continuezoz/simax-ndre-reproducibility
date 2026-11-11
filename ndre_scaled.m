function flag=ndre_scaled(newFlag)
%NDRE_SCALED Switch between the scaled equation of Section 2.1 (default) and
% the original one. With ndre_scaled true, the problem constructors return
% the scaled problem with its constant nu (ndre_setup), rksm_controller and
% the block controllers stop on the a posteriori error estimate
% (rksm_bound_controller), and results go to results_scaled/.
% ndre_scaled(false) restores the earlier behavior and results/.
persistent state
if isempty(state), state=true; end
if nargin>0, state=logical(newFlag); end
flag=state;
end
