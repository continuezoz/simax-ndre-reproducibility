function require_accepted(status,label)
%REQUIRE_ACCEPTED Stop when the dimension search returned without passing
% its stopping test (Algorithm 1 returns REJECTED). STATUS is the status of
% rksm_controller or block_controller; an invariant space is exact and is
% accepted.
status=string(status);
if ~(startsWith(status,"OK") || status=="INVARIANT")
    error('rksm:rejected','%s: dimension search rejected (%s).',label,status);
end
end
