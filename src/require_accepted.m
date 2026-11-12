function require_accepted(status,label)
%REQUIRE_ACCEPTED Stop when the dimension search returned without passing
% its stopping test (the poles were exhausted). STATUS is the status of
% rksm_bound_controller; an invariant space is exact and is accepted.
status=string(status);
if ~(startsWith(status,"OK") || status=="INVARIANT")
    error('rksm:rejected','%s: dimension search rejected (%s).',label,status);
end
end
