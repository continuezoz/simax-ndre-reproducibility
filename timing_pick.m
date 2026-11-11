function s=timing_pick(seconds,stat)
%TIMING_PICK Reported time from repeated runs: with one run its time, with
% several the median (or min/max) excluding the first (warm-up) run.
if nargin<2, stat="median"; end
x=seconds(:); if numel(x)>1, x=x(2:end); end
switch string(stat)
    case "min", s=min(x);
    case "max", s=max(x);
    otherwise, s=median(x);
end
end
