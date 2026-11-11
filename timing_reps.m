function r=timing_reps(r)
%TIMING_REPS Number of timed repetitions. Single run (2026-11-11): every
% solve is run once and timed once; set SINGLE=false to restore r runs
% (first run excluded as warm-up).
SINGLE=true;
if SINGLE, r=1; end
end
