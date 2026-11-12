function s=pole_order(s,order)
%POLE_ORDER Reorder a fixed pole set so that every prefix samples its range.
% "monotone" keeps the given order, and "leja" uses discrete Leja points of
% log|s| starting at s(1).
n=numel(s);
switch string(order)
    case "monotone"
    case "leja"
        x=log(abs(s)); idx=zeros(1,n); idx(1)=1; used=false(1,n); used(1)=true;
        logdist=log(abs(x-x(1))+realmin);
        for j=2:n
            logdist(used)=-Inf; [~,k]=max(logdist); idx(j)=k; used(k)=true;
            logdist=logdist+log(abs(x-x(k))+realmin);
        end
        s=s(idx);
    otherwise
        error('pole:order','Unknown pole order %s.',order);
end
end
