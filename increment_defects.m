function z=increment_defects(red,prep,times,RA,RD)
%INCREMENT_DEFECTS Lifted-residual indicator at the stopping nodes, from
% composed increment evaluation. The derivative is taken from the error
% equation, dE = -Ac E - E Dc + E S E with E = Y - Ystar, so the reduced
% defect equals -N(Ystar) up to rounding (Proposition reduced-backward-error)
% and the Sylvester kernel L is not needed. Rows of z: [total, inner, left,
% right] relative to ||Q||_F, one row per node in times (increasing, > 0).
sol=incremental_centered_nodes(red,prep,times,"composed");
z=zeros(numel(times),4);
for j=1:numel(times)
    Y=sol.Yhistory(:,:,j); E=Y-prep.Ystar;
    dY=-prep.Ac*E-E*prep.Dc+E*red.S*E;
    inner=norm(dY+red.A*Y+Y*red.D-Y*red.S*Y-red.Q,'fro')/red.normQ;
    left=norm(RA*Y,'fro')/red.normQ; right=norm(Y*RD','fro')/red.normQ;
    z(j,:)=[sqrt(inner^2+left^2+right^2),inner,left,right];
end
end
