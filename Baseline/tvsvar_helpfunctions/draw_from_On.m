function Qt = draw_from_On(n)
% draw from On

% Using the loop
%     [Qt,Rt] = qr(randn(n));
%     for jcol=1:n
%         if Rt(jcol,jcol)<0
%             Qt(:,jcol)=-Qt(:,jcol);
%         end
%     end
    
% Without loop    
[Qt,Rt] = qr(randn(n));
sel_Rt=(diag(Rt)<0)';
Qt(:,sel_Rt) = -Qt(:,sel_Rt);

