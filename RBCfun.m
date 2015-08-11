function dP = RBCfun(t,P,fV,fH,L,J,birth,death,nV,nH,vGridU,vGridL,hGridL,hGridU,DH,DV)
% size(-(repmat(JV(:)',[nH.^2,1])*(fV(:).*P(:)) + repmat(JH(:)',[nV.^2,1])*(fH(:).*P(:))))
% size((repmat(Lv(:)',[nH.^2,1])*P(:)+repmat(Lh(:)',[nV.^2,1])*P(:)))
% size(birth(:)-death(:))
% P=reshape(P,[nV,nH]);
% JD = J';
% keyboard
dP = J*P(:)+L*P(:)+birth(:).*sum(death(:).*P(:))/sum(birth(:))-death(:).*P(:); % .*sum(death(:).*P(:))/sum(birth(:))
% dP1 = -1*(dP<0).*min(abs(dP),P)+dP.*(dP>=0);
%  dP = -(JV.*(fV.*P) + JH.*(fH.*P))...
%      +(Lv.*P+Lh.*P)+birth-death.*sign(P);
%  dP = dP(:);
[sum(P),t]%,sum(death(:).*P(:))/sum(birth(:))]
if sum(dP)<1e-2
   return
end
%  keyboard
end

