function dP = RBCfun(t,P,fV,fH,Lv,Lh,JV,JH,birth,death,nV,nH)
% size(-(repmat(JV(:)',[nH.^2,1])*(fV(:).*P(:)) + repmat(JH(:)',[nV.^2,1])*(fH(:).*P(:))))
% size((repmat(Lv(:)',[nH.^2,1])*P(:)+repmat(Lh(:)',[nV.^2,1])*P(:)))
% size(birth(:)-death(:))
dP = -(repmat(JV(:)',[nH.^2,1])*(fV(:).*P(:)) + repmat(JH(:)',[nV.^2,1])*(fH(:).*P(:)))...
    +(repmat(Lv(:)',[nH.^2,1])*P(:)+repmat(Lh(:)',[nV.^2,1])*P(:))+birth(:)-death(:);
t
end

