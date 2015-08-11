clear
clc
close all
load ('RETdist.mat')
%%
x = [0.05 26 16 2.3 0.025 72]; %set containing [alpha, betaV, betaH, DV, DH, Vc];
% x = [0.05 26 16 0.014 0.0014 0.8]; %set containing [alpha, betaV, betaH, DV, DH, Vc];
alpha = x(1);
betaV = x(2);
betaH = x(3);
DV = x(4);
DH = x(5);
Vc = x(6);
%%
MCV = 90; %fl
MCH = 30; %g/dl
vGridU = 140;
hGridU = 50;
vGridL = 30;
hGridL = 10;
nV = 62;
nH = 23;

disV = linspace(vGridL,vGridU,nV); %1-D grid in Vol
disH = linspace(hGridU,hGridL,nH); %1-D grid in Hgb
[Vgrid,Hgrid]=meshgrid(disV,disH); %Creating a 2-D grid for V and H
% Hgrid = flipud(Hgrid);
%%
P0(:,1) = P0(:,1);
P0(:,2) = P0(:,2).*P0(:,1)/MCV;
P0 = fliplr(P0);
P0fit = hist3(P0,{fliplr(disH),disV});
% P1 = hist3(P0,{disV,fliplr(disH)});
P0fit = flipud(P0fit)./sum(P0fit(:));
% Birth
birth = P0fit;%mvnpdf([Vgrid(:),Hgrid(:)],[1.3,1.2],[30,10]); %this should be the RET distribution
birth=reshape(birth,length(disH),length(disV));
% Death
VVbar = Vgrid.*MCV;
HHbar = Hgrid.*MCH;
theta = atan(MCH/MCV)*ones(length(disH),length(disV))-atan(HHbar./VVbar);
Delta = 100 * ( cos(theta) .* sqrt(VVbar.^2+HHbar.^2)-Vc.* sqrt(MCH^2.*ones(length(disH),length(disV))+...
    MCV^2.*ones(length(disH),length(disV))) )./ (sqrt(MCH^2.*ones(length(disH),length(disV))+...
    MCV^2.*ones(length(disH),length(disV))));
death = 1./(1+exp(Delta)); %zeros(nV,nH);
% death = -(sign(Delta)-1)/2;
% Jacobian
DisMat = -1*eye(length(disH),length(disV));
dummy = [1*eye(length(disH)-1,length(disV)-1);zeros(1,length(disV)-1)];
dummy = [zeros(length(disH),1),dummy];
% J = (DisMat+dummy);%sparse((nV-1)/(vGridU-vGridL)*(DisMat+dummy));
% JH = ((nH-1)/(hGridU-hGridL)*(DisMat+dummy));%sparse((nH-1)/(hGridU-hGridL)*(DisMat+dummy));
Po = birth(:);
% clear J
% load('/Users/anweshachaudhury/Desktop/Anwesha research/RBCcontrol/JJvJh.mat')
fV = alpha.*(MCV./50).*exp(betaV.*((Vgrid./MCV)-(Hgrid./MCH))).*((nV-1)/(vGridU-vGridL));
fH = alpha.*(MCH./50).*exp(betaH.*((Hgrid./MCH)-(Vgrid./MCV))).*((nH-1)/(hGridU-hGridL));
% fV = -1*mVDot;
% fH = -1*mHDot;
%Laplacian
dummy1 = [zeros(1,length(disV)-1);1*eye(length(disH)-1,length(disV)-1)];
dummy1 = [dummy1,zeros(length(disH),1)];
% Lv = (DV.*(((nV-1)/(vGridU-vGridL)).^2).*(2*DisMat + dummy1 + dummy));
% L = (2*DisMat + dummy1 + dummy);
I1 = eye(nV,nV);
I2 = eye(nH,nH);
I = eye(nH,nV);
E = sparse(2:nV,1:nV-1,1,nV,nV);
Ep = sparse(2:nH,1:nH-1,1,nH,nH);
D1 = E+E'-2*I1;
D2 = Ep+Ep'-2*I2;
L1 = kron(I1,D2.*DH.*((nH-1)/(hGridU-hGridL)).^2)+kron(D1.*DV.*((nV-1)/(vGridU-vGridL)).^2,I2);
% J1j = kron(I,J'.*(nH-1)./(hGridU-hGridL)).*repmat(fH(:)',[nH.*nV,1])+kron(J.*(nV-1)./(vGridU-vGridL),I').*repmat(fV(:)',[nH.*nV,1]);
J1 = -1*eye(nH*nV);
% Jj = -1*eye(nV,nH);
J11 = [eye(nV*nH-1);zeros(1,nV*nH-1)];
J11 = [zeros(nV*nH,1),J11];
% JJ = J';
% for i = 1:nH-1
%  JJ =  blkdiag(JJ,J);
% end

Jh = (J1+J11').*repmat(fH(:)',[nH.*nV,1]);
for i = 1:nV-1
   Jh(nH*i+1,nH*i) = 0;
end
J12 = zeros(nV*nH);
for i = 1:nH*nV-nH
J12(i,nH+i) = 1;  
end

Jv = (J12+J1).*repmat(fV(:)',[nH.*nV,1]);
J1 = Jv+Jh;
%Solution using inbuilt MATLAB function 
options = odeset('AbsTol',1e-6);
sol = ode15s(@(t,P) RBCfun(t,P,fV,fH,L1,J1,birth,death,nV,nH,vGridU,vGridL,hGridL,hGridU,DH,DV),[0 1000],Po,options);
%% Plotting
figure
for i = 1:numel(sol.y(1,:))
  P = (reshape(sol.y(:,i),nH,nV));
   surf(disV,disH,P)
%    axis tight manual
 axis([vGridL vGridU hGridL hGridU 0 max(max(sol.y))]);
  F(i) = getframe;
%   keyboardf
end
% movie(F)