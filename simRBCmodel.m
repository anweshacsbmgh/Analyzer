clear
clc
close all
x = [0.06 25 15 0.012 1e-3 0.8]; %set containing [alpha, betaV, betaH, DV, DH, Vc];
alpha = x(1);
betaV = x(2);
betaH = x(3);
DV = x(4);
DH = x(5);
Vc = x(6);
vGridU = 200;
hGridU = 50;
vGridL = 40;
hGridL = 20;
nV = 50;
nH = 50;
MCV = 100; %fl
MCH = 35; %g/dl
disV = linspace(vGridL,vGridU,nV)./MCV; %1-D grid in Vol
disH = linspace(hGridL,hGridU,nH)./MCH; %1-D grid in Hgb
[Vgrid,Hgrid]=meshgrid(disV,disH); %Creating a 2-D grid for V and H

% Birth
birth = mvnpdf([Vgrid(:),Hgrid(:)],[1.3*MCV,1.2*MCH],[30,10]); %this should be the RET distribution
birth=reshape(birth,length(disV),length(disH));
% Death
VVbar = MCV.*Vgrid;
HHbar = MCH.*Hgrid;
theta = atan(MCH/MCV)*ones(length(disV),length(disH))-atan(HHbar./VVbar);
Delta = 100 * ( cos(theta) .* sqrt(VVbar.^2+HHbar.^2)-Vc.* sqrt(MCH^2.*ones(length(disV),length(disH))+...
    MCV^2.*ones(length(disV),length(disH))) )./ (Vc.*sqrt(MCH^2+MCV^2));
death = 1./(1+exp(Delta));
% Jacobian
DisMat = -1*eye(length(disV),length(disH));
dummy = [1*eye(length(disV)-1,length(disH)-1);zeros(1,length(disH)-1)];
dummy = [zeros(length(disV),1),dummy];
JV = sparse((nV-1)/(vGridU-vGridL)*(DisMat+dummy));
JH = sparse((nH-1)/(hGridU-hGridL)*(DisMat+dummy));
P0 = birth(:);
fV = alpha.*exp(betaV.*(Vgrid-Hgrid));
fH = alpha.*exp(betaH.*(Hgrid-Vgrid));
%Laplacian
dummy1 = [zeros(1,length(disH)-1);1*eye(length(disV)-1,length(disH)-1)];
dummy1 = [dummy1,zeros(length(disV),1)];
Lv = sparse(DV.*(((nV-1)/(vGridU-vGridL)).^2).*(2*DisMat + dummy1 + dummy));
Lh = sparse(DH.*(((nH-1)/(hGridU-hGridL)).^2).*(2*DisMat + dummy1 + dummy));
%Solution using inbuilt MATLAB function
sol = ode45(@(t,P) RBCfun(t,P,fV,fH,Lv,Lh,JV,JH,birth,death,nV,nH),[0 0.01],P0);
