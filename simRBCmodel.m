clear
clc
close all
load ('RETdist.mat')
%%
x = [0 25 16 2.3 0.025 72]; %set containing [alpha, betaV, betaH, DV, DH, Vc];
% x = [0.05 26 16 0.014 0.0014 0.8]; %set containing [alpha, betaV, betaH, DV, DH, Vc];
alpha = x(1);
betaV = x(2);
betaH = x(3);
DV = x(4);
DH = x(5);
Vc = x(6);
%%
MCV = 100; %fl
MCH = 35; %g/dl
vGridU = 140/MCV;
hGridU = 45/MCH;
vGridL = 5/MCV;
hGridL = 0.1/MCH;
nV = 50;
nH = 50;

disV = linspace(vGridL,vGridU,nV); %1-D grid in Vol
disH = linspace(hGridL,hGridU,nH); %1-D grid in Hgb
[Vgrid,Hgrid]=meshgrid(disV,disH); %Creating a 2-D grid for V and H
%%
P0(:,1) = P0(:,1);
P0(:,2) = P0(:,2); 
P0fit = hist3(P0,{disV.*MCV,disH.*MCH});
% Birth
birth = P0fit;%mvnpdf([Vgrid(:),Hgrid(:)],[1.3,1.2],[30,10]); %this should be the RET distribution
birth=reshape(birth,length(disV),length(disH));
% Death
VVbar = MCV.*Vgrid;
HHbar = MCH.*Hgrid;
theta = atan(MCH/MCV)*ones(length(disV),length(disH))-atan(HHbar./VVbar);
Delta = 100 * ( cos(theta) .* sqrt(VVbar.^2+HHbar.^2)-Vc.* sqrt(MCH^2.*ones(length(disV),length(disH))+...
    MCV^2.*ones(length(disV),length(disH))) )./ (sqrt(MCH^2.*ones(length(disV),length(disH))+...
    MCV^2.*ones(length(disV),length(disH))));
death = 1./(1+exp(Delta)); %zeros(nV,nH);
% death = -(sign(Delta)-1)/2;
% Jacobian
DisMat = -1*eye(length(disV),length(disH));
dummy = [1*eye(length(disV)-1,length(disH)-1);zeros(1,length(disH)-1)];
dummy = [zeros(length(disV),1),dummy];
J = (DisMat+dummy);%sparse((nV-1)/(vGridU-vGridL)*(DisMat+dummy));
% JH = ((nH-1)/(hGridU-hGridL)*(DisMat+dummy));%sparse((nH-1)/(hGridU-hGridL)*(DisMat+dummy));
Po = birth(:);
fV = alpha.*MCV.*exp(betaV.*(Vgrid-Hgrid));
fH = alpha.*MCH.*exp(betaH.*(Hgrid-Vgrid));
%Laplacian
dummy1 = [zeros(1,length(disH)-1);1*eye(length(disV)-1,length(disH)-1)];
dummy1 = [dummy1,zeros(length(disV),1)];
% Lv = (DV.*(((nV-1)/(vGridU-vGridL)).^2).*(2*DisMat + dummy1 + dummy));
% L = (2*DisMat + dummy1 + dummy);
I = eye(nV,nH);
E = sparse(2:nV,1:nH-1,1,nV,nH);
D = E+E'-2*I;
L1 = kron(D.*DH.*((nH-1)/(hGridU-hGridL)),I)+kron(I,D.*DV.*((nV-1)/(vGridU-vGridL)));
J1 = kron(J.*(nV-1)./(vGridU-vGridL).*fV,I)+kron(I,J.*(nH-1)./(hGridU-hGridL).*fH);

%Solution using inbuilt MATLAB function
sol = ode45(@(t,P) RBCfun(t,P,fV,fH,L1,J1,birth,death,nV,nH,vGridU,vGridL,hGridL,hGridU,DH,DV),[0 1],Po);
figure
for i = 1:numel(sol.y(1,:))
  P = reshape(sol.y(:,i),nV,nH);
   surf(disV,disH,P)
%    axis tight manual
 axis([vGridL vGridU hGridL hGridU 0 max(max(sol.y))]);
  F(i) = getframe;
%   keyboard
end
% movie(F)