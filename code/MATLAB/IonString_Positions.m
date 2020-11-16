function positions = IonString_Positions(Nions,freq,ionmass)
% calculates the positions of a linear N-ion string of ions in micrometers
% Nions: number of ions
% freq: axial frequency (Hz)
% ionmass: ion mass in atomic units

% constants
PhysC.el  =  1.60217653e-19;    % elementary charge
PhysC.eps0=  8.854187817e-12; % electric constant
PhysC.amu =  1.66053886e-27;    % atomic mass unit

% find equilibrium positions of the ions (in dimensionless units)
u = equilibriumpositions_fun(Nions);

 % length scale in which the equilibrium positions are given
lscale = (PhysC.el^2/(4*pi*PhysC.eps0*PhysC.amu*ionmass*(2*pi*freq)^2))^(1/3); 
positions = u*lscale*1e6; % positions in micrometers


%--------------------------------------------------------------------------
% 8.2.2011
% This function calculates the equilibrium positions of a linear ion string
% of N ions.
% It gives as an output the equilibrium positions in units of 
% l=((Z^2e^2)/(4*pi*epsilon_0*m*omega^2))^(1/3)
% cf. D. James, Appl. Phys. B 66, 181 (1998)


function y = equilibriumpositions_fun(N) 

options = optimset('Display','off');
%y0 = -(N-1)/2:(N-1)/2;
y0 = (-(N-1)/2:(N-1)/2)*2/sqrt(N); % this seems to be quite a reasonable starting point

y = fsolve(@forces_fun,y0,options);

%--------------------------------------

function y = forces_fun(x)

N = length(x);
U = ones(N,1)*x;

% calculate the differences between all ion positions
UDiff = U-U';

% Coulomb force acting on ion m caused by ion n
cforcematrix = (tril(ones(N))'-tril(ones(N)))./UDiff.^2;
cforcematrix(eye(N)==1)=0;
y = -x + sum(cforcematrix); % add trap confining force and coulomb force


