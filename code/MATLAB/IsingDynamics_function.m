function [lscale eigfreqx,eigfreqy,eigfreqz,ionpos] = IsingDynamics(Params)
%
% 15.10.2013, CR
% This program calculates the equilibrium positions and normal modes of a
% linear string of trapped ions (all having the same mass and charge).
% Based on this information, the spin-spin coupling constants induced by a
% bichromatic field off-resonantly coupling to the red and the blue sidebands 
% are calculated. In the last step, the program calculates the spin
% dynamics induced by such an Ising interaction in the presence of a
% transverse field

%clear all
global PhC

CALCULATESPINSPINCOUPLINGS = 0;   % set to '0' if Jij has already been calculated
SAVEALLDATA = 1;                  % set to '1' if you want to save all data

%--------------------------------------------
% Parameters 
%--------------------------------------------
%these experimental parameters are now loaded in when the function is
%called

%--------------------------
% Relevant physical constants
PhC = PhysConst_fun;


    %% --------------------------------------------
    % Step 1: calculate the equilibrium positions
    %%--------------------------------------------

        [ionpos,lscale] = CalculateEquilibriumPositions_fun(Params);

        if(max(max(abs(ionpos(:,1:2))))>1e-6)
            'THE ION CRYSTAL IS NOT LINEAR !!!'
            disp(['Maximum distance from center along x/y = ' sprintf('%e',lscale*max(max(abs(ionpos(:,1:2)))))])
            return  % stop further execution of the program
        end

%         % display ion positions
%         disp(sprintf('Ion positions in micrometers:'))
%         pos_string=[];
%         for j = Params.N:-1:1
%             pos_string = [sprintf('%2.1f ',ionpos(j,3)'*lscale*1e6) pos_string];
%         end
%         disp(pos_string)



    %% --------------------------------------------------------------------------
    % Step 2: Use ion positions to calculate normal modes
    %%--------------------------------------------------------------------------

        [eigfreqx,eigfreqy,eigfreqz,Vx,Vy,Vz] = CalculateNormalModeFreqs(Params,ionpos);


    %% --------------------------------------------------------------------------
    % Step 3: Calculate spin-spin coupling matrix
    %%--------------------------------------------------------------------------
% 
%         Jij_x = CalculateSpinSpinCouplings(Params,Vx,eigfreqx);
%         Jij_y = CalculateSpinSpinCouplings(Params,Vy,eigfreqy);
% 
%         Jij = Jij_x + Jij_y;
% 
%         for j=1:Params.N-1
%             Ji(j) = mean(diag(Jij,j));
%         end    


%% -------------------------------------------------------------------------
% SUBROUTINES
%%--------------------------------------------------------------------------
function [ionpos,lscale] = CalculateEquilibriumPositions_fun(Params)
%% 22.2.2010
% This program is used to calculate equilibrium positions of an ion crystal
% in an anisotropic harmonic potential
% For this, the trajectories of ions interacting with each other in this
% potential are calculated. A damping force is introduced to relax the ions
% to their equilibrium positions.

        global PhC
        SHOWTRAJECTORIES = 0;   % set to 1 to display the relaxation of the ions to their equilibrium positions
        
        % time interval for simulation
        tstart = 0;
        tstop = 200; 200;
        
        % damping constant
        gam = .3;
        
        % mass of ion in atomic mass units:
        m_amu = Params.m*PhC.amu;   

        % length scale in which the equilibrium positions are given
        lscale = (PhC.el^2/(4*pi*PhC.eps0*m_amu*Params.referencefreq^2))^(1/3); 

        % spring constants of the harmonic part of the external potential
        % (normalized to reference frequency
        Dx = (Params.omx/Params.referencefreq)^2*sign(Params.omx);  % the sign is also sensed to allow for deconfining harmonic potentials (if there is a confining quadrupole potential present)
        Dy = (Params.omy/Params.referencefreq)^2*sign(Params.omy);
        Dz = (Params.omz/Params.referencefreq)^2*sign(Params.omz);

        N = Params.N; % number of ions
        
        %--- start with random ion positions ---

        % ion positions and velocities: column vector
        % entries: 1 ... N    x-positions
        %          N+1 ...2N  y-positions
        %         2N+1 ...3N  z-positions
        %         3N+1 ...4N  x-velocities
        %         4N+1 ...5N  y-velocities
        %         5N+1 ...6N  z-velocities

        %--- choose starting positions:we try to make a good guess of the initial conditions
        % the approach is to start with an equidistant ion chain with an ion-ion distance minimizing the potential energy
        rrdot0 = zeros(N,6);     
        a0_2 = (log(N)/N^2*6/Dz)^(1/3);
        a0_4 = (log(N)/N^4*40/Params.gamma4)^(1/5);
        if(a0_4==Inf)
            aa = a0_2;
        elseif(a0_2==Inf)
            aa = a0_4;
        else
            aa = min([a0_2 a0_4]);
        end
        rrdot0(:,1) = rand(N,1)*aa/100;
        rrdot0(:,2) = rand(N,1)*aa/100;
        rrdot0(:,3) = [(1-N):2:(N-1)]*aa/2;
        rrdot0 = reshape(rrdot0,6*N,1);
        %---

        options = odeset('RelTol',1e-6); 
        options = odeset(options,'AbsTol',1e-8);

        %--- solve differential equation ---
        [t,rrdot] = ode113(@equilibriumpositions_x4_fun,[tstart tstop],rrdot0,options,[Dx Dy Dz Params.gamma4 gam]);
        %-----------------------------------

        tN = length(t);
        rrdotFinal = rrdot(tN,:);
        ionpos = reshape(rrdotFinal(1:3*N),N,3);

        %--- plot ion trajectories relaxing to their final positions? ---
        if(SHOWTRAJECTORIES)
            figure(111)
            clf
            subplot(3,1,1)
            hold on
            for j=1:N
               plot(t/2/pi,rrdot(:,j),'b')
            end
            hold off

            subplot(3,1,2)
            hold on
            for j=1:N
               plot(t/2/pi,rrdot(:,N+j),'b')
            end
            hold off

            subplot(3,1,3)
            hold on
            for j=1:N
               plot(t/2/pi,rrdot(:,2*N+j),'b')
            end
            hold off
        end
        
  

        
%%--------------------------------------------------------------------------
function timederivs = equilibriumpositions_x4_fun(t,xyzv,par)
%%
% Dynamical equations of harmonically trapped ions coupled by Coulomb force
% The goal is to calculate the equilibrium positions of the ion crystal
%
% par(1) : trap frequency omx in terms of the reference frequency
% par(2) : trap frequency omy in terms of the reference frequency
% par(3) : trap frequency omz in terms of the reference frequency
% par(4) : parameter specifying strength of quartic part along z-axis
% par(5) : damping constant

% ion positions and velocities: column vector xyzv
% entries: 1 ... N    x-positions
%          N+1 ...2N  y-positions
%         2N+1 ...3N  z-positions
%         3N+1 ...4N  x-velocities
%         4N+1 ...5N  y-velocities
%         5N+1 ...6N  z-velocities

    
    Dx = par(1);               % spring constant in x-direction
    Dy = par(2);               % spring constant in y-direction
    Dz = par(3);               % spring constant in z-direction 
    gamma4 = par(4);              % parameter quantifying the strength of the quartic part of the potential along the axis of the trap  
    gam = par(5);                 % damping constant 

    N = length(xyzv)/6;
    ind=1:N;

    %-------------------------
    % calculate velocities
    drdt = xyzv(3*N+1:6*N); 

    % calculate accelerations (more complicated)
    x = xyzv(ind);
    y = xyzv(ind+N);
    z = xyzv(ind+2*N);
    vx = xyzv(ind+3*N);
    vy = xyzv(ind+4*N);
    vz = xyzv(ind+5*N);

    % calculation of difference of coordinates
    xx = x*ones(1,N);
    xx = xx'-xx;

    yy = y*ones(1,N);
    yy = yy'-yy;

    zz = z*ones(1,N);
    zz = zz'-zz;

    rr = (xx.^2 + yy.^2 + zz.^2).^(-1.5);
    rr(rr==Inf)=0;

    xcoul=sum(xx.*rr)';  % Coulomb force for x-component
    ycoul=sum(yy.*rr)';  % Coulomb force for y-component
    zcoul=sum(zz.*rr)';  % Coulomb force for z-component

    % Coulomb force contributions:
    xyzcoul = [xcoul;ycoul;zcoul];

    % trap confinement force
    xyztrap = -[Dx*x;Dy*y;Dz*z]; % now the transverse restoring force depends on the ion mass
    xyztrap4= -[0*ones(size(x));0*ones(size(y));2*gamma4*z.^3];     % quartic part of the potential; FOR THE MOMENT, WE LOOK JUST AT ITS INFLUENCE ALONG THE AXIS OF THE TRAP!!!

    % damping force
    xyzdamp = -gam*[vx;vy;vz]; % damp all motion equally
    %xyzdamp = -gam*[vx;vy;vz].*[0*damping;0*damping;damping]; % only z-motion damped

    % everything added up:
    dvdt = (xyzcoul + xyztrap + xyztrap4 + xyzdamp);

    %--------------------------
    timederivs = [drdt; dvdt];


    
    
    
    
%%--------------------------------------------------------------------------
function [eigfreqx,eigfreqy,eigfreqz,Vx,Vy,Vz] = CalculateNormalModeFreqs(Params,ionpos)
%% Use the ion positions, trap frequencies and ion mass to calculate the
% normal modes of the ion string
%
%  4.5.2012, Christian
% 24.4.2012: Can we do a calculation where we calculate the normal modes of
% a (mixed) crystal which is transversally confined by an rf-field and
% longitudinally by a static field composed of a quadratic and a quartic
% part?
% WORK IN PROGRESS!!! 
% For the moment, I've added a z^4 term without caring about a
% violation of the Laplace equation. It seems to me that if we a assume a
% radial symmetry in the transverse plane, the potential would have to be
% ~ z^4-3z^2(x^2+y^2)+3/8(x^2+y^2)^2
%
% WORK IN PROGRESS!!!
% Also for the normal mode calculation, only the curvature along z due to
% the quartic term enters in to the calculation 
%
% 1.3.2012
% Calculation of the normal mode spectrum of a mixed crystal of
% arbitrary shape in a linear trap
% 1. step: for a given potential, calculate the equilibrium positions
% 2. step: use the calculated positions to reduce the Lagrangian to second
% order in ion excursions and diagonalize the resulting set of differential
% equations.
% The way the program is currently written assumes that omega_z~1/sqrt(m)
% and omega_x~1/m, omega_y~1/m where 'm' is the mass of the ion 
%


        N = Params.N;
        % relative frequencies (stored in 'relfreqs')
        freqs = zeros(N,3);
        freqs(:,1) = Params.omx;    %in the transverse direction, the frequency scales with 1/m (we assume a pure rf-potential)
        freqs(:,2) = Params.omy;    %in the transverse direction, the frequency scales with 1/m (we assume a pure rf-potential)
        freqs(:,3) = Params.omz;    %in the longitudinal direction, the frequency scales with 1/sqrt(m)  (we assume a pure static potential)
        relfreqs = freqs/Params.referencefreq;           % normalized frequencies

        % --- rename array of equilibrium positions ---
        r0 = ionpos;  % first index: ion number (1-N), second index: coordinate (1-3)

        % --- calculate the difference of ion positions ---
        r0diff = zeros(N,N,3);  % contains the differences of the three coordinates between all pairs of ions
        r0dist = zeros(N,N);    % contains the distances between all pairs of ions
        for j=1:N
            for k=1:N
                r0diff(j,k,:) = r0(j,:)-r0(k,:);
                r0dist(j,k) = norm(squeeze(r0diff(j,k,:)));
            end
        end

        %--- for convenience, create also N x N arrays containing r0dist^-5 without
        %the need of removing infinities on the diagonal
        r0dist_min5 = (1./(r0dist-eye(N))+eye(N)).^5;    % array |r_i - r_j|^(-5)
        r0dist_2 = r0dist.^2;                            % array |r_i - r_j|^2

        %--- now calculate the N x N x 3 x 3 array that looks like a dipole-dipole
        %interaction ---
        tildA = zeros(N,N,3,3);
        elem = zeros(3,3);
        for j=1:N
            for k=1:N
                for al=1:3
                    for be=1:3
                        elem(al,be) = 3*r0diff(j,k,al)*r0diff(j,k,be) - (al==be)*r0dist_2(j,k);
                    end
                end
                tildA(j,k,:,:) = (r0dist_min5(j,k)*ones(3)).*elem;
            end
        end


        %--- sum up elements of tildA ---
        tildAsummed = zeros(N,3,3);
        for j=1:N
            for al=1:3
                for be=1:3
                    tildAsummed(j,al,be) = sum(squeeze(tildA(j,:,al,be)));   % we do not need to explicitly take into account that we should not include the element tildA(j,j,al,be) in the sum as this element is zero anyway
                end
            end
        end


        %--- calculate the oscillation frequencies with respect to the reference oscillation frequency ---

        %--- calculate matrix to be diagonalized ---
        inds = zeros(N,3);
        for j=1:N
            for al=1:3;
                inds(j,al) = j+ (al-1)*N;    % this array help putting together the matrix to be diagonalized
            end
        end

        B = zeros(N*3,N*3);
        for j=1:N
            for k=1:N
                for al=1:3
                    for be=1:3
                        B(inds(j,al),inds(k,be)) =  (j==k)*(al==be)*relfreqs(j,al)^2 + ((j==k)*tildAsummed(j,al,be)-tildA(j,k,al,be));
                    end
                end
            end
        end

         %--- Here, we include the curvature along z and x/y due to the quartic
         %term ~(z^4-3z^2(x^2+y^2)+3/8(x^2+y^2)^2). This treatment is valid
         %only for a linear crystal along the z-axis of the trap !!!!!!!!!!!!!
         for j=1:N
             B(inds(j,1),inds(j,1)) = B(inds(j,1),inds(j,1)) - 3*Params.gamma4*r0(j,3)^2;
             B(inds(j,2),inds(j,2)) = B(inds(j,2),inds(j,2)) - 3*Params.gamma4*r0(j,3)^2;
             B(inds(j,3),inds(j,3)) = B(inds(j,3),inds(j,3)) + 6*Params.gamma4*r0(j,3)^2;
         end

        %--- Diagonalize the matrix B ---
        [V,D] = eig(B);
        eigfreqs = sqrt(diag(D));   % eigenfrequencies in units of the reference frequency

        %figure(3)
        %imagesc(abs(B))

        %'eigenfrequencies in MHz'
        %eigfreqs'*referencefreq/(2*pi*1e6)

    %     figure(4)
    %     bar(sort(real(eigfreqs))*referencefreq/(2*pi*1e6))
    %     title('distribution of eigenfrequencies')
    %     ylabel('frequency')    


        sumx = sum(abs(V(0*N+1:1*N,:)).^2);
        sumy = sum(abs(V(1*N+1:2*N,:)).^2);
        sumz = sum(abs(V(2*N+1:3*N,:)).^2);

        indx = find(sumx>.9);
        indy = find(sumy>.9);
        indz = find(sumz>.9);

        eigfreqx = sqrt(diag(D(indx,indx)))';
        eigfreqy = sqrt(diag(D(indy,indy)))';
        eigfreqz = sqrt(diag(D(indz,indz)))';

        Vx=V(0*N+1:1*N,indx);
        Vy=V(1*N+1:2*N,indy);
        Vz=V(2*N+1:3*N,indz);

        % sort eigenfrequencies 
        [eigfreqx,ind_eigx] = sort(eigfreqx);
        [eigfreqy,ind_eigy] = sort(eigfreqy);
        [eigfreqz,ind_eigz] = sort(eigfreqz);

        % sort positions of the ions along z
        [z0,ind_z0] = sort(r0(:,3));

        % sort eigenvectors and their coordinates according to new order or
        % eigenvalues and positions
        Vx = Vx(ind_z0,ind_eigx);
        Vy = Vy(ind_z0,ind_eigy);
        Vz = Vz(ind_z0,ind_eigz);
        
%%--------------------------------------------------------------------------
function PhC = PhysConst_fun;
%% some physical constants that we need

    global PhC
      
    % Physical constants: CODATA 2002 values
    PhC.eps0 = 8.854187817e-12; % electric constant
    PhC.el = 1.60217653e-19;    % elementary charge
    PhC.amu = 1.66053886e-27;    % atomic mass unit
    PhC.hbar = 1.05457168e-34;  % hbar 

    
%%--------------------------------------------------------------------------
function Jij = CalculateSpinSpinCouplings(Params,V,eigfreq)
%% calculate spin spin coupling matrix    

    global PhC
    
    % Carrier Rabi frequency
    omega=pi./(Params.tpi); %going for the 2 x pi-pulse definition of rabi frequency.
    laser_det =  max([Params.omx Params.omy])/Params.referencefreq + Params.laser_off; % detuning of the laser from the carrier

    N = Params.N;
    Jij = zeros(N,N);

    ionmass = Params.m*PhC.amu;
    kvec = 2*pi/Params.lam;
    if(abs(max(eigfreq)-max([Params.omx Params.omy]/Params.referencefreq))<1e-6)
        mode_overlap_squared = Params.geo^2;  % cosine squared of the angle between the direction of highest confinement and the laser beam
    else
        mode_overlap_squared = 1-Params.geo^2;
    end
    prefactor = (omega'*omega) * mode_overlap_squared*PhC.hbar*kvec^2/(4*ionmass);
    for j=1:N
        Jij_mode = prefactor .* (V(:,j)*V(:,j)')/((laser_det^2-eigfreq(j)^2)*Params.referencefreq^2);
        Jij = Jij + Jij_mode;  ; %in units angular freq (s^-1}
    end
    Jij = Jij.*(1-eye(N)); % get rid of diagonal elements
    xxx=1;


        
        

