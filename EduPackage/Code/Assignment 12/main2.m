clear;
%---------------------------------------------------------------
%---------------------------------------------------------------
% SETUP
%---------------------------------------------------------------
%---------------------------------------------------------------


%---------------------------------------------------------------
% initial atom parameters
%---------------------------------------------------------------

% average speed of He-3 atom through solenoid (along x-axis) in m/s
V0=768;

% maximum fractional difference of atom speeds from average
DV=0.1;

% number of velocities simulated at each point
nV=20;

% initialise vector of velocities
V = linspace(V0-(V0*DV),V0+(V0*DV), nV);

% initial spin polarisation (parallel to the y-axis)
spin=[0 1 0]'; Smag = sum(spin.^2);

% scattering angle at surface
sctAngle = pi/4;

%---------------------------------------------------------------
%---------------------------------------------------------------






%---------------------------------------------------------------
% solenoid parameters
%---------------------------------------------------------------

% length of solenoid /m
L = 0.005;

% width/height of rectangular solenoid /m
Lx = 0.05;
Ly = 0.05;

%---------------------------------------------------------------
%---------------------------------------------------------------



%---------------------------------------------------------------
% solenoid dimensions of interest
%---------------------------------------------------------------

% SETUP:

% length of the region of interest (symmetrical about solenoid)
Zmax = 2*L;

% number of Z-points to simulate
nZ = 750;


% maximum allowed beam radius
a = min(Lx, Ly)/2;

% number of radii to simulate
nr = 5;


% number of angles to simulate -- in range [0,2*pi)
ntheta = 8;


% CALCULATION:

% initialise vector of Z-points
Z = linspace(0,Zmax,nZ)';

% calculate spacing between Z-points (we wish to keep this consistent throughout the simulation
dZ = Z(2)-Z(1);


% set of radii around central beam in range [0,a)
r = linspace(0,a,nr+1); r(end) = [];


% set of angles around central beam in range [0,2*pi)
theta = linspace(0,2*pi,ntheta+1); theta(end) = [];



% meshgrid of all lateral points on the beamline/cylinder of atoms (the meshgrid here isn't actually important, just each meshgrid index corresponds to the field value in B at the same index (see below)
[R, Theta] = meshgrid(r, theta);

nP = numel(R); % number of points in meshgrid



% extract Cartesian grid for use with Solenoid class
X = reshape(R.*cos(Theta), [nP 1]); Y = reshape(R.*sin(Theta), [nP 1]);

%---------------------------------------------------------------
%---------------------------------------------------------------



%---------------------------------------------------------------
% solenoid B-field
%---------------------------------------------------------------

s = Solenoid(L, Lx, Ly, X', Y', Z'-Zmax/2);

smoothLen = (Z(end) - L)/2;

[Bx,By] = s.box_xy();

%Bz = s.box_z();

% the smoothing length is chosen to maximise smoothing over the considered region Z
Bz = s.smoothstep_z((Z(end) - L)/2);

% the standard deviation is chosen to be a third of the radius of the solenoid so that the field strength is close to zero at the edges
radial_decay_const = 3 * 1/(2*(min(Lx,Ly)/2)^2);
for ij = 1:nP
    for k = 1:nZ
        Bz(ij,k) = Bz(ij,k) * exp(-(X(ij)^2)*radial_decay_const - (Y(ij)^2)*radial_decay_const);
    end
end

% intialise set of vectors for B-field at each point along the length
B=zeros([3 nP nZ]);
B(1,:,:) = Bx; B(2,:,:) = By; B(3,:,:) = Bz;

B = 0.15 .* B;

%---------------------------------------------------------------
%---------------------------------------------------------------







%---------------------------------------------------------------
% spin rotator field parameters
%---------------------------------------------------------------

% this variable is: - 0 if there is no spin rotator field
%                   - 1 if the spin rotation is ideal (i.e. uniform rotation for all speeds)
%                   - 2 if the spin rotator field is (fast-transitioning) uniform (i.e. rotation depends on speed)
spinRotator = 2;

% distance from end of each solenoid to the spin rotator field
rotBDist= 0.05;

% length of the spin rotator field
rotBLen = 0.05;

%---------------------------------------------------------------
%---------------------------------------------------------------







%---------------------------------------------------------------
%---------------------------------------------------------------
% ACTUAL EXECUTION BEGINS
%---------------------------------------------------------------
%---------------------------------------------------------------

%---------------------------------------------------------------
% first solenoid
%---------------------------------------------------------------

% propagate short distance up to solenoid
Spins = repmat(spin, [1 nP 5 nV]);
Disps = [-0.01:0.002:-0.002]';

% precess spin through first solenoid
Si = repmat(spin, [1 nP nV]);
Sf1 = precess(Si, Z, V, B);

% update stored spins/distances
Spins = cat(3, Spins, Sf1);
Disps = cat(1, Disps, Z);

% continue propagating up to SR field
numStepsTorotB = rotBDist/dZ - 1;
Disps = cat(1, Disps, linspace(Disps(end)+dZ, Disps(end)+rotBDist, numStepsTorotB)');
Spins = cat(3, Spins, repmat(Spins(:,:,end,:), [1 1 numStepsTorotB 1]));

%---------------------------------------------------------------
%---------------------------------------------------------------






%---------------------------------------------------------------
% spin rotator field
%---------------------------------------------------------------
% Rotates spin about the x-axis through sctAngle rad without affecting trajectory (and same for ALL speeds
%  - The field is uniform along x, but the spread in angle due to different speeds is negligible, so is neglected here
%  - Since we take a projection of the final spin in measurement, the effect is further reduced

Nsr = rotBLen/dZ - 1;
Zsr  = linspace(Disps(end)+dZ, Disps(end)+rotBLen, Nsr)';

if spinRotator == 1
    % for an ideal uniform rotation for all speeds...
    
    % define small rotation for each step along axis
    RotMsr = [1 0 0; 0 cos(sctAngle/Nsr) -sin(sctAngle/Nsr); 0 sin(sctAngle/Nsr) cos(sctAngle/Nsr)];
    
    % initialise temporary spin variable
    Sfsr = zeros([3 nP Nsr nV]);
    
    % first rotation
    for iP = 1:nP
        Sfsr(:,iP,1,:) = RotMsr * squeeze(Spins(:,iP,end,:));
    end
    
    % iteratively rotate through small angle until sctAngle rad in total
    for iV=2:Nsr
        for iP = 1:nP
            Sfsr(:,iP,iV,:) = RotMsr * squeeze(Sfsr(:,iP,iV-1,:));
        end
    end

    % update stored spins
    Spins = cat(3, Spins, Sfsr);
    
elseif spinRotator == 2
    % for a (fast-transitioning) uniform SR field with variance in rotation for diff speeds...
    
    % CREATE SR FIELD:
    % intialise set of vectors for B-field at each point along the length
    Bsr = zeros([3 nP length(Zsr)]);    

    % for each point along the length...
    Bsr(1,:,:) = ones([1 nP, length(Zsr)]);
    
    % gyromagnetic ratio for helium 3
    gamma=2.0378e8; %rad/sec;
    
    % scale B-field for phi=sctAngle/Nsr -> <phi>=l*gamma*B/<v>
    Bsr = sctAngle*V0/(rotBLen*gamma) .* Bsr;
    
    % PRECESS SPIN THROUGH FIELD
    Sfsr = precess(Spins(:,:,end,:), Zsr, V, Bsr);

    % update stored spins
    Spins = cat(3, Spins, Sfsr);
    
else
    % for no SR field at all...
    
    % update stored spins
    Spins = cat(3, Spins, repmat(Spins(:,:,end,:), [1 1 Nsr 1]));
end

% update stored distances
Disps = cat(1, Disps, Zsr);

%---------------------------------------------------------------
%---------------------------------------------------------------







%---------------------------------------------------------------
% resolve Sf into new direction ready for second solenoid
%---------------------------------------------------------------

% rotation matrix between two bases
RotM = [1 0 0; 0 cos(sctAngle) sin(sctAngle); 0 -sin(sctAngle) cos(sctAngle)];

% rotate frame for final spins
Sfrot = zeros(size(Spins(:,:,end,:)));
for iP = 1:nP
    Sfrot(:,iP,:) = RotM * squeeze(Spins(:,iP,end,:));
end
 
% continue propagating up to next solenoid
Disps = cat(1, Disps, linspace(Disps(end)+dZ, Disps(end)+rotBDist, numStepsTorotB)');
Spins = cat(3, Spins, repmat(Spins(:,:,end,:), [1 1 numStepsTorotB 1]));

%---------------------------------------------------------------
%---------------------------------------------------------------



%---------------------------------------------------------------
% second solenoid
%---------------------------------------------------------------

% precess final spins through second solenoid (opposite field direction)
Sf2rot = precess(Sfrot, Z, V, -B);

%---------------------------------------------------------------
%---------------------------------------------------------------



%---------------------------------------------------------------
% resolve Sf2 back
%---------------------------------------------------------------

% rotate frame back to original
Sf2 = zeros(size(Sf2rot));
for iP=1:nP
    for iV=1:nV
        Sf2(:,iP,:,iV) = RotM \ squeeze(Sf2rot(:,iP,:,iV));
    end
end

% update stored spins/distances
Spins = cat(3, Spins, Sf2);
Disps = cat(1,Disps,Disps(end)+dZ+Z);

% propagate short distance after solenoid
Spins = cat(3, Spins, repmat(Spins(:,:,end,:), [1 1 5 1]));
Disps = cat(1, Disps, [Disps(end)+0.002:0.002:Disps(end)+0.01]');

%---------------------------------------------------------------
%---------------------------------------------------------------






%---------------------------------------------------------------
% plot graphs of spin for diff speeds on central axis
%---------------------------------------------------------------

figure();

% lateral point we wish to plot components for (central)
iP = 1;

% plot x-component of spin against distance through solenoid
subplot(3,1,1); 
plot(Disps, squeeze(Spins(1,iP,:,:)));
title('$x$ component'); xlim([Disps(1) Disps(end)]); ylim([-Smag Smag]);

% plot y-component of spin against distance through solenoid
subplot(3,1,2);
plot(Disps, squeeze(Spins(2,iP,:,:)));
title('$y$ component'); xlim([Disps(1) Disps(end)]); ylim([-Smag Smag]);

% plot z-component of spin against distance through solenoid
subplot(3,1,3);
plot(Disps, squeeze(Spins(3,iP,:,:)));
title('$z$ component'); xlim([Disps(1) Disps(end)]); ylim([-Smag Smag]);

%---------------------------------------------------------------
%---------------------------------------------------------------



%---------------------------------------------------------------
% plot graphs of spin for diff axes at a single speed
%---------------------------------------------------------------

figure();

% lateral point we wish to plot components for (central)
iV = 1;

% plot x-component of spin against distance through solenoid
subplot(3,1,1); 
plot(Disps, squeeze(Spins(1,:,:,iV)));
title('$x$ component'); xlim([Disps(1) Disps(end)]); ylim([-Smag Smag]);

% plot y-component of spin against distance through solenoid
subplot(3,1,2);
plot(Disps, squeeze(Spins(2,:,:,iV)));
title('$y$ component'); xlim([Disps(1) Disps(end)]); ylim([-Smag Smag]);

% plot z-component of spin against distance through solenoid
subplot(3,1,3);
plot(Disps, squeeze(Spins(3,:,:,iV)));
title('$z$ component'); xlim([Disps(1) Disps(end)]); ylim([-Smag Smag]);

%---------------------------------------------------------------
%---------------------------------------------------------------



