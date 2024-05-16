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
nV=4;

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

% width/height of solenoid /m
Lx = 0.05;
Ly = 0.05;

%---------------------------------------------------------------
%---------------------------------------------------------------



%---------------------------------------------------------------
% solenoid dimensions of interest
%---------------------------------------------------------------

% SETUP:

% length of the region of interest (symmetrical about solenoid)
Zmax = 5/4*L;

% number of Z-points to simulate
nZ = 750;


% maximum allowed beam radius
a = min(Lx, Ly)/2;

% number of radii to simulate
nr = 5;


% number of angles to simulate -- in range [0,2*pi)
ntheta = 6;


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

[Bx,By] = s.box_xy();

Bz = s.box_z();
% the smoothing length is chosen to maximise smoothing over the considered region Z
%Bz = s.smoothstep_z((Z(end) - L)/2);

% the standard deviation is chosen to be a third of the radius of the solenoid so that the field strength is close to zero at the edges
%radial_decay_const = 3 * 1/(2*(min(Lx,Ly)/2)^2);
%for ij = 1:nP
%    for k = 1:nZ
%        Bz(ij,k) = Bz(ij,k) * exp(-(X(ij)^2)*radial_decay_const - (Y(ij)^2)*radial_decay_const);
%    end
%end

% intialise set of vectors for B-field at each point along the length
B=zeros([3 nP nZ]);
B(1,:,:) = Bx; B(2,:,:) = By; B(3,:,:) = Bz;

B = 0.15 .* B;

%---------------------------------------------------------------
%---------------------------------------------------------------




%---------------------------------------------------------------
% precess spins through B-field
%---------------------------------------------------------------

% precess spin through solenoid
Si = repmat(spin, [1 nP nV]);
Sf1 = precess(Si, Z, V, B);

%---------------------------------------------------------------
%---------------------------------------------------------------



%---------------------------------------------------------------
% plot spin components against z-distance for diff speeds
%---------------------------------------------------------------
figure();
iP = 1;

% plot x-component of spin against distance through solenoid
q(1) = subplot(6,1,3); 
plot(Z, squeeze(Sf1(1,iP,:,:)));
xlim([Z(1) Z(end)]); ylim([-Smag Smag]);
title('$S_x$');
set(gca,'xtick',[]); set(gca,'xticklabel',[]); yticks([-1 1])

% plot y-component of spin against distance through solenoid
q(2) = subplot(6,1,4);
plot(Z, squeeze(Sf1(2,iP,:,:)));
xlim([Z(1) Z(end)]); ylim([-Smag Smag]);
title('$S_y$');
set(gca,'xtick',[]); set(gca,'xticklabel',[]); yticks([-1 1])

% plot z-component of spin against distance through solenoid
q(3) = subplot(6,1,5);
plot(Z, squeeze(Sf1(3,iP,:,:)));
xlim([Z(1) Z(end)]); ylim([-Smag Smag]);
title('$S_z$');
set(gca,'xtick',[]); set(gca,'xticklabel',[]); yticks([-1 1])

% plot z-component of spin against distance through solenoid
q(4) = subplot(6,1,6);
%plot(Z, squeeze(B(3,iP,:))', 'color', '#483493');
%xlim([Z(1) Z(end)]); ylim([0 0.15]);
[x,z] = meshgrid(Z, X);
s = surf(x, z, squeeze(B(3,:,:))); view([0 90]); s.EdgeColor = 'none';
title('$B_z$'); xlim([Z(1) Z(end)]); ylim([min(X) max(X)]);
set(gca,'xtick',[]); set(gca,'xticklabel',[]); yticks([0 0.15])

%---------------------------------------------------------------
%---------------------------------------------------------------



%---------------------------------------------------------------
% plot spin components against z-distance for diff trajectories
%---------------------------------------------------------------
figure();
iV = 1;

% plot x-component of spin against distance through solenoid
p(1) = subplot(6,1,3); 
plot(Z, squeeze(Sf1(1,:,:,iV)));
title('$S_x$'); xlim([Z(1) Z(end)]); ylim([-Smag Smag]);
set(gca,'xtick',[]); set(gca,'xticklabel',[]); yticks([-1 1])

% plot y-component of spin against distance through solenoid
p(2) = subplot(6,1,4);
plot(Z, squeeze(Sf1(2,:,:,iV)));
title('$S_y$'); xlim([Z(1) Z(end)]); ylim([-Smag Smag]);
set(gca,'xtick',[]); set(gca,'xticklabel',[]); yticks([-1 1])

% plot z-component of spin against distance through solenoid
p(3) = subplot(6,1,5);
plot(Z, squeeze(Sf1(3,:,:,iV)));
title('$S_x$'); xlim([Z(1) Z(end)]); ylim([-Smag Smag]);
set(gca,'xtick',[]); set(gca,'xticklabel',[]); yticks([-1 1])

% plot z-component of spin against distance through solenoid as colormap
p(4) = subplot(6,1,6);
[x,z] = meshgrid(Z, X);
s = surf(x, z, squeeze(B(3,:,:))); view([0 90]); s.EdgeColor = 'none';
title('$B_z$'); xlim([Z(1) Z(end)]); ylim([min(X) max(X)]);
set(gca,'xtick',[]); set(gca,'xticklabel',[]); yticks([0 0.15])


%---------------------------------------------------------------
%---------------------------------------------------------------


pPos = get(p, 'Position');
qPos = get(p, 'Position');
for i = [1,2,3,4]
    % define new posn
    newpPos = pPos{i}(1:2);
    newqPos = qPos{i}(1:2);
    
    % modify y-coord
    newpPos(2) = 1.5 * newpPos(2);
    newqPos(2) = 1.5 * newqPos(2);
    
    % set new posn
    set(p(i), 'Position', [newpPos, pPos{i}(3:end)]);
    set(q(i), 'Position', [newqPos, qPos{i}(3:end)]);
end

%exportgraphics(gcf,'boxUniformPrec.pdf','ContentType','vector')

