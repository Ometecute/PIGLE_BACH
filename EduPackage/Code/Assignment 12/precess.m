function [Sf] = precess(Si, Z, V, B)
% precondition: Si = initial spin of all atoms at each of the na different speeds and different points P
%               Z = the set of points along the axis of the solenoid
%               B = matrix of the magnetic field at each point (P,Z) in the solenoid
%               V0 = average speed of incident atoms
%               DV = fractional spread in speeds of incident atoms
%               na = number of atoms considered
%
% postcondition: Sf = set of final spins for each speed and lateral point P
%
    
    %---------------------------------------------------------------
    % set up dynamic variables of system
    %---------------------------------------------------------------
    
    % number of Z-points within solenoid
    nZ = length(Z);
    nP = size(B,2);
    nV = length(V);

    % produce a linearly-varying set of speeds for each atom in the given range
    tmp(:,1,1) = V;
    % repeat this over the extent of the beam, i.e. constant speed since particles are uncharged
    V = repmat(tmp(:,1,1), [1 nP nZ]);

    % calculate the time taken for each atom to travel each considered distance => function t(d)
    t = zeros(size(V));
    t(:,1,:) = reshape(repmat(Z', [nV 1]), [nV 1 nZ])./V(:,1,:);
    % repeat this over the extent of the beam, i.e. constant speed since particles are uncharged
    t = repmat(t(:,1,:), [1 nP 1]);

    % differentiate t(Z) once wrt. Z to obtain the time interval dt taken to travel between each considered distance along the length
    dt=diff(t,1,3);
    dt(:,:,end+1)=dt(:,:,end);

    % normalise spin to get initial orientation (unit vector):
    % initial spins should all have the same magnitude...
    Simag = sqrt(sum(Si(:,1).^2, 1));
    
    % ...but possibly different directions for each speed
    Sidir=Si./(sqrt(sum(Si.^2))+eps);
    
    %---------------------------------------------------------------
    %---------------------------------------------------------------

    
    
    
    
    
    %---------------------------------------------------------------
    % calculate Larmor frequency
    %---------------------------------------------------------------

    % gyromagnetic ratio for helium 3
    gamma=2.0378e8; %rad/sec;

    % calculate Larmor frequency for atoms of all speeds b/c independent of speed
    w=gamma*sqrt(sum(B.^2,1));

    % magnitude of B
    Bmag=sqrt(sum((B.^2),1));

    % normalise B since we only needed the magnitude for the Larmor frequency (eps rounds to the nearest floating-point number)
    Bdir=B./(eps+repmat(Bmag, [3 1 1]));
    
    %---------------------------------------------------------------
    %---------------------------------------------------------------
    
    
    
    
    
    
    
    %---------------------------------------------------------------
    % precess spin vectors of atoms along the length of the field
    %---------------------------------------------------------------

    % rotation matrix for precessing spin at point in the solenoid
    R=zeros(3,3,nP,nZ);

    % final spin direction for each atom
    spinvec1=zeros(3,nP,nZ,nV);

    % matrix of *change* in precession angle in travelling between each considered distance, for each atom
    dphimat=repmat(w, [nV 1 1]).*dt;

    % for each velocity iV...
    for iV=1:nV
        % change in angle phi in each step along the length for atom in1
        dphi = (dphimat(iV,:,:));

        % calculate the rotation matrix for atoms of speed V(iV) that acts on the spin at each point for every step along the solenoid

        R(1,1,:,:)=(Bdir(1,:,:).^2).*(1-cos(dphi))+cos(dphi);
        R(1,2,:,:)=(Bdir(1,:,:).*Bdir(2,:,:)).*(1-cos(dphi))-Bdir(3,:,:).*sin(dphi);
        R(1,3,:,:)=(Bdir(1,:,:).*Bdir(3,:,:)).*(1-cos(dphi))+Bdir(2,:,:).*sin(dphi);


        R(2,1,:,:)=(Bdir(1,:,:).*Bdir(2,:,:)).*(1-cos(dphi))+Bdir(3,:,:).*sin(dphi);
        R(2,2,:,:)=(Bdir(2,:,:).^2).*(1-cos(dphi))+cos(dphi);
        R(2,3,:,:)=(Bdir(2,:,:).*Bdir(3,:,:)).*(1-cos(dphi))-Bdir(1,:,:).*sin(dphi);


        R(3,1,:,:)=(Bdir(1,:,:).*Bdir(3,:,:)).*(1-cos(dphi))-Bdir(2,:,:).*sin(dphi);
        R(3,2,:,:)=(Bdir(2,:,:).*Bdir(3,:,:)).*(1-cos(dphi))+Bdir(1,:,:).*sin(dphi);
        R(3,3,:,:)=(Bdir(3,:,:).^2).*(1-cos(dphi))+cos(dphi);

        % initial spin for velocity iV at Z-point Z(iZ=1)=0
        spinvec1(:,:,1,iV)=Sidir(:,:,iV);

        % for points iP at each further step Z(iZ) along the length, apply the corresponding rotation operator
        for iP = 1:nP
            for iZ=2:nZ
                spinvec1(:,iP,iZ,iV)=R(:,:,iP,iZ)*spinvec1(:,iP,iZ-1,iV);         
            end    
        end
    end
    
    % output spin, making sure dimensions are correct (so no squeezing of 1-dim axes)
    Sf = Simag .* reshape(spinvec1, [3 nP nZ nV]);
end