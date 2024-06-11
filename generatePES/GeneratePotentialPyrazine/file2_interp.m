clear;

%--------------------------------------------------------------------------
% load in Potential.mat
%--------------------------------------------------------------------------

load('Potential_AQ_base_file.mat');
% contains: - matrix of grid point (x,y)-positions
%           - vector of AQ angles
%           - vector of weights for each point
%           - matrix of sampled potentials (i.e. Fourier components at each point and each angle)

%--------------------------------------------------------------------------
% plot the locations of the 157 grid points
%--------------------------------------------------------------------------

% Plotting of 157 pts.
figure; hold on; axis equal;
x_157=xypos(:,1); y_157=xypos(:,2);
for i=1:size(xypos,1)
    plot (x_157(i),y_157(i),'.')
    text(x_157(i),y_157(i),num2str(i));
end



% Number of Interpolation Points for rectangular grid & angle interpolation
nx=120;
ny=200;
ntheta=2; % to request 30 theta vals., starting at 0, e.g., set ntheta=31 | 2-1=1 if no angular dependency is used

%% define range over which potential will be calculated
% potential is calculated across region (x,y) \in ([0,ax],[0,ay])

% the maximum x-displacement of the WS cell is half the basic horizontal lattice spacing
ax=2*max(xypos(:,1));

% the maximum y-displacement of the WS cell is half the basic vertical lattice spacing
ay=4*max(xypos(:,2));

% axes are typically chosen such that dir(x) || dir(a1)
a=ax;

% number of points in non-interpolated WS cell
np=19;

V3Dinterp=V3Dinterp_AQ;

%% --------------------------------------------------------------------------
% interpolating in theta
%--------------------------------------------------------------------------

% intitial number of angles
ntheta0=length(thetapos);

% generate new linearly-spaced vector from [0,2*pi) with ntheta-1 entries
thetavec=linspace(0,2*pi,ntheta); thetavec=thetavec(1:(ntheta-1));

ntheta=ntheta-1;

% use cubic interpolation with spline() function to interpolate in theta dimension for each position (x,y)
% for ixy=1:(size(xypos,1))
%     V3Dinterp_temp(ixy,:)=spline(thetapos,V3Dinterp(ixy,:),thetavec); % spline() uses cubic interpolation
% end
% V3Dinterp=V3Dinterp_temp;

%% --------------------------------------------------------------------------
% interpolating in x and y
%--------------------------------------------------------------------------

% original integer factor N as defined in the theoretical bckgd
nscale=2;

% use halvespacing() method to double the resolution of V3Dinterp
[V3Dnew, xyposnew, weightnew, npnew, nscalenew] = halvespacing(V3Dinterp, xypos, weight, np,nscale,a,12);

%[V3Dnew, xyposnew, weightnew, npnew, nscalenew] = halvespacing(V3Dnew, xyposnew, weightnew, npnew,nscalenew,a,16);


%% --------------------------------------------------------------------------
% interpolating in x and y
%--------------------------------------------

% basic vertical reciprocal lattice spacing
G1=4*pi/sqrt(3)/a;

% basic vertical reciprocal lattice spacing
Gx0=2*pi/ax;

% basic vertical reciprocal lattice spacing
Gy0=2*pi/ay;

%this gives twice the point spacing in the original grid

% iGp is the number of mesh cells in real space to lattice cell along one direction
iGp=0;

% iterate through desired no of points
for ix=1:nx
    for iy=1:ny
        if(mod(ix+iy,2)==0)
            
            % place origin at centre of grid by limiting range to +/- nx or ny:
            
            % limit x
            ixind=ix-1;
            if(ixind>nx/2)
                ixind=ixind-nx;
            end
            
            % limit y
            iyind=iy-1;
            if(iyind>ny/2)
               iyind=iyind-ny;
            end
            
            % find corresponding reciprocal vector for given posn, Gtemp
            Gtemp(1)=ixind*Gx0; Gtemp(2)=iyind*Gy0;
            
            % weight that position based on its position within the range
            weighttemp=InHexRecip(Gtemp,G1,nscalenew);
            
            % if the point is within the basic reciprocal WS cell... (weight=1,1/2,1/3)
           if(weighttemp~=0)
                % add the point to the set of considered reciprocal vectors
                iGp=iGp+1;
                
                Gxp(iGp)=ixind*Gx0; Gyp(iGp)=iyind*Gy0;
                
                weightp(iGp)=weighttemp;
            end
        end
    end
end

% total no of reciprocal grid vectors in basic reciprocal WS cell
nGp=iGp;

% set of moduli of reciprocal grid lattice vectors
Gmag=(Gxp.*Gxp+Gyp.*Gyp);

% sort Gmag and store the permutatiion of the old order to the new order in IX
[Gmag, IX]=sort(Gmag);

% we can now use IX to reorder the points (Gxp,Gyp) by their modulus (reorder all their properties)
Gdum=Gxp(IX); Gxp=Gdum;
Gdum=Gyp(IX); Gyp=Gdum;

Gdum=weightp(IX); weightp=Gdum;

% Normalise the (ixGp,iyGp) to integers defining the reciprocal vectors
ixGp=round(Gxp./Gx0);
iyGp=round(Gyp./Gy0);

% Switch to 1-indexing for use in the Fourier transform
ix0Gp=ixGp+1;
iy0Gp=iyGp+1;

% add multiples of nx,ny to map all points (ix0Gp,iy0Gp)<1, i.e. points outside the range [1,nx] & [1,ny] to 1,2,...
for iGp=1:nGp
    if(ix0Gp(iGp)<1)
        ix0Gp(iGp)=ix0Gp(iGp)+nx;
    end
    if(iy0Gp(iGp)<1)
        iy0Gp(iGp)=iy0Gp(iGp)+ny;
    end
end

% Gxp,Gyp are the x and y components of the nGp reciprocal vectors used to describe the potential. 

% ixGp and iyGp give these vectors as integer multiples of Gx0 and Gy0 which are the basic G vectors of the rectangular unit cell

% ix0Gp and iy0Gp are these integers, but mapped on to the range 1 to nx and 1 to ny where (1,1) is the (Gx,Gy)=(0,0) - i.e. these are the integers to use in the fourier transforms

% so now ix0Gp and iy0Gp have the integers that will locate each G vector in the matrix for a iff2'ed potential

%% --------------------------------------------------------------------------
% calculate potential  Fourier components for interpolated potential
%--------------------------------------------------------------------------
%VG(iGp,iz)
VGnew=zeros(nGp,ntheta);
for iGp=1:nGp
    for itheta=1:ntheta
        dum=0.0;
        for ip=1:npnew  % sum weighted contributions of each Fourier component
            dum=dum+V3Dnew(ip,itheta)*exp(-1i*(Gxp(iGp)*xyposnew(ip,1)+Gyp(iGp)*xyposnew(ip,2)))*weightnew(ip);
        end
        VGnew(iGp,itheta)=dum/nscalenew/nscalenew;
    end
end

%% --------------------------------------------------------------------------
% weight the potential
%--------------------------------------------------------------------------

pe3D=zeros(nx,ny,ntheta);

for itheta=1:ntheta

    clear A
    A=zeros(nx,ny);
    for iGp=1:nGp
        % weight the potential at each point to account for double/treble-counting
        A(ix0Gp(iGp),iy0Gp(iGp))=VGnew(iGp,itheta)*weightp(iGp);
    end
    % take the real (cosine) part so that there is a top site at the origin
    B=real(fft2(A));
    pe3D(:,:,itheta)=B(:,:);
end

%% --------------------------------------------------------------------------
% produce grid of positions we wish to sample in our plot, as defined by nx and ny
%--------------------------------------------------------------------------

% Note we calculate it across a rectangle cell rather than a WS cell for ease of use
X=zeros(nx,ny);
Y=X;
for ix=1:nx
    for iy=1:ny
        X(ix,iy)=ax*(ix-1)/nx;
        Y(ix,iy)=ay*(iy-1)/ny;
    end
end

%% --------------------------------------------------------------------------
% plot potentials for each calculated theta in [0,180] degrees
%--------------------------------------------------------------------------

% define plot parameters
color_scale_min=0; color_scale_max=1050;
asp = [1 1 1/color_scale_max];
Lx = [0 max(X,[],'all')]; Ly = [0 max(Y,[],'all')];

% plot surface/contour map for each calculated angle theta in [0,180] degrees
for i = 1:ntheta
    figure(100+i);
    %stest = surf(X,Y,pe3D(:,:,i)); set(stest,'LineStyle','none'); caxis([color_scale_min color_scale_max]); view([0 90]); daspect(asp); xlim(Lx); ylim(Ly);
    stest = contourf(X,Y,pe3D(:,:,i)); axis equal; title(sprintf('%d',round(thetavec(i)*(180/pi))));
    
    %exportgraphics(gcf,sprintf('deg%d.png',round(thetavec(i)*3*(180/pi))),'Resolution','300')
end

%% --------------------------------------------------------------------------
% prepare potential object for implementation in MD simulation
%--------------------------------------------------------------------------

% Including 3rd z dimension BUT only with length 1 because we defined a single surface z(x,y)
PES4D_AQ = zeros(nx,ny,ntheta-1,1);
    
PES4D_AQ(:,:,1,:)=pe3D;

%permute axes t: y,x,theta,z
PES4D_AQ=permute(PES4D_AQ,[2 1 3 4]);

save('loadPES', 'PES4D_AQ');
