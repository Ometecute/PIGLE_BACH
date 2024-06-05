function [ V3Dnew, xyposnew, weightnew, npnew, nscalenew ] = halvespacing( V3D, xypos, weight, np,nscale,a,nneb)
% For a hexagonal cell it makes sense to sample it with a grid of points 
% of the same symmetry as the basic lattice - and this will be computed
% originally on a lattice with spacing a/nscale if the original lattice
% constant is a. This rouine interpolates between this original grid,
% doubling the point spacing, using a cubic interpolation scheme
% nneb is the number of neighbours to use in the intepolation
%find list of new points needed
%clear xyposreq weightreq
dum10=a;
dumvec0=size(V3D);
nz=dumvec0(2);
del=0.0001*a*a;
a2=a*a/nscale/nscale;
annn2=a2*3;
npnew=0;
radius=a*(1/sqrt(3.0)+2/nscale);
V3Dnew=V3D;
V3Dextra=V3D;
for ip1=1:np-1
    for ip2=ip1+1:np
        dum=(xypos(ip1,1)-xypos(ip2,1))^2+(xypos(ip1,2)-xypos(ip2,2))^2;
        
        if(abs(dum-annn2)<del) && (weight(ip1)+weight(ip2)<1.1)
            % we have a new edge point
            npnew=npnew+1;
            xyposreq(npnew,1)=0.5*(xypos(ip1,1)+xypos(ip2,1));
            xyposreq(npnew,2)=0.5*(xypos(ip1,2)+xypos(ip2,2));
            weightreq(npnew)=0.5;
        elseif( abs(dum-annn2)<del) 
            npnew=npnew+1;
            xyposreq(npnew,1)=0.5*(xypos(ip1,1)+xypos(ip2,1));
            xyposreq(npnew,2)=0.5*(xypos(ip1,2)+xypos(ip2,2));
            weightreq(npnew)=1;
        end
    end
end
xyposnew=xypos;
weightnew=weight;
xyposnew(np+1:np+npnew,:)=xyposreq(:,:);
weightnew(np+1:np+npnew)=weightreq(:,:);
nscalenew=nscale*2;
figure(2000)
clf
hold on; axis equal
plot(xypos(:,1),xypos(:,2),'rx')
plot(xyposreq(:,1),xyposreq(:,2),'go')
sumofrealspaceweightinhalfspacingcell=sum(weight)+sum(weightreq)
title(['orig. pts redx.  new pts green o, interpolation grid blue square'])

% this total should be equal to nscalenew^2
%
% we need points outside the original unit cell - so need to add points
% from copies of that cell along the 6 sides of the original ones
%
% first translate to left and right - one needs all the points except the
% corner points along y=0
xyposextra=xypos;
npextra=np;
del=a/nscale*0.01;
for ip=1:np
   if(abs(xypos(ip,1)+a/2)>del)
        % the point is not on the left row of the original cell so copy
        % right
        npextra=npextra+1;
        xyposextra(npextra,1)=xypos(ip,1)+a;
        xyposextra(npextra,2)=xypos(ip,2);
        if(xyposextra(npextra,1)^2+xyposextra(npextra,2)^2 >radius^2)
            npextra=npextra-1;
        else
            V3Dextra(npextra,:)=V3D(ip,:);
        end
    end
    if(abs(xypos(ip,1)-a/2)>del)
        % the point is not on the right row of the original cell so copy
        % left
        npextra=npextra+1;
        xyposextra(npextra,1)=xypos(ip,1)-a;
        xyposextra(npextra,2)=xypos(ip,2);
        if(xyposextra(npextra,1)^2+xyposextra(npextra,2)^2 >radius^2)
            npextra=npextra-1;
            else
            V3Dextra(npextra,:)=V3D(ip,:);
        end
    end
    if(weight(ip)==1)
        % the point is not on an edge or corner, so copy other 4
        npextra=npextra+1;
        xyposextra(npextra,1)=xypos(ip,1)+a/2;
        xyposextra(npextra,2)=xypos(ip,2)+a*sqrt(3)/2;
                if(xyposextra(npextra,1)^2+xyposextra(npextra,2)^2 >radius^2)
            npextra=npextra-1;
            else
            V3Dextra(npextra,:)=V3D(ip,:);
                end

        npextra=npextra+1;
        xyposextra(npextra,1)=xypos(ip,1)+a/2;
        xyposextra(npextra,2)=xypos(ip,2)-a*sqrt(3)/2;
                if(xyposextra(npextra,1)^2+xyposextra(npextra,2)^2 >radius^2)
            npextra=npextra-1;
            else
            V3Dextra(npextra,:)=V3D(ip,:);
                end

        npextra=npextra+1;
        xyposextra(npextra,1)=xypos(ip,1)-a/2;
        xyposextra(npextra,2)=xypos(ip,2)+a*sqrt(3)/2;
                if(xyposextra(npextra,1)^2+xyposextra(npextra,2)^2 >radius^2)
            npextra=npextra-1;
            else
            V3Dextra(npextra,:)=V3D(ip,:);
                end


        npextra=npextra+1;
        xyposextra(npextra,1)=xypos(ip,1)-a/2;
        xyposextra(npextra,2)=xypos(ip,2)-a*sqrt(3)/2;
                if(xyposextra(npextra,1)^2+xyposextra(npextra,2)^2 >radius^2)
            npextra=npextra-1;
            else
            V3Dextra(npextra,:)=V3D(ip,:);
                end

    end
    if((abs(xypos(ip,1)-a/2)<del) && weight(ip)>0.4)
        % the point is  on the right row of the original cell so copy
        % and is not in a corner - so copy up and down diagonally 
        npextra=npextra+1;
        xyposextra(npextra,1)=xypos(ip,1)-a/2;
        xyposextra(npextra,2)=xypos(ip,2)+a*sqrt(3)/2;
                if(xyposextra(npextra,1)^2+xyposextra(npextra,2)^2 >radius^2)
            npextra=npextra-1;
            else
            V3Dextra(npextra,:)=V3D(ip,:);
                end
                
               npextra=npextra+1;
        xyposextra(npextra,1)=xypos(ip,1)-a/2;
        xyposextra(npextra,2)=xypos(ip,2)-a*sqrt(3)/2;
                if(xyposextra(npextra,1)^2+xyposextra(npextra,2)^2 >radius^2)
            npextra=npextra-1;
            else
            V3Dextra(npextra,:)=V3D(ip,:);
                end

    end
end
xyposextrasize=size(xyposextra);
V3Dextrasize=size(V3Dextra);
plot(xyposextra(:,1),xyposextra(:,2),'bs')
% we now have a candidate neighbours - npextra of them, xyposextra in
% position with values of potential stored in V3Dextra and we are looking
% to find potential values at the npnew points listed in xyposreq. V3Dnew
% already has the existing np points loaded into it.
%
% now find nearest neighbours of all required points
distance2=zeros(npnew,npextra);
distance=zeros(npnew,npextra,2);
dumvec=linspace(0,0,npextra);
ipextraindex=linspace(0,0,npextra);
distance0=linspace(0,0,nneb);
distancetemp=linspace(0,0,nneb);
D=linspace(0,0,nneb);
sigma=linspace(1,1,nneb);
sigma(5:nneb)=100;
tol=a/nscale*0.01;
iflag=0;
ni=10;
f=zeros(ni,nneb);
for ipnew=1:npnew
    for ipextra=1:npextra
        distance(ipnew,ipextra,:)=xyposextra(ipextra,:)-xyposreq(ipnew,:);
        distance2(ipnew,ipextra)=distance(ipnew,ipextra,1)^2+distance(ipnew,ipextra,2)^2;            
    end
    [dumvec,ipextraindex]=sort(squeeze(distance2(ipnew,:)));
    distance2(ipnew,:)=dumvec;
    distance(ipnew,:,1)=squeeze(distance(ipnew,ipextraindex,1));
    distance(ipnew,:,2)=squeeze(distance(ipnew,ipextraindex,2));
    % so now have points in order - ipextraindex has indices of the
    % ipextras 
    % now calculate
    distancetemp=sqrt(squeeze(distance2(ipnew,1:nneb)));
     if(iflag==0)
         iflag=1;
         distance0=distancetemp;
     else
         if (sum(abs(distancetemp-distance0))>tol)
             xohhelp=1
         end
     end
% pause on
% figure(2000)
% clf
% hold on
% plot(xyposextra(:,1),xyposextra(:,2),'rx')
% plot(xyposreq(:,1),xyposreq(:,2),'go')
% plot(xyposreq(ipnew,1),xyposreq(ipnew,2),'ro')
% for ij=1:nneb
%     plot(xyposextra(ipextraindex(ij),1),xyposextra(ipextraindex(ij),2),'bs')
% end
% pause;
% pause off
     for ij=1:nneb
         dx=distance(ipnew,ij,1);
         dy=distance(ipnew,ij,2);
         % index of neighbour
         f(1,ij)=1;
         f(2,ij)=dx;
         f(3,ij)=dy;
         f(4,ij)=dx*dx;
         f(5,ij)=dx*dy;
         f(6,ij)=dy*dy;
         f(7,ij)=dx*dx*dx;
         f(8,ij)=dx*dx*dy;
         f(9,ij)=dx*dy*dy;
         f(10,ij)=dy*dy*dy;
     end

     M=zeros(ni,ni);
     for ik=1:ni
         for ii=1:ni
             for ij=1:nneb
                 M(ik,ii)=M(ik,ii)+f(ik,ij)*f(ii,ij)/sigma(ij)^2;
             end
         end
     end
     pause on
     for iz=1:nz
% figure(2001)
% clf
% hold on
% % xlim([xyposreq(ipnew,1)-2*a/nscale,xyposreq(ipnew,1)+2*a/nscale])
% % ylim([xyposreq(ipnew,2)-2*a/nscale,xyposreq(ipnew,2)+2*a/nscale])
% plot3(xyposextra(:,1),xyposextra(:,2),V3Dextra(:,iz),'rx')
%     plot(xyposreq(ipnew,1),xyposreq(ipnew,2),'go')
% 
%     for ij=1:nneb
%     plot3(xyposreq(ipnew,1)+distance(ipnew,ij,1),xyposreq(ipnew,2)+distance(ipnew,ij,2),V3Dextra(ipextraindex(ij),iz),'bs')
%     end
% pause
         for ij=1:nneb
             D(ij)=V3Dextra(ipextraindex(ij),iz);
         end
         B=linspace(0,0,ni);
         for ik=1:ni
             for ij=1:nneb
                 B(ik)=B(ik)+D(ij)*f(ik,ij)/sigma(ij)^2;
             end
         end
         Avec=B/M;
         V3Dnew(np+ipnew,iz)=Avec(1);
     end
     pause off
end
npnew=npnew+np;


end
