clear all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% For Pigle by rewriting EduPackage Assignment 13 file1/file2_interp
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Potential Values
a=2.461; ac=2.461/2/cos(pi/6); % lattice vectors (in Angström)

% coronene 9x9 graphite (tl-9x9-co-30.x)
% pot_vals =[-2.272  -2.228  -2.169  -1.986]*(1000);

% coronene 9x9 bilayer (bl-9x9-co-30.x)
% pot_vals = [-2.244  -2.209  -2.110  -1.952]*(1000);


% coronene 9x9 graphite incl. rot [tl-9x9-co-1.x;tl-9x9-co-30.x]
pot_vals =[-2.124  -2.136  -2.164  -2.208;  -2.272  -2.228  -2.169  -1.986]*(1000);
% coronene 9x9 bilayer incl. rot [bl-9x9-co-1.x;bl-9x9-co-30.x]
% pot_vals =[-2.081  -2.089  -2.153  -2.182;  -2.244  -2.209  -2.110  -1.952]*(1000);


pot_vals = pot_vals'-min(pot_vals(:)); % for differences ( in meV)


%% parameters

num_ang= 2; % number of angles
thetavec = [0,30];
num_pos=19; % number of potential positions
xypos=zeros(num_pos,2);
potential_mat=zeros(num_pos,num_ang);
new_file='Potential.mat';
N = [200,120,9];%x,y,theta result


if mod(N(3),2) ~= 1, error('theta dim should be a non-even Integer'),end
%%  map positions to new system (given pot_vals_pos to posittion: 1->1, 2->2, 5->3, 7->4)
position=[7,5,3,8,2,9, 8,2,9, 5,3,6, 4,2,1, 5,3,2, 6];
%       3
%     2   2
%   9   5   9
%     8   8
%   6   7   6
%     5   5
%   3   4   3
%     2   2
%       1

%position=[1,4,7,2,5,3, 2,5,3, 8,7,5, 6,5,3, 8,7,5, 5];
%       7
%     5   5
%   3   4   3
%     2   2
%   5   1   5
%     8   8
%   7   6   7
%     5   5
%       3

% remapping to match to pot_vals
old = [1,2,3,4,5,6,7,8,9];
new = [1,2,1,5,5,2,7,5,1];
position = changem(position,new, old);

%remapping to match to new numbering
old = [1,2,5,7];
if ~all(sum(position==old',1)), error('Position remapping is incomplete'), end
new = [1,2,3,4];
position = changem(position,new, old);


for temp_ang=1:num_ang
    for temp_pos=1:num_pos
       potential_mat(temp_pos,temp_ang)=pot_vals(position(temp_pos),temp_ang);
    end
end
%%
ynew = [0:4]*ac/4;
x_new = [0:2]*a/4;
xypos = [x_new(1),ynew(1);x_new(1),ynew(3);x_new(1),ynew(5);x_new(2),ynew(2);...
    x_new(2),ynew(4);x_new(3),ynew(3);...
    -x_new(2),ynew(2);-x_new(2),ynew(4);-x_new(3),ynew(3);...
    -x_new(2),-ynew(2);-x_new(3),-ynew(3);-x_new(3),-ynew(1);...
    -x_new(1),-ynew(3);-x_new(2),-ynew(4);-x_new(1),-ynew(5);...
    x_new(2),-ynew(2);x_new(3),-ynew(3);x_new(2),-ynew(4);...
    x_new(3),-ynew(1)];

tmp = zeros([size(xypos),num_ang]);
for i = 1:num_ang
    tmp(:,:,i) =xypos;
end
xypos = tmp;
%filling out a full unitcell
x_tmp = [0,1,2]*a*sin(pi/6);
y_tmp = [0,1,2]*a*cos(pi/6);
[x_m,y_m] = meshgrid(x_tmp,y_tmp);
x_trans = x_m([1,3,5,7,9])';
y_trans = y_m([1,3,5,7,9])';
trans_vec = [x_trans,y_trans]; 

xypos_tmp = xypos; potential_mat_tmp = potential_mat;
for trans = 2:length(trans_vec)
    xypos_tmp = [xypos_tmp;xypos+trans_vec(trans,:)];
    potential_mat_tmp = [potential_mat_tmp;potential_mat];
end
xypos = xypos_tmp; potential_mat = potential_mat_tmp;
% plot xypos
figure(4)
plot(xypos(:,1),xypos(:,2),'rx')
hold on
for i=1:length(position)
    text(xypos(i,1),xypos(i,2),num2str(position(1,i)))
end
hold off


%%



x_n = linspace(0,a,N(2));
y_n = linspace(0,2*a*cos(pi/6),N(1));
theta_n = linspace(thetavec(1),thetavec(2)*2,N(3));

potential_mat = interp1(thetavec,potential_mat',theta_n(1:(N(3)+1)/2))';

for i = 1:(N(3)+1)/2
    [X,Y]=meshgrid(x_n,y_n);
    pe3D(:,:,i)=griddata(xypos(:,1,1),xypos(:,2,1),potential_mat(:,i),X,Y,'cubic');
end
pe3D = cat(3,pe3D,flip(pe3D(:,:,1:end-1),3));
for i = 1:N(3)
    figure(100+i);
    %stest = surf(X,Y,pe3D(:,:,i)); set(stest,'LineStyle','none'); caxis([color_scale_min color_scale_max]); view([0 90]); daspect(asp); xlim(Lx); ylim(Ly);
    stest = contourf(X,Y,pe3D(:,:,i),linspace(0,300,25)); axis equal; title(num2str(theta_n(i)));
    clim([0,300])
    %exportgraphics(gcf,sprintf('deg%d.png',round(thetavec(i)*3*(180/pi))),'Resolution','300')
end
PES4D_AQ = zeros([N(1),N(2),1,N(3)]);
PES4D_AQ(:,:,1,:) = pe3D;
save_file='loadPES_cor.mat';
save(save_file, 'PES4D_AQ');

% %%
% weight = [1,1,1/3,1,1/2,1/3, 1,1/2,1/3, 1,1/3,1/2, 1,1/2,1/3, 1,1/3,1/2, 1/2]';
% thetapos = [0];
% 
% potential_mat(:,((num_ang/2)+1):num_ang)=potential_mat(:,1:(num_ang/2));
% 
% %%%%% Create new file
% 
% V3Dinterp_AQ=potential_mat;
% 
% xn=linspace(min(xypos(:,1)),max(xypos(:,1)),201);
% yn=linspace(min(xypos(:,2)),max(xypos(:,2)),201);
% [xnew,ynew]=meshgrid(xn,yn);
% znew=griddata(xypos(:,1),xypos(:,2),V3Dinterp_AQ,xnew,ynew,'cubic');
% 
% figure
% mesh(xnew,ynew,znew)
% general_file='Potential_AQ_base_file.mat';
% save(general_file, 'V3Dinterp_AQ',"xypos","thetapos","weight");
% clear;


