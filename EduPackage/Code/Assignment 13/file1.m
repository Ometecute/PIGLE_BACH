% Potential for 157 pt. script
%% parameters

num_ang=24; % number of angles
num_pos=157; % number of potential positions
potential_mat=zeros(num_pos,num_ang);
new_file='Potential.mat';

%% position values

% top1=600;
% top2=470;
% top3=370;
% top4=270;

top1new=1000;
top2new=900;
top3new=800;
top4new=700;
top5new=600;
top6new=550;
top7new=500;
top8new=450;
top9new=400;
% top10new=600;
% top11new=600;
% top12new=600;

nearfccmin=85;
nearfccmax=250;

nearhcpmin=5;
nearhcpmax=65;
nearhcphighmax=100;

hcpmin=0;
hcpmax=60;
hcphighmax=80;

fccmin=75;
fccmax=210;

bridgemin=35;
bridgemax=300;

bridgebarmin=90;
bridgebarmax=340;

nearbridgebarmin=65; % positions near bridge sites
nearbridgebarmax=360; % positions near bridge sites

otherpoints=190;

%% positions

% top1pos=[1];
% top2pos=[2,6,136,114,91,84,62,32];
% top3pos=[3,7,10,13,140,110,118,124,100,95,88,58,66,39,36,33];
% top4pos=[4,8,11,14,17,146,143,137,115,121,127,131,106,103,98,92,85,63,69,72,43,40,37,34];
% 
% topnewpos=[1,   6,136,114,84,62,32,   2,13,110,91,58,39,   10,140,118,88,66,36,   7,17,146,137,115,124,95,85,63,72,43,33,   3,14,22,143,111,121,100,92,59,69,48,40,   11,20,149,141,119,127,98,89,67,75,46,37,   8,18,25,153,147,138,116,125,131,103,96,86,64,73,79,51,44,34,   4,15,23,28,151,144,112,122,129,106,101,93,60,70,77,54,49,41,   12,21,27,155,150,142,120,128,133,105,99,90,68,76,81,53,47,38,   9,19,26,30,157,154,148,139,117,126,132,135,108,104,97,87,65,74,80,83,56,52,45,35,   5,16,24,29,31,42,50,55,57,61,71,78,82,94,102,107,109,113,123,130,134,145,152,156]; 

topnewpos1=[1];
topnewpos2=[6,136,114,84,62,32];
topnewpos3=[2,13,110,91,58,39];
topnewpos4=[10,140,118,88,66,36];
topnewpos5=[7,17,146,137,115,124,95,85,63,72,43,33];
topnewpos6=[3,14,22,143,111,121,100,92,59,69,48,40];
topnewpos7=[11,20,149,141,119,127,98,89,67,75,46,37];
topnewpos8=[8,18,25,153,147,138,116,125,131,103,96,86,64,73,79,51,44,34];
topnewpos9=[4,15,23,28,151,144,112,122,129,106,101,93,60,70,77,54,49,41];
% topnewpos10=[12,21,27,155,150,142,120,128,133,105,99,90,68,76,81,53,47,38];
% topnewpos11=[9,19,26,30,157,154,148,139,117,126,132,135,108,104,97,87,65,74,80,83,56,52,45,35];
% topnewpos12=[5,16,24,29,31,42,50,55,57,61,71,78,82,94,102,107,109,113,123,130,134,145,152,156];

nearhcppos=[35,9,65,87,117,139];
nearfccpos=[30,157,83,56,108,135];

hcppos=[5,61,113];
fccpos=[31,57,109];

bridge1pos=[50,130];
bridge2pos=[24,102];
bridge3pos=[78,152];

bridge1barpos=[42,55,123,134];
bridge2barpos=[16,29,94,107];
bridge3barpos=[71,82,145,156];

nearbridge1pos=[45,47,52,126,128,132];
nearbridge2pos=[19,21,26,97,99,104];
nearbridge3pos=[74,76,80,148,150,154];

for i=1:num_ang/2
    for j=1:size(topnewpos1,2)
        pos_temp=topnewpos1(j);
        potential_mat(pos_temp,i)=top1new;
    end
    for j=1:size(topnewpos2,2)
        pos_temp=topnewpos2(j);
        potential_mat(pos_temp,i)=top2new;
    end
    for j=1:size(topnewpos3,2)
        pos_temp=topnewpos3(j);
        potential_mat(pos_temp,i)=top3new;
    end
    for j=1:size(topnewpos4,2)
        pos_temp=topnewpos4(j);
        potential_mat(pos_temp,i)=top4new;
    end
    for j=1:size(topnewpos5,2)
        pos_temp=topnewpos5(j);
        potential_mat(pos_temp,i)=top5new;
    end
    for j=1:size(topnewpos6,2)
        pos_temp=topnewpos6(j);
        potential_mat(pos_temp,i)=top6new;
    end
    for j=1:size(topnewpos7,2)
        pos_temp=topnewpos7(j);
        potential_mat(pos_temp,i)=top7new;
    end
    for j=1:size(topnewpos8,2)
        pos_temp=topnewpos8(j);
        potential_mat(pos_temp,i)=top8new;
    end
    for j=1:size(topnewpos9,2)
        pos_temp=topnewpos9(j);
        potential_mat(pos_temp,i)=top9new;
    end
%     for j=1:size(topnewpos10,2)
%         pos_temp=topnewpos10(j);
%         potential_mat(pos_temp,i)=top10new;
%     end
%     for j=1:size(topnewpos11,2)
%         pos_temp=topnewpos11(j);
%         potential_mat(pos_temp,i)=top11new;
%     end
%     for j=1:size(topnewpos12,2)
%         pos_temp=topnewpos12(j);
%         potential_mat(pos_temp,i)=top12new;
%     end
end

% topnewpos=[1,   6,136,114,84,62,32,   2,13,110,91,58,39,   10,140,118,88,66,36,  3,14,22,143,111,121,100,92,59,69,48,40,   11,20,149,141,119,127,98,89,67,75,46,37,   8,18,25,153,147,138,116,125,131,103,96,86,64,73,79,51,44,34,   4,15,23,28,151,144,112,122,129,106,101,93,60,70,77,54,49,41,   12,21,27,155,150,142,120,128,133,105,99,90,68,76,81,53,47,38,   9,19,26,30,157,154,148,139,117,126,132,135,108,104,97,87,65,74,80,83,56,52,45,35]; otherpoints=[]; 
% for i=1:157 
%     if sum(ismember([i],topnewpos))==1;
%         continue        
%     else
%         sizeotherpoints=size(otherpoints,2);
%         otherpoints(sizeotherpoints+1)=i; 
%     end
% end

%% angles

% all angles - top sites same for all angles

% for i=1:num_ang/2
%     for j=1:size(top1pos,2)
%         pos_temp=top1pos(j);
%         potential_mat(pos_temp,i)=top1;
%     end
%     for j=1:size(top2pos,2)
%         pos_temp=top2pos(j);
%         potential_mat(pos_temp,i)=top2;
%     end
%     for j=1:size(top3pos,2)
%         pos_temp=top3pos(j);
%         potential_mat(pos_temp,i)=top3;
%     end
%     for j=1:size(top4pos,2)
%         pos_temp=top4pos(j);
%         potential_mat(pos_temp,i)=top4;
%     end
% end

%%%%% 0, 60, 120 degrees (first inputting 0 degrees values)

%%% hcp, fcc

for i=1:size(hcppos,2)
    potential_mat(hcppos(i),1)=hcpmin;
end

for i=1:size(fccpos,2)
    potential_mat(fccpos(i),1)=fccmin;
end

%%% near hcp, fcc

for i=1:size(nearhcppos,2)
    potential_mat(nearhcppos(i),1)=nearhcpmin;
end

for i=1:size(nearfccpos,2)
    potential_mat(nearfccpos(i),1)=nearfccmin;
end

%%% bridge sites, barriers, near bridge sites

% bridge 1 

for i=1:size(bridge1pos,2)
    potential_mat(bridge1pos(i),1)=bridgemax;
end

for i=1:size(bridge1barpos,2)
    potential_mat(bridge1barpos(i),1)=bridgebarmax;
end

for i=1:size(nearbridge1pos,2)
    potential_mat(nearbridge1pos(i),1)=nearbridgebarmax;
end

% bridge 2

for i=1:size(bridge2pos,2)
    potential_mat(bridge2pos(i),1)=bridgemax;
end

for i=1:size(bridge2barpos,2)
    potential_mat(bridge2barpos(i),1)=bridgebarmax;
end

for i=1:size(nearbridge2pos,2)
    potential_mat(nearbridge2pos(i),1)=nearbridgebarmax;
end

% bridge 3

for i=1:size(bridge3pos,2)
    potential_mat(bridge3pos(i),1)=bridgemax;
end
    
for i=1:size(bridge3barpos,2)
    potential_mat(bridge3barpos(i),1)=bridgebarmax;
end

for i=1:size(nearbridge3pos,2)
    potential_mat(nearbridge3pos(i),1)=nearbridgebarmax;
end

%%% Same values for 60, 120 degrees

potential_mat(:,5)=potential_mat(:,1);
potential_mat(:,9)=potential_mat(:,1);

%%%%% 30, 90, 150 degrees (first inputting 30 degrees values)

%%% hcp, fcc

for i=1:size(hcppos,2)
    potential_mat(hcppos(i),3)=hcphighmax;
end

for i=1:size(fccpos,2)
    potential_mat(fccpos(i),3)=fccmax;
end

%%% near hcp, fcc

for i=1:size(nearhcppos,2)
    potential_mat(nearhcppos(i),3)=nearhcphighmax;
end

for i=1:size(nearfccpos,2)
    potential_mat(nearfccpos(i),3)=nearfccmax;
end

%%% bridge sites, barriers, near bridge sites

% bridge 1 

for i=1:size(bridge1pos,2)
    potential_mat(bridge1pos(i),3)=bridgemax;
end

for i=1:size(bridge1barpos,2)
    potential_mat(bridge1barpos(i),3)=bridgebarmax;
end

for i=1:size(nearbridge1pos,2)
    potential_mat(nearbridge1pos(i),3)=nearbridgebarmax;
end

% bridge 2

for i=1:size(bridge2pos,2)
    potential_mat(bridge2pos(i),3)=bridgemax;
end

for i=1:size(bridge2barpos,2)
    potential_mat(bridge2barpos(i),3)=bridgebarmax;
end

for i=1:size(nearbridge2pos,2)
    potential_mat(nearbridge2pos(i),3)=nearbridgebarmax;
end

% bridge 3

for i=1:size(bridge3pos,2)
    potential_mat(bridge3pos(i),3)=bridgemax;
end
    
for i=1:size(bridge3barpos,2)
    potential_mat(bridge3barpos(i),3)=bridgebarmax;
end

for i=1:size(nearbridge3pos,2)
    potential_mat(nearbridge3pos(i),3)=nearbridgebarmax;
end

%%% Same values for 90, 150 degrees

potential_mat(:,7)=potential_mat(:,3);
potential_mat(:,11)=potential_mat(:,3);

%%%%% 15, 45 degrees (first inputting 15 degrees values)

%%% hcp, fcc

for i=1:size(hcppos,2)
    potential_mat(hcppos(i),2)=hcpmax;
end

for i=1:size(fccpos,2)
    potential_mat(fccpos(i),2)=fccmax;
end

%%% near hcp, fcc

for i=1:size(nearhcppos,2)
    potential_mat(nearhcppos(i),2)=nearhcpmax;
end

for i=1:size(nearfccpos,2)
    potential_mat(nearfccpos(i),2)=nearfccmax;
end

%%% bridge sites, barriers, near bridge sites

% bridge 1 

for i=1:size(bridge1pos,2)
    potential_mat(bridge1pos(i),2)=bridgemin;
end

for i=1:size(bridge1barpos,2)
    potential_mat(bridge1barpos(i),2)=bridgebarmin;
end

for i=1:size(nearbridge1pos,2)
    potential_mat(nearbridge1pos(i),2)=nearbridgebarmin;
end

% bridge 2

for i=1:size(bridge2pos,2)
    potential_mat(bridge2pos(i),2)=bridgemax;
end

for i=1:size(bridge2barpos,2)
    potential_mat(bridge2barpos(i),2)=bridgebarmax;
end

for i=1:size(nearbridge2pos,2)
    potential_mat(nearbridge2pos(i),2)=nearbridgebarmax;
end

% bridge 3

for i=1:size(bridge3pos,2)
    potential_mat(bridge3pos(i),2)=bridgemax;
end
    
for i=1:size(bridge3barpos,2)
    potential_mat(bridge3barpos(i),2)=bridgebarmax;
end

for i=1:size(nearbridge3pos,2)
    potential_mat(nearbridge3pos(i),2)=nearbridgebarmax;
end

%%% Same values for 45 degrees

potential_mat(:,4)=potential_mat(:,2);

%%%%% 135, 165 degrees (first inputting 135 degrees values)

%%% hcp, fcc

for i=1:size(hcppos,2)
    potential_mat(hcppos(i),10)=hcpmax;
end

for i=1:size(fccpos,2)
    potential_mat(fccpos(i),10)=fccmax;
end

%%% near hcp, fcc

for i=1:size(nearhcppos,2)
    potential_mat(nearhcppos(i),10)=nearhcpmax;
end

for i=1:size(nearfccpos,2)
    potential_mat(nearfccpos(i),10)=nearfccmax;
end

%%% bridge sites, barriers, near bridge sites

% bridge 1 

for i=1:size(bridge1pos,2)
    potential_mat(bridge1pos(i),10)=bridgemax;
end

for i=1:size(bridge1barpos,2)
    potential_mat(bridge1barpos(i),10)=bridgebarmax;
end

for i=1:size(nearbridge1pos,2)
    potential_mat(nearbridge1pos(i),10)=nearbridgebarmax;
end

% bridge 2

for i=1:size(bridge2pos,2)
    potential_mat(bridge2pos(i),10)=bridgemin;
end

for i=1:size(bridge2barpos,2)
    potential_mat(bridge2barpos(i),10)=bridgebarmin;
end

for i=1:size(nearbridge2pos,2)
    potential_mat(nearbridge2pos(i),10)=nearbridgebarmin;
end

% bridge 3

for i=1:size(bridge3pos,2)
    potential_mat(bridge3pos(i),10)=bridgemax;
end
    
for i=1:size(bridge3barpos,2)
    potential_mat(bridge3barpos(i),10)=bridgebarmax;
end

for i=1:size(nearbridge3pos,2)
    potential_mat(nearbridge3pos(i),10)=nearbridgebarmax;
end

%%% Same values for 165 degrees

potential_mat(:,12)=potential_mat(:,10);

%%%%% 75, 105 degrees (first inputting 75 degrees values)

%%% hcp, fcc

for i=1:size(hcppos,2)
    potential_mat(hcppos(i),6)=hcpmax;
end

for i=1:size(fccpos,2)
    potential_mat(fccpos(i),6)=fccmax;
end

%%% near hcp, fcc

for i=1:size(nearhcppos,2)
    potential_mat(nearhcppos(i),6)=nearhcpmax;
end

for i=1:size(nearfccpos,2)
    potential_mat(nearfccpos(i),6)=nearfccmax;
end

%%% bridge sites, barriers, near bridge sites

% bridge 1 

for i=1:size(bridge1pos,2)
    potential_mat(bridge1pos(i),6)=bridgemax;
end

for i=1:size(bridge1barpos,2)
    potential_mat(bridge1barpos(i),6)=bridgebarmax;
end

for i=1:size(nearbridge1pos,2)
    potential_mat(nearbridge1pos(i),6)=nearbridgebarmax;
end

% bridge 2

for i=1:size(bridge2pos,2)
    potential_mat(bridge2pos(i),6)=bridgemax;
end

for i=1:size(bridge2barpos,2)
    potential_mat(bridge2barpos(i),6)=bridgebarmax;
end

for i=1:size(nearbridge2pos,2)
    potential_mat(nearbridge2pos(i),6)=nearbridgebarmax;
end

% bridge 3

for i=1:size(bridge3pos,2)
    potential_mat(bridge3pos(i),6)=bridgemin;
end
    
for i=1:size(bridge3barpos,2)
    potential_mat(bridge3barpos(i),6)=bridgebarmin;
end

for i=1:size(nearbridge3pos,2)
    potential_mat(nearbridge3pos(i),6)=nearbridgebarmin;
end

%%% Same values for 105 degrees

potential_mat(:,8)=potential_mat(:,6);


%%%%% All other points on matrix

for i=1:num_pos
    for j=1:num_ang
        if potential_mat(i,j)==0
            potential_mat(i,j)=otherpoints;
        end
    end
end

potential_mat(:,((num_ang/2)+1):num_ang)=potential_mat(:,1:(num_ang/2));

%%%%% Create new file

V3Dinterp_AQ=potential_mat;
general_file='Potential_AQ_base_file.mat';
copyfile(general_file, new_file);
save(new_file, 'V3Dinterp_AQ', '-append');
clear;
