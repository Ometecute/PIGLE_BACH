clear all; close all; fclose all

a=2.461; ac=2.461/2/cos(pi/6);

%0 deg rotation, down
En1=[-3.042324 -3.075258 -3.099149 -3.081047 -3.102110 -3.080679 -3.099965 -3.073353 -3.038074 ]; zticks1=[-3.12:0.01:-3.04];
%30 deg rotation, down
En2=[-3.076117 -3.071919 -3.059282 -3.086208 -3.080722 -3.080666 -3.099143 -3.082242 -3.080720 ]; zticks2=[-3.12:0.01:-3.04];
%60 deg rotation, down
En3=[-3.111292 -3.078379 -3.042995 -3.091814 -3.066175 -3.084313 -3.103534 -3.099570 -3.112259 ]; zticks3=[-3.12:0.01:-3.04]; 
Envec=[En1;En2;En3]; ztickvec=[zticks1;zticks2;zticks3];

Enrot=[En3-En2,En3-En1,En2-En1];
dErot=mean(abs(Enrot)*1000);


%En=En1;
% %60 deg rotation, up
Enup1=[-2.925250 -2.926801 -2.912945 -2.958340 -2.952534 -2.931609 -2.993975 -2.950563 -2.925241];
dEinv1=[Enup1-En1]*1000;
Enup2=[-2.904734 -2.934680 -2.957077 -2.939878 -2.960218 -2.932259 -2.963362 -2.939877 -2.901346];
dEinv2=[Enup2-En2]*1000;
Enup3=[-2.911527 -2.926229 -2.935856 -2.953167 -2.951907 -2.930252 -2.998072 -2.946639 -2.908063]; %zticks=[-3:0.01:-2.9]; 
dEinv3=[Enup3-En3]*1000;
    
x1=[-1,0,1,2,3].*ac*cos(pi/6);
x2=[-1,+1,+3,+5,+7]*ac/2*cos(pi/6);

for kk=1:16
    if mod(kk,2)==1
       xx(kk,:)=x1;
    else
       xx(kk,:)=x2;
    end
    yy(kk,:)=ones(size(x1))*ac/2*sin(pi/6)*(kk-1)-ac;
end
yy=yy-ac/2; xx=xx-ac*cos(pi/6);


figure1=figure('units','centimeters','position',[5,3,36,14],'color','white','DefaultTextInterpreter','LaTex','DefaultLineLineWidth',1.5);

for kk=1:3
    clear zz;
    En=Envec(kk,:);
    zticks=ztickvec(kk,:);
        
    zz(1,:)=[En(5),En(1),En(5)];
    zz(2,:)=[En(2),En(2),En(2)];
    zz(3,:)=[En(3),En(4),En(3)];
    zz(4,:)=[En(5),En(5),En(5)];
    zz(5,:)=[En(6),En(7),En(6)];
    zz(6,:)=[En(8),En(8),En(8)];
    zz(7,:)=[En(9),En(5),En(9)];
    zz(8,:)=[En(2),En(2),En(2)];
    zz(9,:)=[En(4),En(3),En(4)];
    zz(10,:)=[En(5),En(5),En(5)];
    zz(11,:)=[En(7),En(6),En(7)];
    zz(12,:)=[En(8),En(8),En(8)];
    zz(13,:)=[En(5),En(9),En(5)];
    zz(14,:)=[En(2),En(2),En(2)];
    zz(15,:)=[En(3),En(4),En(3)];
    zz(16,:)=[En(5),En(5),En(5)];
    
    zz=[zz,zz(:,2:3)];
    
    xn=linspace(min(xx(:)),max(xx(:)),201);
    yn=linspace(min(yy(:)),max(yy(:)),201);
    [xnew,ynew]=meshgrid(xn,yn);
    znew=griddata(xx,yy,zz,xnew,ynew,'cubic');
    
    subplot(1,4,kk)
    set(gca,'LineWidth',1.5,'FontSize',14); box(gca,'on'); hold(gca,'all'); set(gca, 'Layer','top');
    set(gca,'TickLabelInterpreter','LaTex');
    %surf(xnew,ynew,znew);
    %shading interp;
    %axis square; view([0,90]);
    contourf(xnew,ynew,znew,'LineStyle','none','LevelStep',1e-3); axis square;
    caxis([min(Envec(:)),max(Envec(:))]);
    if kk==1
            ylabel('$y$ (\AA)','Interpreter','LaTex');
    elseif kk==3
        cb=colorbar;  %cb.TickLabels = '%.3f';
        set(cb,'YTick',zticks,'TickLabelInterpreter', 'LaTex');
        for nn=1:length(cb.Ticks)
            tlabel{nn}=num2str(cb.Ticks(nn),'%.2f');
        end
        cb.TickLabels=tlabel;
        set(cb,'LineWidth',1.5,'FontSize',16,'ticklength',1.5*get(cb,'ticklength'));
        ax=gca; cpos=cb.Position; axpos = ax.Position;
        cpos(3) = 0.7*cpos(3); cpos(1)=1.1*cpos(1); 
        cpos(2)=cpos(2)*0.82;  cpos(4)=cpos(2)*1.35;
        cb.Position = cpos; ax.Position=axpos;
        
        set(gca,'TickDir','in','TickLabelInterpreter', 'LaTex');
        set(gca,'ticklength',1.5*get(gca,'ticklength'))
        ylabel(cb,'Adsorption energy $E_{a}$ (eV)','Interpreter','LaTex');
    end
    xlabel('$x$ (\AA)','Interpreter','LaTex');
    title([num2str((kk-1)*30),'$^{\circ}$ rotation']);

    
    
    %plot the corresponding hexagonal grid
    r=ac; c=r;
    v=30:60:390;
    
    cv=r*cosd(v); sv=r*sind(v)+ac/2;
    for y=-1:2:1
        for w=0:2:2
            line(w*sqrt(3)/2*c+cv,y*1.5*c+sv,'tag','h','Color',[1,0.5,0]);
        end
    end
    for y=0:2:2
        for w=-1:2:1
            line(w*sqrt(3)/2*c+cv,y*1.5*c+sv,'tag','h','Color',[1,0.5,0]);
        end
    end
    
    cv=r*cosd(v)+ac*cos(pi/6); sv=r*sind(v)-ac/2;
    for y=0:2:2
        for w=-2:2:0
            line(w*sqrt(3)/2*c+cv,y*1.5*c+sv,'tag','h','Color',[1,0,0]);
        end
    end
    for y=-1:2:1
        for w=-1:2:1
            line(w*sqrt(3)/2*c+cv,y*1.5*c+sv,'tag','h','Color',[1,0,0]);
        end
    end
    
    
    xlim([min(xx(:)),a]);
    ylim([min(yy(:)),2*ac]);
    
end

%set(gcf,'PaperPositionMode','auto')
%print('-dpng','DFT_PES','-r300');
%print('-depsc','DFT_PES_30up');
% % %print('-dsvg','DFT_PES');