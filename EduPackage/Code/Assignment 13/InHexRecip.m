function weight=InHexRecip(G,G1,n)
x=G(1);
y=G(2);
del=G1/1e4;
r3=sqrt(3);
weight=0;
if ( (y+del>-n*G1/2) && (y-del<n*G1/2) && (y-del<n*G1-r3*x) && (y+del>r3*x-n*G1) && (y+del>-r3*x-n*G1) && (y-del<r3*x+n*G1))
% we have a point inside the Wigner-Seitz cell
    weight=1;
% now need to check if its on an edge or corner

    dum=0;
    if( abs(y+n*G1/2)<del)
        dum=dum+1;
    end
    
    if( abs(y-n*G1/2)<del)
        dum=dum+1;
    end
    
    if( abs(y-n*G1+r3*x)<del)
        dum=dum+1;
    end
    
    if( abs(y-r3*x+n*G1)<del)
        dum=dum+1;
    end
    
    if( abs(y+r3*x+n*G1)<del)
        dum=dum+1;
    end
    
    if( abs(y-r3*x-n*G1)<del)
        dum=dum+1;
    end
    if(dum==1)
        % its on an edge
        weight=1/2;
    end
    if(dum==2)
        % its at a corner
        weight=1/3;
    end

end

