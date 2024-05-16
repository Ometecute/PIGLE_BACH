classdef auxCls
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties(Constant)

    end
 
    methods(Static)
        
        function [a1, a2, b1, b2] = unitCellVecs(surfaceType,lttcCnst, rotation)
            % UNITCELLVECS creates the vectors of a unit cell (Real & Reciprocal space)
            % Currently, two types of symmetries are supported, hexagonal and fcc. Also. tt is possible to define rotation of the unit cell.
            % Input:   surfaceType - 1-hexagonal 2-fcc
            %         lttcCnst - Lattice constant
            %         rotation - either a 2x2 rotational Matrix, or the rotation in degrees. 
            % Output: {a1,a2}  - real space vectors
            %         {b1,b2}  - reciprocal vectors
            
            % Substrate vectors
            if length(surfaceType) == 3
                angle=surfaceType(3);
                Rmat = [cosd(angle) sind(angle); -sind(angle) cosd(angle)];                
                a1s = [1 0]*surfaceType(1);
                a2s = [Rmat*[1 0]']'*surfaceType(2);
            elseif surfaceType == 1 % hexagonal
                a1s=[1 0]*lttcCnst;
                a2s=[-0.5 sqrt(3)/2]*lttcCnst;    
            elseif surfaceType == 2 % fcc(100)
                a1s=[1 0]*lttcCnst;
                a2s=[0 1]*lttcCnst;
            end
            
            if length(rotation)==1,%assume rotation angle in degrees
                Rmat = [cosd(rotation) sind(rotation); -sind(rotation) cosd(rotation)];
                a1 = [Rmat*a1s']';
                a2 = [Rmat*a2s']';
            else
                % Adsorbate vectors
                a1 = rotation(1,1)*a1s + rotation(1,2)*a2s;
                a2 = rotation(2,1)*a1s + rotation(2,2)*a2s;
            end
            
            b1=2*pi*cross([a2 0],[0 0 1])/([a1 0]*[cross([a2 0],[0 0 1])]');
            b1=b1(1:2);
            b2=2*pi*cross([0 0 1],[a1 0])/([a2 0]*[cross([0 0 1],[a1 0])]');
            b2=b2(1:2);
        end
        
        function [d1,d2]=real2recip(c1,c2,toRecip)
            if toRecip
                d2=2*pi*cross([0 0 1],[c1 0])/norm(cross([c1 0],[c2 0]));
                d2=d2(1:2);
                d1=2*pi*cross([c2 0],[0 0 1])/norm(cross([c1 0],[c2 0]));
                d1=d1(1:2);
            else
                disp('real2recip: not implemented yet ...')
            end
        end
        
        function  R = sampleUnitCell(len,a1,a2)
            % Used to be called uniSpacedPointsInUnitCell
            % Uniformly distributed over the unit cell
            vec=0:1/len:1-1/len;
            len=length(vec);
            [N,M]=meshgrid(vec);
            for n=1:len
                for m=1:len
                    R((n-1)*len+m,1:2)=N(n,m)*a1+M(n,m)*a2;
                end
            end
        end   
        
        function G = spanVecs(b1,b2,nVec1,nVec2)

            [N_g,M_g]=meshgrid(nVec2,nVec1);
            Gx = N_g * b1(1) + M_g * b2(1);
            Gy = N_g * b1(2) + M_g * b2(2);
            G = [Gx(:) Gy(:)];

        end
        
    end
    
end