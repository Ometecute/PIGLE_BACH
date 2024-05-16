classdef Solenoid < handle
    %Producing a set of simple solenoid fields to be used for spin precession simulations
    
    %Detailed explanation goes here
    
    properties
        l %Length of solenoid in z-direction
        lx %Length of solenoid in x-direction
        ly %Length of solenoid in y-direction
        
        x % vectors of considered x-coords
        y % vectors of considered y-coords
        z % vectors of considered z-coords
    end
    
    methods
        
        function obj = Solenoid(l, lx, ly, x, y, z)
            if nargin > 0
                obj.l = l;
                obj.lx = lx;
                obj.ly = ly;
                
                % Creating three 3 row vectors (called 'x', 'y' and 'z') that 
                % contain the x-, y- and z-coordinates of all of the points at which 
                % we are going to generate the B-field
                obj.x = x;
                obj.y = y;
                obj.z = z;
                
                %When we generate a matrix, B, containing the magnitude of one of
                %the components of the B field, the value B(i,j,k) will give the B
                %field at position (x(i), y(j), z(k))
            end
        end
        
        function [Bz_box] = box_z(obj)
            %Creating a Bz field that is equal to 1 inside a rectangular solenoid and 0 elsewhere
            
            %The centre of the solenoid is at x=y=z=0
            
            %Using grid function (above) to obtain range of x, y and z coordinates
            %at which to evaluate B field
          
            
            %A 3D matrix which will be populated with the size of the
            %Bz-field at each point in space
            Bz_box = zeros([size(obj.x, 2) size(obj.z, 2)]);
            
            
            %Assigning a value of 1 to those points in the Bz_box grid that are inside the solenoid
            
            for ij = 1:size(obj.x,2)
                for k = 1:size(obj.z,2)
                    if obj.x(ij) < (obj.lx/2) && obj.x(ij) > -(obj.lx/2) && obj.y(ij) < (obj.ly/2) && obj.y(ij) > -(obj.ly/2) && obj.z(k) < (obj.l/2) && obj.z(k) > -(obj.l/2)
                        Bz_box (ij,k) = 1;
                    end
                end
            end
        end
        
        function [Bz_axial_decay] = axial_decay(obj, axial_decay_const)
            %Creating a Bz field that is equal to 1 in the middle of the
            %solenoid and decays exponentially towards the edges with a decay
            %constant given by the property axial_decay_const

            %The field is still 0 when you go outside of the solenoid in the x
            %and y directions

            %The centre of the solenoid is at x=y=z=0

            %Using grid function (above) to obtain range of x, y and z coordinates
            %at which to evaluate B field


            %A 3D matrix which will be populated with the size of the
            %Bz-field at each point in space
            Bz_axial_decay = zeros([size(obj.x, 2) size(obj.z, 2)]);

            %Adding the exponential decay

            for ij = 1:size(obj.x,2)
                for k = 1:size(obj.z,2)
                   if obj.x(ij) < (obj.lx/2) && obj.x(ij) > -(obj.lx/2) && obj.y(ij) < (obj.ly/2) && obj.y(ij) > -(obj.ly/2)
                        Bz_axial_decay (ij,k) = exp(-abs(obj.z(k))*axial_decay_const);

                   end 
                end
            end


        end

        function [Bz_radial_decay] = radial_decay(obj, radial_decay_const)
            %Creating a Bz field that is equal to 1 on axis and decays
            %exponentially as you move away from the axis in either the x
            %or y directions

            %The decay constant in both directions is given by the property
            %radial_decay_const

            %The field is still 0 when you go outside the solenoid in the z
            %direction

            %The centre of the solenoid is at x=y=z=0

            %Using grid function (above) to obtain range of x, y and z coordinates
            %at which to evaluate B field


            %A 3D matrix which will be populated with the size of the
            %Bz-field at each point in space
            Bz_radial_decay = zeros([size(obj.x, 2) size(obj.z, 2)]);

            %Adding the exponential decay

            for ij = 1:size(obj.x,2)
                for k = 1:size(obj.z,2)
                   if obj.z(k) < (obj.l/2) && obj.z(k) > -(obj.l/2)
                        Bz_radial_decay (ij,k) = exp(-abs(obj.x(ij))*radial_decay_const - abs(obj.y(ij))*radial_decay_const);
                   end
                end
            end

        end
        
        function [Bz_twodecay] = twodecay(obj, axial_decay_const, radial_decay_const)
            %Combining radial and axial decay
            
            %The centre of the solenoid is at x=y=z=0
            
            %Using grid function (above) to obtain range of x, y and z coordinates
            %at which to evaluate B field

            %A 3D matrix which will be populated with the size of the
            %Bz-field at each point in space
            Bz_twodecay = zeros([size(obj.x, 2) size(obj.z, 2)]);

            %Adding the exponential decay

            for ij = 1:size(obj.x,2)
                for k = 1:size(obj.z,2)
                    Bz_twodecay (ij,k) = exp(-abs(obj.x(ij))*radial_decay_const - abs(obj.y(ij))*radial_decay_const - abs(obj.z(k))*axial_decay_const);
                end
            end
            
        end
        
        
        function [Bx_box,By_box] = box_xy(obj)
            %Creating a Bx and By field that is zero everywhere
            
            %The centre of the solenoid is at x=y=z=0
            
            %3D matrices which will be populated with the size of the
            %Bx and By-field at each point in space
            Bx_box = zeros([size(obj.x, 2) size(obj.z, 2)]);
            By_box = zeros([size(obj.x, 2) size(obj.z, 2)]);
        end
        
        function [Bx_box,By_box] = xy_tapered(obj)
            %Creating a Bx and By field that is tapered
            %Bx and By components are tapered in opposite senses
            
            %The centre of the solenoid is at x=y=z=0    
            
            %3D matrices which will be populated with the size of the
            %Bx and By fields at each point in space
            Bx_box = zeros([size(obj.x, 2) size(obj.z, 2)]);
            By_box = zeros([size(obj.x, 2) size(obj.z, 2)]);
            
            %Adding the taper
            
            %Defining k_taper, where k_taper = 2pi/wavelength of taper
            %Using 2x the length of the solenoid as the taper wavelength
            
            k_taper = (2*pi)/(2*obj.l);
            
            for ij = 1:size(obj.x,2)
                for k = 1:size(obj.z,2)
                    if obj.x(ij) < (obj.lx/2) && obj.x(ij) > -(obj.lx/2) && obj.y(ij) < (obj.ly/2) && obj.y(ij) > -(obj.ly/2) && obj.z(k) < (obj.l/2) && obj.z(k) > -(obj.l/2)
                        Bx_box(ij,k) = (cos(obj.z(k)*k_taper))^2;
                        By_box(ij,k) = (cos(obj.z(k)*k_taper + pi/2))^2;
                    end
                end
            end
        end
        
       function [Bz_radial_decay] = gaussian_radial_decay(obj, radial_decay_const)
            %Creating a Bz field that is equal to 1 on axis and decays with
            %a gaussian distribtuion in both the x and y directions

            %The decay constant in both directions is given by the property
            %radial_decay_const

            %The field is still 0 when you go outside the solenoid in the z
            %direction

            %The centre of the solenoid is at x=y=z=0

            %Using grid function (above) to obtain range of x, y and z coordinates
            %at which to evaluate B field


            %A 3D matrix which will be populated with the size of the
            %Bz-field at each point in space
            Bz_radial_decay = zeros([size(obj.x, 2) size(obj.z, 2)]);

            %Adding the exponential decay

            for ij = 1:size(obj.x,2)
                for k = 1:size(obj.z,2)
                   if obj.z(k) < (obj.l/2) && obj.z(k) > -(obj.l/2)
                        Bz_radial_decay(ij,k) = exp(-(obj.x(ij)^2)*radial_decay_const - (obj.y(ij)^2)*radial_decay_const);
                    end 
                end
            end

       end
        
       
       function [Bz] = smoothstep_z(obj, smoothLen)
            %Creating a Bz field that is equal to 1 inside a rectangular
            %solenoid and 0 elsewhere with a smooth transtion between
            
            %The centre of the solenoid is at x=y=z=0
            
            %A 3D matrix which will be populated with the size of the
            %Bz-field at each point in space
            Bz = zeros([size(obj.x, 2) size(obj.z, 2)]);
            
            %Assigning a value of 1 to those points in the Bz_box grid that are inside the solenoid, smoothing out around the edges
            
            for ij = 1:size(obj.x,2)
                for k = 1:size(obj.z,2)

                    if obj.x(ij) < (obj.lx/2) && obj.x(ij) > -(obj.lx/2) && obj.y(ij) < (obj.ly/2) && obj.y(ij) > -(obj.ly/2)
                        if obj.z(k) < (obj.l/2) && obj.z(k) > -(obj.l/2)
                            Bz(ij,k) = 1;
                        elseif obj.z(k) > -(obj.l/2 + smoothLen) && obj.z(k) <= -(obj.l/2)
                            Bz(ij,k) = 3*(1 + (obj.z(k)+obj.l/2)/smoothLen)^2 - 2*(1 + (obj.z(k)+obj.l/2)/smoothLen)^3;
                        elseif obj.z(k) >= (obj.l/2) && obj.z(k) < (obj.l/2 + smoothLen)
                            Bz(ij,k) = 3*(1 - (obj.z(k)-obj.l/2)/smoothLen)^2 - 2*(1 - (obj.z(k)-obj.l/2)/smoothLen)^3;
                        end
                    end

                end
            end
       end
        
    end
end

