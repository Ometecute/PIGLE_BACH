classdef Grid < handle
    
    properties(Constant)
        % scaled Boltzmann const. /kCal K^{-1} so timestep =~ 2.5 fs
        k = 0.5;

        % activation energy for diffusion at 0 coverage /kcal mol^{-1}
        E0 = 5.2;

        % nearest neighbour interaction energy /kcal mol^{-1}
        e = 5.0;

        % Arrhenius prefactor for hopping frequency
        nu = 1;
	end
    
    
    
    properties
        % substrate lattice parameter /A
        a = 3.912;

        % adsorbate adlayer depth /A
        c = 1;

        % temperature /K
        T = 300;

        % MC grid horizontal dimensions, same for u and v
        len = 100;

        % MC grid vertical dimensios, w
        hgt = 100;

        % x,y,z grid of lattice sites
        Rxyz = Grid.xyGrid(3.912, 1, 100, 100);

        % INDEX (u,v,w) positions of adsorbates (1 if present, 0 if not)
        rs = zeros([100 100 100]);
    end
    
    
    
    methods(Static)
        
        function Rxyz = xyGrid(a, c, len, hgt)
            % generates a hexagonal grid in orthonormal x,ycoords

            % Coordinates in hexagonal basis:
            [U,V,W] = meshgrid(a.*[0:1:len-1], a.*[0:1:len-1], c.*[0:1:hgt-1]);

            % Coordinates in orthonormal rectangular basis:
            X = U+V.*1/2; Y = V.*sqrt(3)/2; Z = W;
            Rxyz = [X(:) Y(:) Z(:)];
        end
        
    end
    
    
    
    methods
        
        function this = Grid(a, c, T, len, hgt)
            % CONSTRUCTOR

            % if the number of arguments is correct, run the following, else just use the default values
            if nargin > 0

            % substrate lattice parameter /A
            this.a = a;

            % adsorbate adlayer depth /A
            this.c = c;

            % temperature /K
            this.T = T;

            % MC grid horizontal dimensions, same for u and
            this.len = len;

            % MC grid vertical dimensios, w
            this.hgt = hgt;

            % x,y,z grid of lattice sites
            this.Rxyz = Grid.xyGrid(a, c, len, hgt);

            % INDEX (u,v,w) positions of adsorbates (1 if present, 0 if not)
            this.rs = zeros([this.len this.len this.hgt]);
            end
        end
            
            
        
        function fz = pbc(this,z)
            % applies PBCs to any scalar/vector/matrix input of indices
            fz = zeros(size(z));

            for i = 1:length(z)
                % indices in range [1, len]
                fz(i) = mod(z(i) - 1, this.len) + 1;
            end
        end
        
        
        
        function addParticle(this, u, v)
            % apply PBCs
            u = this.pbc(u); v = this.pbc(v);

            % find lowest empty site:
            w=1;
            while this.rs(u,v,w)>0 && w<this.hgt
                w = w+1;
            end

            % i.e. run out of vertical space
            if (w==this.hgt && this.hgt~=1) || (this.rs(u,v,this.hgt)==1 && this.hgt==1)
                fprintf('Location already fully occupied.');
                % space left -> add it
            else
                this.rs(u,v,w) = 1;
            end
        end
        
        
        
        function addRandParticle(this)
            % tries=current xy attempt, w = current z attempt
            tries=1; w=1;

            
            while true
                % randomly pick horizontal position
                u = randi([1 this.len]);
                v = randi([1 this.len]);
                % start at bottom layer

                w=1;

                % Try all vertical sites
                while this.rs(u,v,w)>0 && w<this.hgt
                    w = w+1;
                end

                % found empty slot or run out of attempts
                if tries>100 || this.rs(u,v,w)==0, break; end
                tries = tries+1;
            end


            % if run out of attempts, print error message, otherwise add particle
            if tries>100
                fprintf('No available site found in 100 random probes.\n');
            else
                this.rs(u,v,w) = 1;
            end
        end
        
        
        
        function [nn, fs] = neighbours(this, u, v, w)
            % nn = count of nearest neighbours of given point
            % fs = matrix of empty neighbours (possible next positions)


            % PBCs
            us = this.pbc([u-1 u u+1]);
            vs = this.pbc([v-1 v v+1]);
            
            fs = [];

            % left
            if this.rs(us(1),vs(2),w) == 0
                fs = [fs; [us(1) vs(2)]];
            end

            % right
            if this.rs(us(3),vs(2),w) == 0
                fs = [fs; [us(3) vs(2)]];
            end

            % down
            if this.rs(us(2),vs(1),w) == 0
                fs = [fs; [us(2) vs(1)]];
            end

            % up
            if this.rs(us(2),vs(3),w) == 0
                fs = [fs; [us(2) vs(3)]];
            end

            % max neighbours - no. empties
            nn = 4 - size(fs, 1);
        end
        
        
        
        function timestep(this)
            % allows each particle to advance one timestep. The length of each timestep is determined by the arbitrary constant nu

            % propagate the state of each layer sequentially, from the bottom up
            for w = 1:this.hgt

                % iterate over horizontal positions
                for u = 1:this.len

                    % iterate over horizontal positions
                    for v = 1:this.len

                        % if an adsorbate is present, we must consider whether it jumps
                        if this.rs(u,v,w) > 0

                            % calculates n_i, the number of nearest neighbours and fs, the index coordinates of the available empty sites
                            [n_i, fs] = this.neighbours(u,v,w);

                            % calculate total energy of adsorbate
                            En_i = Grid.E0 + n_i*Grid.e;

                            % calculate hopping probability
                            w_if = Grid.nu * exp(-En_i/(Grid.k* this.T));

                            % If it can move...
                            if ~isempty(fs) && binornd(1,w_if)
                                % no. of empty sites
                                frows = size(fs,1);

                                % sample empty sites
                                f = fs(randperm(frows,1),:);

                                % Delete old position
                                this.rs(u,v,w) = 0;

                                % Try same layer @ new u,v posn
                                w_new = w;
                                while w_new~=1 && this.rs(f(1),f(2),w_new-1)==0
                                    % If there's room, 'fall' down
                                    w_new = w_new - 1;
                                end

                                % Concretely move to lower position @ new u,v
                                this.rs(f(1),f(2),w_new) = 1;
                                
                            % If it can't but the site below is empty...
                            elseif w~=1 && this.rs(u,v,w-1)~=1
                                % Delete old position
                                this.rs(u,v,w) = 0;

                                % Start at current layer
                                w_new = w;
                                while w_new~=1 && this.rs(f(1),f(2),w_new-1)==0
                                    % If there's room, 'fall' down
                                    w_new = w_new - 1;
                                end

                                % Concretely move to lower position @ same u,v
                                this.rs(u,v,w_new) = 1;
                                
                            end
                        end
                    end
                end
            end
        end
        
        
        
        function plot(this)
            % Reshape rs to match the meshgrid index convention
            rs_flat = reshape(this.rs, [this.len^2*this.hgt 1]);

            % Use find() to return all the indices for which there is an adsorbate
            ind = find(rs_flat > 0);
            figure();

            % Plot surface
            fill(10.*[0 this.a 3*this.a/2 this.a/2]', 10.*[0 0 this.a*sqrt(3)/2 this.a*sqrt(3)/2]', [0.9 0.9 0.9], 'LineWidth', 0.01); hold on;

            % Plot points at locations of adsorbates
            scatter3(this.Rxyz(ind,1), this.Rxyz(ind,2), this.Rxyz(ind,3), '.'); axis equal;


            % Set limits of plot to encompass the entire surface we're simulating, as well as all possible adlayer heights
            xlim([0 this.len*this.a*3/2]); ylim([0 this.len*this .a*sqrt(3)/2]);
            zlim([0 this.len*this.a]);
        end
        
    end
end