classdef (Abstract) MCSimu < handle
    
    properties
        % length of each timestep /microseconds
        timestep = 1;

        % Probability of a particle jumping one space in one timestep
        hopProb = 0.01;

        % lattice parameter /A
        lattConst = 2.5;

        % 3 vectors recording the trajectory of the particle
        x = [];
        y = [];
        t = [];

        % maximum no. timesteps allocated for the simulation, including t=0 state
        numSteps = 10;
    end
    
    
    
    methods(Static, Abstract)
        % abstract − samples range of integers corresponding to jump directions
        hopDirec()
    end

    
    
    methods (Abstract)
        % abstract − propagates the state of the simulation through n timesteps of length this.timestep
        nStates_loop(obj, n)
    end
    
    
    
    methods
        
        function obj = MCSimu(hopProb, lattConst, timestep, numSteps)
            % CONSTRUCTOR

            % If the correct no. of arguments are passed in...
            if nargin > 0

                % Copy passed in arguments of constructor to their corresponding properties:

                % jump probability (in any direction) at each timestep
                obj.hopProb = hopProb;

                % lattice spacing of grid /A
                obj.lattConst = lattConst;

                % length of each timestep /ps
                obj.timestep = timestep;

                % maximum number of timesteps that simulation can store
                obj.numSteps = numSteps;

                % Initialise the x,y,t vectors with max length of numSteps:
                obj.x = zeros(numSteps,1);
                obj.y = zeros(numSteps,1);
                obj.t = zeros(numSteps,1);
            end

            % Otherwise, just use the default values
        end
        
        
        
        function h = hopOrNot(obj)
            % Single Bernoulli trial => gives 1 w/. prob p:
            h = binornd(1, obj.hopProb);
        end
        
        
        
        function plotAll(obj)
            % Plot all three graphs determining the trajectory.
            figure();
            plot(obj.t, obj.x);
            xlabel('t /micro s'); ylabel('x /A');

            figure();
            plot(obj.t, obj.y);
            xlabel('t /micro s'); ylabel('y /A');

            figure();
            plot(obj.x, obj.y); axis equal;
            xlabel('x /A'); ylabel('y /A');
        end
        
    end
    
end