classdef MCSimu_sqr < MCSimu
    
    methods(Static)
        
        function d = hopDirec()
            % Define: 1 up, 2 right, 3 down, 4 left
            d = randi(4);
        end
        
    end
    
    
    
    methods    
        
        function nStates_loop(obj, n)
            % Propagates the state of the simulation through ntimesteps of length this.timestep

            % for n timesteps...
            for i = 2:n

                % if the particle jumps...
                if obj.hopOrNot() == 1

                    % randomly sample amongst 4 jump directions
                    d = obj.hopDirec();

                    % execute jump
                    switch d
                        
                        % Move (x,y) -> (x,y+a)
                        case 1
                            obj.x(i) = obj.x(i-1);
                            obj.y(i) = obj.y(i-1) + obj.lattConst;

                        % Move (x,y) -> (x+a,y)
                        case 2
                            obj.x(i) = obj.x(i-1) + obj.lattConst;
                            obj.y(i) = obj.y(i-1);

                        % Move (x,y) -> (x,y-a)
                        case 3
                            obj.x(i) = obj.x(i-1);
                            obj.y(i) = obj.y(i-1) - obj.lattConst;

                        % Move (x,y) -> (x-a,y)
                        otherwise
                            obj.x(i) = obj.x(i-1) - obj.lattConst;
                            obj.y(i) = obj.y(i-1);
                        
                    end

                % if it doesn't jump...
                else
                    
                    % record that its position is unaltered
                    obj.x(i) = obj.x(i-1);
                    obj.y(i) = obj.y(i-1);
                    
                end

            	% increment time by one timestep
            	obj.t(i) = obj.t(i-1) + obj.timestep;
            end
        end
        
    end
    
end