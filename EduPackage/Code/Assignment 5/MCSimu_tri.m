classdef MCSimu_tri < MCSimu
    
    methods(Static)
        
        function d = hopDirec()
            % Diff directions for different integers
            d = randi(6);
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
                        case 1
                            obj.x(i) = obj.x(i-1) + obj.lattConst;
                            obj.y(i) = obj.y(i-1);
                        case 2
                            obj.x(i) = obj.x(i-1) + 0.5*obj.lattConst;
                            obj.y(i) = obj.y(i-1) - 0.5*sqrt(3)* obj.lattConst;
                        case 3
                            obj.x(i) = obj.x(i-1) - 0.5*obj.lattConst;
                            obj.y(i) = obj.y(i-1) - 0.5*sqrt(3)* obj.lattConst;
                        case 4
                            obj.x(i) = obj.x(i-1) - obj.lattConst;
                            obj.y(i) = obj.y(i-1);
                        case 5
                            obj.x(i) = obj.x(i-1) - 0.5*obj.lattConst;
                            obj.y(i) = obj.y(i-1) + 0.5*sqrt(3)* obj.lattConst;
                        otherwise
                            obj.x(i) = obj.x(i-1) + 0.5*obj.lattConst;
                            obj.y(i) = obj.y(i-1) + 0.5*sqrt(3)* obj.lattConst;
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