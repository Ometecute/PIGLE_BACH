classdef MCSimu_hex < MCSimu
    
    methods(Static)
        
        function d = hopDirec()
            % Diff directions for different integers
            d = randi(3);
        end
        
    end
    
    
    
    methods    
        
        function nStates_loop(obj, n)
            
            % defined to be 0/1 to represent the type of site we are on
            siteType = 0;

            for i = 2:n
                if obj.hopOrNot() == 1
                    d = obj.hopDirec();
                    switch d
                        case 1
                            obj.x(i) = obj.x(i-1) + (-1)^siteType*obj.lattConst;
                            obj.y(i) = obj.y(i-1);
                        case 2
                            obj.x(i) = obj.x(i-1) - (-1)^siteType*0.5*obj.lattConst;
                            obj.y(i) = obj.y(i-1) - 0.5*sqrt(3)* obj.lattConst;
                        otherwise
                            obj.x(i) = obj.x(i-1) - (-1)^siteType*0.5*obj.lattConst;
                            obj.y(i) = obj.y(i-1) + 0.5*sqrt(3)* obj.lattConst;
                    end

                    % switch siteType if adsorbate moves, otherwise leave it alone
                    if siteType == 0
                        siteType = 1;
                    else
                        siteType = 0;
                    end
                    
                else
                    
                    % record that its position is unaltered
                    obj.y(i) = obj.y(i-1);
                    obj.x(i) = obj.x(i-1);
                    
                end
                
                % increment time by one timestep
                obj.t(i) = obj.t(i-1) + obj.timestep;
            end
            
        end
        
    end
    
end