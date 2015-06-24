classdef List
    %UNTITLED5 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        head
        at_end
    end
    
    methods
        function [b] = is_empty(obj)
            b = isempty(obj.head);
        end
        function [obj1] = append(obj, item)
            obj1 = obj;
            obj1.head = [obj1.head, item];
        end
        function [obj1] = remove(obj)
            obj1 = obj;
            obj1.head = [];
        end
    end
    
end

