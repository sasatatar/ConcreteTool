classdef Point
    %POINT Tacka u poprecnom presjeku 
    %   Tacka predstavljena x i y koordinatom sa koordinatnim pocetkom u
    %   gornjem lijevom uglu poprecnog presjeka
    
    properties
        x           % x koordinata [mm]
        y           % y koordinata [mm]
    end
    
    methods (Access = 'public')
        function obj = Point(x,y)
            if nargin == 0
                x = 0;
                y = 0;
            elseif nargin == 1
                y = x(2);
                x = x(1);
            end
            obj.x = x;
            obj.y = y;
        end
        
        function dx = dx(a,b) 
            dx = b.x - a.x;
        end
        
        function dy = dy(a,b)
            dy = b.y - a.y;
        end
        
        function dist = dist(a,b)
            dist = sqrt(dx(a,b)^2 + dy(a,b)^2);
        end        
    end
    methods (Static)
        function points = defineSection(p)
            [m n] = size(p);
            if n ~= 4
                error('Nepravilan unos: matrica mora imati 4 kolone.');
            end
            for i = m:-1:1
                points(i,2) = Point(p(i,3),p(i,4));
                points(i,1) = Point(p(i,1),p(i,2));
            end
        end
    end
    
end

