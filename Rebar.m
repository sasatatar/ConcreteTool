classdef Rebar < dynamicprops
    %REBAR Šipka armature
    %   Definiše pre?nik (Diameter [mm]) i površinu poprecnog presjeka (Area [cm^2])
    %   sipke, te x i y koordinata u poprecnom presjeku grede.
    
    properties
        d;              % precnik sipke [mm]
        Position;       % [x y] koordinate (Position)
        section;        % CrossSection handle
        %fyk = 500;      % 500 MPa za armaturu B500 
        Es;             % modul elasticnosti celika 200000 [MPa]
        %gamma_s = 1.15;  % koef. sigurnosti za fyd = fyk/gamma_s
        fyd;            % racunski dopusteni napon u celiku [MPa]
        x;              % x koordinata [mm]
        y;              % y koordinata [mm]
        row;            % red u kom se nalazi sipka
        column;         % kolona u kojoj se nalazi sipka
        zone;           % zona u kojoj sipka stoji (1 ili 2)
        Fs;             % sila u sipki
        euk = 50/1000;  % stvarna vrijednost epsilon uk
        k = 1.08;
    end
    properties (SetAccess = private)
        Area % povrsina poprecnog presjeka [mm^2]
    end
    properties (Dependent)
        strainHardening;
        eud;
    end
    
    methods
        %% constructor
        function obj = Rebar(cs,d,x,y,row,column,zone)
            if nargin == 7
                correction = (cs.ds_max(zone)-d)/2;
                if cs.ds_max(zone) ~= d && column == 1
                    y = y - correction;
                elseif cs.ds_max(zone) ~= d && column == cs.RPR(cs.ds_max(zone), zone)
                    y = y + correction;
                end
                if row == 1 && d ~= cs.ds_max(zone)
                    if zone == 1
                        x = x + correction;
                    else
                        x = x - correction;
                    end
                end
                
                obj.d = d;
                obj.section = cs;
                obj.x = x;
                obj.y = y;
                obj.row = row;
                obj.column = column;
                obj.zone = zone;
                % ovi podaci se vuku iz presjeka
%                 obj.Es = Es;
%                 obj.fyk = fyk;
%                 obj.gamma_s = gammas;
            elseif nargin == 1
                obj.section = cs;
            end
        end
        function out = get.strainHardening(this)
            out = this.section.strainHardening;
        end
        
        % dopustene dilatacije u celiku
        function eud = get.eud(this)
            eud = 0.9*this.euk;
        end
        
        %% Polozaj sipke (x i y koordinata u mm)
        function p = get.Position(this)
            p = [this.x this.y];
        end
        
        function set.Position(this,p)
            this.x = p(1);
            this.y = p(2);
        end
        
        %% povrsina poprecnog presjeka sipke
        function area = get.Area(this)
            area = this.Area;
        end
        
        %% precnik sipke [mm] (setter)
        function set.d(obj, d)
            obj.d = d;
            obj.Area = d^2*pi/4; % [mm^2]
        end
        
        %% racunski dopusteni napon celika, na granici tecenja
        function fyd = get.fyd(this)
            fyd = this.section.fyd;
        end
        
        function Es = get.Es(this)
            Es = this.section.Es;
        end
        
        %% napon u armaturi
        function sigmas = sigmas(this, x)
            s = this.section;
            if nargin < 2
                x = this.x;
            end
            % ide apsolutna vrijednost zbog negativnih vrijednosti dilatacija u
            % pritisnutoj zoni
            es1 = abs(s.strain(x)); 
            Es = this.Es;
            fyd = this.fyd;
            % granica plasticnosti
            epl = fyd / Es;
            if es1 >= epl
                if this.strainHardening == 0
                    sigmas = fyd;
                else
                    sigmas = fyd + (this.k*fyd-fyd)/(this.euk-epl)*(es1-epl);
                end
            else
                sigmas = es1*Es;
            end
        end
        
        %% sila Fs u sipki
        function Fs = get.Fs(this)
            Fs = this.Area * this.sigmas; %[N]
        end
        
    end
    % metode za crtanje sipke
    methods
        function rect = draw(this, ax)
            cs = this.section;
            x = this.x;
            y = this.y;
            d = this.d;
            row = this.row;
            column = this.column;
            zone = this.zone;
            rectX = x-d/2;
            rectY = y-d/2;
            
            rect = rectangle('Parent', ax,...
                'Position', [rectX rectY d d],...
                'Curvature', [1 1],...
                'LineStyle', 'none',...
                'FaceColor', [0 0.6 0.6],... %[0.2 0 1]
                'UserData', [row column zone]);
        end
        
        function showTag(this, ax, mousePosition)
            mouseX = mousePosition(1,1);
            mouseY = mousePosition(1,2);
            x = mouseX;
            y = mouseY;
            tag = text('String', [num2str(this.d) ' mm'],...
                'Units', 'data',...
                'Position', [x y 0],...
                'Tag', 'rebarTag',...
                'Parent', ax,...
                'FontSize', 11,...
                'BackgroundColor', [0.8706 0.4902 0],...
                'EdgeColor', 'none',...
                'Color', [1 1 1]);
            % mijenjam units u pixels da bih mogao definisati apsolutan polozaj taga
            % nezavisno od razmjere grafika
            tag.Units = 'pixels';
            % pomjera tag udesno za 5 px i nagore za 20 px u odnosu na
            % polozaj misa
            tag.Position(1:2) = [tag.Position(1)+5 tag.Position(2)+20];
        end
        
        function out = mouseOver(this, mousePosition)
        %%% MOUSEOVER provjerava da li je mis iznad sipke armature
            mouseX = mousePosition(1,1);
            mouseY = mousePosition(1,2);
            x = this.x-mouseX;
            y = this.y-mouseY;
            d = this.d;
            out = (x^2+y^2<=max([(d/2)^2 15^2]));
        end
    end
end

