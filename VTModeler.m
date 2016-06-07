classdef VTModeler < handle
    %VTModeler models the interaction of Ved and Ted loads
    %   It uses the CrossSection object to calculate the amount of shear
    %   and torsional reinforcement and graphs the dependence of
    %   reinforcement area to degree of utilization of shear and torsional
    %   capacity of the cross section
    
    properties
        cs;         % CrossSection object
        d;          % priblizna stat. visina presjeka d=0.9h
        z;          % priblizna visina resetke
        bw;         % sirina rebra
        Asw;        % povrsina jednog kraka uzengija
        A;          % povrsina rebra
        u;          % obim rebra
        tef;        % debljina efektivnog tankozidnog nosaca
        Ak;         % povrsina Ak
        uk;         % obim povrsine Ak
        q;          % tok smicanja
        v1;         % koef. redukcije cvrstoce na pritisak u betonu
        Asw_max;    % maksimalna kolicina poprecne armature
        Asw_min;    % minimalna kolicina poprecne armature
        s1=0;       % razmak konturnih uzengija
        s2=0;       % razmak unutrasnjih uzengija
        s_max;      % maksimalan razmak izmedju uzengija
        m;          % sjecnost
        ax;         % axes object for plotting
        
        cotTheta;   % cot. ugla theta
        Asw_sum;    % ukupna kolicina uzengija
        Asw_req;    % potrebna kolicina poprecne armature
        Asl;        % potrebna kolicina poduzne armature
        al;         % duzina prepustanja poduzne armature za dFtd
        VTcombined; % kombinovan uticaj V i T
        Vrd;        % maksimalna sila V koju presjek moze da prenese
        Trd;        % maks. momtent torzije koji presjek moze da prenese
    end

    % set i get metode
    methods
        % funkcija koja racuna Vrd_max za dato cotTheta
        function Vrd = calcVrd(this, cotTheta)
            cs = this.cs;
            bw = cs.dims.bw;
            Vrd = bw*this.z*this.v1*cs.fcd./(cotTheta + 1./cotTheta);
        end
        function Vrd = get.Vrd(this)
            Vrd = this.calcVrd(this.cotTheta);
        end
        % funkcija koja racuna Trd_max za dato cotTheta
        function Trd = calcTrd(this, cotTheta);
            cs = this.cs;
            Trd = this.v1*cs.fcd*2*this.Ak*this.tef./(cotTheta + 1./cotTheta);
        end
        function Trd = get.Trd(this)
            Trd = this.calcTrd(this.cotTheta);
        end
        
        % funkcija koja racuna izraz VTcombined = Ved/Vrd + Ted/Trd
        function VTcombined = calcVT(this, cotTheta)
            cs = this.cs;
            VTcombined = cs.Ved./this.calcVrd(cotTheta) + cs.Ted./this.calcTrd(cotTheta);
        end
        function VTcombined = get.VTcombined(this)
            VTcombined = this.calcVT(this.cotTheta);
        end
            
        % vektorizovana funkcija za proracun potrebne poprecne armature
        function [Asw_req, Asw_Ted, Asw_Ved] = calcAsw(this, cotTheta)
            if cotTheta == 0
                Asw_req = 0;
                Asw_Ted = 0;
                Asw_Ved = 0;
                return;
            end
            cs = this.cs;
            % potrebna kolicina poprecne armature za Ved po metru duznom
            % grede
            Asw_Ved = cs.Ved./(cs.fywd*this.z*cotTheta);
            % potrebna kolicina poprecne armature za Ted po metru duznom
            % grede - izraz se mnozi sa 2 jer je formula izvedena za jednu stranicu,
            % tj. jedan krak uzengija
            Asw_Ted = 2*cs.Ted./(cs.fywd*2*this.Ak*cotTheta);
            Asw_req = Asw_Ved + Asw_Ted;
        end
        
        function set.s1(this, s1)
            this.s1 = s1;
            % ako je sjecnost veca od 2, setuj s2 = s1
            if this.m > 2
                this.s2 = s1;
            else
                this.s2 = 0;
            end
            cs = this.cs;   % CrossSeciton - objekat poprecnog presjeka
            s2 = this.s2;   % razmak unutrasnjih uzengija
            Asw = this.Asw; % povrsina jednog kraka uzengije
            m = cs.m;       % sjecnost
            if m>2
                this.Asw_sum = 2*Asw/s1 + (m-2)*Asw/s2;
            else
                this.Asw_sum = 2*Asw/s1;
            end
            
            if this.Asw_sum > this.Asw_max
                message = ['Presjek je prearmiran, Asw/s(max) = ' num2str(round(this.Asw_max, 2)) '. '...
                    'Potrebno je smanjiti sjecnost ili usvojiti veci razmak uzengija.'];
                msgbox(message, 'Upozorenje', 'help');
            elseif s1 > this.s_max
                message = ['Usvojeni razmak s1 je veci od maks. dopustenog (s(max) = ' ...
                    num2str(round(this.s_max, 2)) ' mm). '...
                    'Potrebno je usvojiti manji razmak uzengija.'];
                msgbox(message, 'Upozorenje', 'help');
            end
        end
        function set.s2(this, s2)
            if s2 < this.s1
                return;
            end
            this.s2 = s2;
            cs = this.cs;   % CrossSeciton - objekat poprecnog presjeka
            s1 = this.s1;   % razmak unutrasnjih uzengija
            Asw = this.Asw; % povrsina jednog kraka uzengije
            m = cs.m;       % sjecnost
            if m>2
                this.Asw_sum = 2*Asw/s1 + (m-2)*Asw/s2;
            else
                this.Asw_sum = 2*Asw/s1;
            end
            % provjera da li unutrasnje uzengije sijeku prit. dijagonale
            cotTheta = this.cotTheta;
            if s2>this.z*cotTheta
                message = ['Razmak s2 je veci od horizontalnog raspona '...
                    'pritisnutih dijagonala (z*ctg(theta) = ' num2str(round(this.z*cotTheta, 2)) ' mm). '...
                    'Potrebno je smanjiti sjecnost ili usvojiti manji razmak uzengija.'];
                msgbox(message, 'Upozorenje', 'help');
            end
        end
        % racuna vrijednost cotTheta za koju je zbir interakcije VTcombined=1
        function cotTheta = calcCotTheta(this)
            cs = this.cs;
            % kada se u jednacinu Ted/Trd+Ved/Vrd=1 uvrste izrazi za Trd i
            % Vrd dobija se kvadratna jednacina u funkciji od cotTheta koja
            % se moze jednoznacno rijesiti
            k = this.v1*cs.fcd*2*this.Ak*this.tef*this.z*this.bw/...
                (this.bw*this.z*cs.Ted + 2*this.Ak*this.tef*cs.Ved);
            % odgovarajuce rjsesenje kvadratne jednacine
            cotTheta = (k+sqrt(k^2-4))/2;
            % provjera da li je dobijeno realno rjesenje jednacine
            if ~isreal(cotTheta)
                % ako nije, dodjeljuje se vrijednost 0
                cotTheta = 0;
            end
            % iz rjesenja se vidi da se minimalna realna vrijednost cotTheta moze
            % dobiti cotTheta = 1, pa nije potrebno vrsiti provjeru za
            % vrijednosti manje od 1
            % ukoliko se dobije cotTheta>2.5, usvaja se cotTheta=2.5
            if cotTheta > 2.5
                cotTheta = 2.5;
            end
        end
        
        function cotTheta = get.cotTheta(this)
            cs = this.cs;
            cotTheta = (cs.Ved/cs.fywd/this.z + cs.Ted/cs.fywd/this.Ak)/this.Asw_sum;
        end
        
        function [s1 s2] = calcS(this, cotTheta)
            % PRORACUN POTREBNE POPRECNE ARMATURE
            cs = this.cs;
            m = cs.m; % zadana sjecnost poprecne armature za Ved
            s1 = 0;
            s2 = 0;
            Asw = this.Asw;
            
            % ako je cotTheta=0 znaci da nema rjesenja, presjek je preslab
            if cotTheta == 0
                return;
            end
            % potrebna kolicina armature 
            [Asw_req Asw_Ted] = this.calcAsw(cotTheta);

            % maksimalna kolicina armature iz uslova cotTheta = 1
            this.Asw_max = this.calcAsw(1);
            
            % maksimalan dopusten razmak izmedju konturnih uzengija
            bw = this.bw; % sirina rebra
            % maksimalan dopusten razmak torzionih uzengija
            st = 2*Asw/Asw_Ted;
            this.s_max = min([0.75*this.d this.u/8 bw st]); % u/8 bw
            
            
            % potreban razmak uzengija za maks. sjecnost
            
            s = m*Asw/Asw_req;                

            if s > this.s_max
                s1 = this.s_max;
                % ako je s1 < s odredjuje se prvo kolicina torzione
                % (konturne) armature Asw_out koja se poredi sa ukupnom Asw_sum
                Asw_out = 2*Asw/s1;
                % ukoliko je m>2, dobice se potencijalno i odredjena kolicina unutrasnje armature
                % Asw_in > 0
                Asw_in = Asw_req - Asw_out;
            
                if Asw_in > 0
                    s2 = (m-2)*Asw/Asw_in;
                else
                    s2 = 0;
                    %Asw_req = Asw_out;
                    Asw_in = 0;
                end
            else
                s1 = s;
                if m>2
                    s2 = s;
                else
                    s2 = 0;
                end
            end
            % provjera da li je razmak unutrasnje armature veci od horizontalnog raspona 
            % kosih dijagonala, u tom slucaju izbaciti upozorenje da je presjek neracionalan, 
            % treba smanjiti sjecnost uzengija ili koristiti neki prikladniji tankozidni presjek.
            if s2 > this.z*cotTheta
                message = ['Nije moguce postici racionalan raspored unutrasnje armature '...
                    'za m>2 (s2 > z cot(Theta)), potrebno je smanjiti sjecnost ili usvojiti racionalniji presjek.'];
                msgbox(message, 'Upozorenje', 'help');
                return;
            end
        end
        
        function s_max = get.s_max(this)
            % maksimalan dopusten razmak izmedju konturnih uzengija
            s_max = min([0.75*this.d this.u/8 this.bw]); % 0.75d u/8 bw
        end
        
        function [Asl, al] = calcAsl(this, cotTheta)
            % dodatna poduzna armatura za torziju
            Asl = this.q*this.uk*cotTheta/this.cs.fywd;
            % duzina prepustanja glavne poduzne armature radi prihvatanja V
            % sila
            al = this.z*(cotTheta)/2;
        end
        
        function Asl = get.Asl(this)
            cotTheta = this.cotTheta;
            Asl = this.calcAsl(cotTheta);
        end
        
        function al = get.al(this)
            cotTheta = this.cotTheta;
            [~, al] = this.calcAsl(cotTheta);
        end
        
        function plotVT(this)
            ax = this.ax;
            % potrebno cotTheta s obzirom na opterecenje
            cotTheta = this.calcCotTheta();
            % generise 100 vrijednosti cotTheta izmedju 1 i 2.5, opadajuci
            cotThetaArray = linspace(2.5, 1, 100);
            Asw_sum_array = this.calcAsw(cotThetaArray);
            VTcombined_array = this.calcVT(cotThetaArray);
            
            % clear axes
            % plotanje glavnog dijagrama
            cla(ax);
            lineColor = 1*[0 0.4471 0.7412];
            ax(1).XGrid = 'on';
            ax(1).YGrid = 'on';
            ax(1).XLabel.String = 'A_{sw}/s [mm]';
            ax(1).YLabel.String = 'T_{Ed}/ T_{Rd} + V_{Ed}/ V_{Rd}';
            line('Color', lineColor,... %0.4*[1 1 1] [0.1529 0.2275 0.3725]
                'Tag', 'interaction_line',...
                'LineWidth', 2,...
                'Parent',ax(1),...
                'XData', Asw_sum_array,...
                'YData', VTcombined_array,...
                'Marker', '*',...
                'MarkerEdgeColor', 'none');
            % plotanje crvene tacke koja oznacava Ved/Vrd+Ted/Trd=1
%             line('XData', Asw_req*[1 1], 'YData', [1 1],...
%                 'Color', 'red', 'Parent', ax(1),...
%                 'Marker', 'o', 'MarkerFaceColor', 'red',...
%                 'MarkerEdgeColor', 'none');
           
            % plotanje crvene tacke koja oznacava Asw_req
            Asw_req = this.calcAsw(cotTheta);
            if Asw_req ~=0
                line('XData', Asw_req*[1 1], 'YData', min([VTcombined_array(1) 1])*[1 1],...
                    'Color', 'red', 'Parent', ax(1),...
                    'Marker', 'o', 'MarkerFaceColor', 'red',...
                    'MarkerEdgeColor', 'none');
            end
            % plotanje vertikalne crvene linije za min Asw/s
            % ima smisla samo ako je min Asw veca od kolicine za cotTheta
            % = 2.5
            if this.Asw_min > Asw_sum_array(1)
                yLim = ax(1).YLim;
                line('XData', this.Asw_min*[1 1], 'YData', yLim,...
                    'Color', 'red', 'Parent', ax(1));
                ax(1).YLim = yLim;
            end
        end
        
        function plotAsw(this)
            ax = this.ax;
            % potrebno cotTheta s obzirom na opterecenje
            cotTheta = this.calcCotTheta();
            % plotanje zelene tacke koja oznacava Asw_sum
            Asw_sum = this.Asw_sum;
            delete(findobj(ax.Children, 'flat', 'Tag', 'Asw_sum'));
            if Asw_sum <= this.Asw_max
                line('XData', Asw_sum*[1 1], 'YData', this.VTcombined*[1 1],...
                    'Color', 'green', 'Parent', ax(1),...
                    'Marker', 'o', 'MarkerFaceColor', 'green',...
                    'MarkerEdgeColor', 'none', 'Tag', 'Asw_sum');
            end
        end
    end
    methods
        % Constructor
        function this = VTModeler(cs, ax)
            this.cs = cs; % cs - CrossSection object
            this.ax = ax; % ax - axes object za plotanje dijagrama Asw - VTcombined
            this.d = 0.9*cs.dims.h; % staticka visina           
            this.z = 0.9*this.d; %this.z;
            % dimenzije rebra
            bw = cs.dims.bw;
            this.bw = bw;
            h = cs.dims.h;
            % povrsina rebra:
            this.A = bw*h;
            % obim rebra
            this.u = 2*bw + 2*h;
            % efektivna debljina zida min 2*c_nom+dsw, 
            this.tef = max([this.A/this.u, 2*cs.c_nom+cs.stirrup]);
            % povrsina Ak
            this.Ak = (bw-this.tef)*(h-this.tef);
            % obim povrsine Ak
            this.uk = 2*(bw-this.tef+h-this.tef);
            % tok smicanja
            this.q = cs.Ted/2/this.Ak;
            % koef. redukcije cvrstoce na pritisak u betonu
            this.v1 = 0.6*(1-cs.fck/250);
            % minimalna kolicina poprecne armature
            Rho_w_min = 0.08*sqrt(cs.fck)/cs.fywk;
            this.Asw_min = Rho_w_min*bw;
            this.Asw = cs.stirrup^2*pi/4; % povrsina jednog kraka uzengija
            % sjecnost
            this.m = cs.m;

            % provjeriti da li je postojeci presjek dovoljan
            % provjerava da li je kombinacija uticaja za theta=45 veca od 1
            cotTheta = this.calcCotTheta();
            if cotTheta == 0
                message = ['Potreban je jaci presjek. Ved/Vrd + Ted/Trd > 1'];
                msgbox(message, 'Upozorenje', 'help');
                return;
            end
        end
    end
    
end
