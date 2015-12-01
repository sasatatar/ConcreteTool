classdef CrossSection < handle
    %CROSSSECTION klasa predstavlja karakteristike poprecnog presjeka
    %nosaca
    %   Sve duzine su u [mm]
    properties 
        dims;               % dimenzije nosaca (struct)   
        %   .bf;           % sirina flanse
        %   .hf;           % visina flanse
        %   .hv;           % visina vute (prelaz imzedju flanse i rebra)
        %   .bw;            % sirina rebra
        %   .h;             % ukupna visina nosaca
        
        x;            % polozaj neutralne ose od vrha presjeka [mm]
        fck = 30;           % karakt. cvrtstoca betona [MPa] (C30/37)
        fctm = 2.9;         % cvrstoca betona na zatezanje (potrebno kod proracuna As1,min)
        alpha = 0.85;       % koef. koji uzima u obzir dugorocne negativne 
                            % faktore na cvrstocu betona
        gammac = 1.5;       % koeficijent sigurnosti za beton
        
        % podaci vezano za armaturu
        c_nom = 35;             % zastitni sloj betona [mm]
        stirrup = 12;           % uzengije [mm]
        Rebars = Rebar.empty;   % matrica sipki armature
        rebarMesh;              % matrica sa 1 na zauzetim mjestima za sipke
        dg = 32;                % najvece zrno agregata  
        ds_max = [19 19];       % maksimalan precnik armature [mm] posebno za zategnutu (1) 
                                % posebno za pritisnutu zonu (2)
        As1_req = 0;            % potrebna povrsina armature [mm^2]
        As2_req = 0;
        fyk = 500;              % karakteristicna cvsrtoca celika [MPa]
        Es = 200000;            % modul elasticnosti celika [MPa]
        gamma_s = 1.15;         % koef. sigurnosti za fyd = fyk / gamma_s
        fyd;                    % racunski dopusteni napon u celiku 
        minRebarSpacing;        % minimalni svijetli razmak izmedju sipki
        xdRatio;                % odnos x/d = 0.45 za C<=50/60, 0.35 za C>=55/67
        Nsd = 0;                % normalna sila koja djeluje u presjeku (pritisak je pozitivan)
    end
    %% Dependent
    properties (Dependent)
        % tacke koje definisu poprecni presjek u sledecem formatu:
        % [x11 y11 x12 y12 % x11,y11 koordinate prve konturne tacke
        %  x21 y21 x22 y22...]
        Points;
        fcd; % design concrete strength fcd = alfacc * fck / gammac
        ecu2;       % strain in concrete / dilatacija u betonu
    end
    %% Read-only properties - pristupa im se preko get metoda
    properties (SetAccess = 'private')
        Fc;                 % sila u betonu Fc [N]
        xFc;                % polozaj sile Fc - x koordinata
        Fs1;                % sila u zategnutoj armaturi [N]
        xFs1;               % polozaj sile Fs1 - x koordinata
        Fs2;                % sila u pritisnutoj armaturi [N]
        xFs2;               % polozaj sile Fs2 - x koordinata
        As1;                % ukupna povrsina ugradjene zat. armature [mm^2]
        As2;                % ukupna povrsina ugradjene pritisnute arm. [mm2]
        Ac;                 % povrsina betona
        Mrd;                % reaktivni moment savijanja (kNm)
        ignoreCRebars = 0;  % boolean flag indicating if compressed rebars should 
                                    % be ignored when calculating Fc
    end%%
    
    %% Metode za armaturu
    methods
        function As1 = get.As1(this)
            rebars = findobj(this.Rebars, 'zone', 1);
            As1 = sum([rebars.Area]);
        end
        function As2 = get.As2(this)
            rebars = findobj(this.Rebars, 'zone', 2);
            As2 = sum([rebars.Area]);
        end
        
        function fyd = get.fyd(this)
            fyd = this.fyk / this.gamma_s; %[MPa tj. N/mm2]
        end
        
        function rect = addRebar(this, ax, ds, mouseX, mouseY, zone)
            %%% ADDREBAR dodaje sipke armature u poprecni presjek, korisnik
            %%% klikom na poprecni presjek odredjuje u koju kolonu presjeka
            %%% da se postavi zadana sipka armature
            dims = this.dims;
            h = dims.h;
            c_nom = this.c_nom;
            stirrup = this.stirrup;
            ds_max = this.ds_max;
            % vertikalni osni razmak izmedju sipki
            rowSpacing = this.minRebarSpacing(zone) + ds_max(zone);
            % rpr - maksimalan broj sipki u redu
            % columnSpacing - horizontalni osni razmak izmedju sipki 
            [rpr, columnSpacing] = this.RPR(ds_max(zone), zone);
            
            % odredjivanje x koordinate sipke
            % x koordinata prvog reda sipki
            if zone == 1
                x0 = h - c_nom - stirrup - ds_max(zone)/2;
                x = x0 - mouseX;
            elseif zone == 2
                x0 = 0 + c_nom + stirrup + ds_max(zone)/2;
                x = mouseX - x0;
                % mnozi se sa -1 zbog minusa u formuli za x
            end
            % odabrani red sipki
            row = round(x/rowSpacing)+1;
            if zone == 1
                % tacna x koordinata sipke
                x = x0 - (row-1)*rowSpacing;
            elseif zone == 2
                % tacna x koordinata sipke
                x = x0 + (row-1)*rowSpacing;
            end
            
            % odredjivanje y koordinate sipke
            % y koordinata prve kolone sipki
            y0 = (dims.bf-dims.bw)/2+c_nom+stirrup+ds_max(zone)/2;
            y = mouseY - y0;
            % odabrana kolona
            column = round(y/columnSpacing)+1;
            % tacna y koordinata sipke
            y = y0 + (column-1)*columnSpacing;
            % provjera da li je klik u okviru presjeka:
            % broj dodatnih kolona u flansi
            fColumns = floor((y0-(c_nom+stirrup+ds_max(zone)/2))/columnSpacing)+1;
            if zone == 1 && (row<1 || column<1 || column>rpr)
                rect = [];
                return;
            elseif zone == 2 && (row<1 || column<1-fColumns || column>rpr+fColumns)
                rect = [];
                return;                
            end
            
            % ukoliko zadani polozaj vec nije zauzet, dodati sipku
            if isempty(findobj(ax, 'UserData', [row column zone]))
                rebar = Rebar(ds, this, x, y, row, column, zone, this.Es, this.fyk);
                this.Rebars(end+1) = rebar;
                % pozvati funkciju za crtanje sipke - dok sredim
                rect = rebar.draw(ax);
            else
                disp('mjesto je zauzeto');
                rect = [];
            end
        end
        % funkcija koja plota pomocnu mrezu za polozaj armature
        function plotGrid(this, ax, zone)
            grids = findobj(ax, 'Type', 'line', 'Tag', 'grid');
            delete(grids);
            dims = this.dims;
            h = dims.h;
            c_nom = this.c_nom;
            stirrup = this.stirrup;
            ds_max = this.ds_max;
            if zone == 1
                % x koordinata prvog reda sipki u zategnutoj zoni
                x0 = h - c_nom - stirrup - ds_max(zone)/2;
            elseif zone == 2
                % x koordinata prvog reda sipki u pritisnutoj zoni
                x0 = 0 + c_nom + stirrup + ds_max(zone)/2;
            end
            % y koordinata prve kolone sipki
            y0 = (dims.bf-dims.bw)/2+c_nom+stirrup+ds_max(zone)/2;
            % vertikalni osni razmak izmedju sipki
            rowSpacing = this.minRebarSpacing(zone) + ds_max(zone);
            % rpr - maksimalan broj sipki u redu
            % columnSpacing - horizontalni osni razmak izmedju sipki 
            [rpr, columnSpacing] = this.RPR(ds_max(zone), zone);
            % vertikalne ose
            x = [0 h];
            if zone == 1
                a = 1;
                b = rpr;
            elseif zone == 2
                % broj dodatnih kolona u flansi
                fColumns = floor((y0-(c_nom+stirrup+ds_max(zone)/2))/columnSpacing);
                % lijevo od pocetnog polozaja dodaju se negativne kolone u
                % sirini flanse
                a = 1-fColumns; 
                % desno od krajnjeg polozaja u rebru se dodaje isti broj
                % kolona, u sirini flanse
                b = rpr+fColumns; % desno se dodaju
            end
            for i=a:b
                % y koordinate vertikalnih osa u pritisnutoj zoni
                y = [y0 y0]+(i-1)*columnSpacing; 
                line('Parent', ax,...
                    'XData', x, 'YData', y, 'Tag', 'grid',...
                    'Color', 0.6*ones(1,3), 'HitTest', 'off',...
                    'LineStyle', ':');
            end
            % horizontalne ose
            if zone == 1
                a = 0;
                b = floor(0.65*h/rowSpacing);
                % y koordinate za horizontalne ose u zategnutoj zoni
                y = [(dims.bf-dims.bw)/2 (dims.bf-dims.bw)/2+dims.bw];
            elseif zone == 2
                a = -floor(0.5*h/rowSpacing);
                b = 0;
                % y koordinate za horizontalne ose u pritisnutoj zoni
                y = [0 dims.bf];
            end
            for i=a:b
                % x koordinate za horizontalne ose
                x = [x0 x0]-i*rowSpacing;
                line('Parent', ax,...
                    'XData', x, 'YData', y, 'Tag', 'grid',...
                    'Color', 0.6*ones(1,3), 'HitTest', 'off',...
                    'LineStyle', ':');
            end
        end
        
        function rspacing = rebarSpacing(this,ds)
            % staviti: Access = 'private'
            dg = this.dg; % najvece zrno agregata     
            spacing = zeros(3, length(ds)); 
            % uslovi za min. razmak izmedju sipki (prema EC2)
            spacing(1,:) = ds;   % precnik armature
            spacing(2,:) = dg+5; % precnik najveceg zrna agregata + 5mm
            spacing(3,:) = 20;   % 20mm
            rspacing = max(spacing); % minimalan svijetli razmak izmedju sipki
        end
        % pomocna funkcija koja odredjuje minimalan razmak izmedju sipki
        % armature u redu, prebaciti u private blok
        function rspacing = get.minRebarSpacing(this)
            rspacing = this.rebarSpacing(this.ds_max);
        end
        function [rpr, rspacing] = RPR(this,ds,zone) 
            %%% RPR(ds) (rebars per row) Odredjuje maksimalan broj sipki u
            %%% jednom redu i osni razmak izmedju njih u [mm] za zadane precnike ds.
            %%% ds - vektor precnika armature u mm
            dims = this.dims;
            bw = dims.bw; % sirina flanse / rebra u zategnutoj zoni [mm]
            dg = this.dg; % najvece zrno agregata
            c_nom = this.c_nom; % ukupan zastitni sloj betona
            stirrup = this.stirrup; % uzengije
            spacing = zeros(3, length(ds)); 
            rebar_wall_spacing = this.minRebarSpacing(zone); % minimalan svijetli razmak izmedju sipki
            % maksimalan broj sipki u jednom redu
            rpr = floor((bw - 2*(c_nom + stirrup) - ds)./(rebar_wall_spacing+ds))+1;
            % osni razmak izmedju sipki u jednom redu [mm]
            rspacing = (bw - 2*(c_nom + stirrup) - ds)./(rpr-1); 
        end
        
        function [As1_pot, As2_pot] = calculateAs(this, Msd, d)
            %%% CALCULATEAS racuna potrebnu kolicinu armature 
            %%% Msd - moment savijanja u presjeku u Nmm
            % ako d nije definisano, usvaja se 0.9h
            %this.ignoreCRebars = 1;
            dims = this.dims;
            if nargin <= 2
                d = 0.9*dims.h;
            end
            % povecanje momenta usljed djelovanja normalne sile
            Nsd = this.Nsd;
            Msd = Msd + Nsd*(d-dims.h/2);
            
            % pretpostavlja se jednostruko armiranje
            As1_pot = 0;
            As2_pot = 0;
            delta = 0.01; % zeljena preciznost (dopustena razlika izmedju Msd i Mrd)
            xmin = 0;
            % gornja granica za x je staticka visina presjeka (od
            % pritisnute ivice do najudaljenijeg reda armature
            xmax = this.xdRatio*(dims.h-this.c_nom-this.stirrup-this.ds_max(1)/2); 
            
            % maksimalni moment savijanja sa jednostrukim armiranjem (Mrd
            % za x = xmax)
            this.x = xmax; % maksimalna vrijednost x dozvoljena prema EC2
            Mrd = getMrd();
            Fc = this.Fc;
            
            % Minimalna i maksimalna kolicina armature:
            
            % minimalna kolicina armature prema EC2:
            As1_min = max([0.26*this.fctm/this.fyk*dims.bw*dims.h...
                0.0013*dims.bw*dims.h]); % ili d umjesto dims.h?
            % ukupna maskimalna kolicina armature (As1+As2)
            As_max = 0.04*this.Ac; % ili 0.04*dims.bw*dims.h
            
            % ako nije definisan Msd, usvaja se As1,min
            if Msd == 0
                As1_pot = As1_min;
                return;
            end

            % provjera da li je potrebno dvostruko armiranje           
            if Mrd < Msd
                % krak sila
                z = d-this.xFc;
                Fs1 = Mrd / z - Nsd; % [N]
                dM = Msd - Mrd;
                if dM > 0.5*Mrd
                    msgbox('Msd > 1.5 Mrd,max pri jednostrukom armiranju, potrebna je veca visina presjeka.',...
                    'Premalen presjek', 'help');
                    return;
                end
                % pretpostavka o polozaju tezista pritisnute armature
                % manja vrijednost od 0.1h i polozaja tezista jednog reda
                % armature precnika 36 mm (najveca moguca armatura)
                d2 = min([0.1*dims.h this.c_nom+this.stirrup+32/2]);
                zd = d - d2;
                Fs2 = dM / zd; % [N]
                As1_pot = (Fs1+Fs2) / this.fyd; % mm2
                As2_pot = Fs2 / this.fyd; % mm2
                return;
            end
            
            % trazenje rjesenja za jednostruko armiranje
            dM = Mrd-Msd;
            while dM > delta
                % ako je Mrd > Msd, neutralna osa se mora nalaziti iznad
                % trenutnog polozaja (dakle x se smanjuje)
                if Mrd > Msd
                    xmax = this.x;
                else % u suprotnom, x se povecava tako sto povecavamo donju granicu (xmin)
                    xmin = this.x;
                end
                this.x = (xmin+xmax)/2;
                Mrd = getMrd();
                dM = abs(Mrd-Msd);
            end
            % nakon postignute zadovoljavajuce preciznosti, racunamo
            % potrebnu povrsinu armature na osnovu poznate sile u betonu (Fc)
            As1_pot = (this.Fc-Nsd)/this.fyd; % [mm2]
            if (As1_pot+As2_pot)>As_max
                disp('As > As_max');
                msgbox('As,uk > As,max. Potrebno je usvojiti veci presjek.',...
                    'Prearmiran presjek', 'help');
            end
            
            % pomocna funkcija koja racuna pribliznu vrijednost reaktivnog
            % momenta Mrd na osnovu poznate sile u betonu (Fc) i kraka
            % unutrasnjih sila (z)
            function Mrd = getMrd()
                Fc = this.Fc;
                xFc = this.xFc;
                z = d-xFc;
                Mrd = Fc*z; % [Nmm]
            end
        end
    end
    
    
    %% SET i GET metode
    methods       
        function xd = get.xdRatio(this)
            %%% odnos x/d = 0.45 za C<=50/60, 0.35 za C>=55/67
            if this.fck <= 50
                xd = 0.45;
            else
                xd = 0.35;
            end
        end
        
        function set.x(this, x)
            %%% SET.X zadaje polozaj neutralne use i poziva funciju za
            %%% azuriranje grafikona (treba sve funkcije za plotanje
            %%% staviti u jednu glavnu funkciju)
            this.x = x;
            %this.plotCompression;
            %this.plotStrain;
        end
        
        function p = get.Points(this)
            %%% POINTS generise matricu tacaka koje definisu poprecni
            %%% presjeka nosaca na osnovu zadanih dimenzija (tacke su
            %%% potrebne za graficki prikaz nosaca). Izlaz je u sledecem
            %%% formatu: % [x11 y11 x12 y12; x21 y21 x22 y22...] x11,y11 koordinate prve
            %%% konturne tacke [mm] na lijevoj ivici, x12, y12 koordinate
            %%% odgovarajuce tacke na desnoj ivici presjeka  
            dim = this.dims;
            p = zeros(4,4);
            p(2) = dim.hf;
            p(3,1:2) = [dim.hf+dim.hv (dim.bf-dim.bw)/2];
            p(4,1:2) = [dim.h (dim.bf-dim.bw)/2];
            p(:,3) = p(:,1);
            p(:,4) = p(:,2) + [dim.bf dim.bf dim.bw dim.bw]';
        end
        
        function fcd = get.fcd(this)
            fcd = this.alpha * this.fck / this.gammac;
        end
        
        function ecu2 = get.ecu2(this)
            % maksimalna dopustena dilatacija u betonu
            % usvaja se da je negativna
            if this.fck < 50
                ecu2 = -3.5/1000;
            else
                % formula iz EC2
                ecu2 = -(2.6+35*((90-this.fck)/100)^4)/1000;
            end
        end
        
        function Fc = get.Fc(this)
            %%% FC racuna intenzitet rezultantne sile u betonu Fc [N]
            x = this.x;
            Fc = abs(integral(@this.sigmacb, 0, x)); % N
%             rebars = findobj(this.Rebars, 'zone', 2);
%             if ~isempty(rebars) %&& ~this.ignoreCRebars
%                 % umanjuje silu Fc za dio koji je nosila povrsina koju sad
%                 % zauzima armatura u pritisnutoj zoni
%                 Fc = Fc - sum([rebars.Area])*this.sigmac(this.xFs2);
%             end
        end
        
        function xFc = get.xFc(this)
            %%% XFC racuna polozaj tezista sile u betonu Fc racunajuci od
            %%% vrha presjeka [mm]
            Fc = this.Fc; % [N]
            xFc = abs(integral(@this.sigmac_moment, 0, this.x)/Fc);  %[mm]
%             rebars2 = findobj(this.Rebars, 'zone', 2);
%             if ~isempty(rebars2)
%                 dFc = sum([rebars2.Area])*this.sigmac(this.xFs2);
%                 xFc = ((Fc+dFc)*xFc - dFc*this.xFs2)/Fc; % provjeriti s prof.
%             end
        end
        
        function Fs2 = get.Fs2(this)
            %%% FS2 racuna rezultantnu silu u celiku u pritisnutoj zoni
            rebars = findobj(this.Rebars, 'zone', 2);
            Fs2 = -sum([rebars.Fs]); % da bi se dobila pozitivna vrijednost u N
        end
        
        function xFs2 = get.xFs2(this)
            %%% XFS2 racuna polozaj sile Fs2 na x osi mjereno od
            %%% vrha presjeka
            rebars = findobj(this.Rebars, 'zone', 2);
            if isempty(rebars)
                xFs2 = 0;
            else
                xFs2 = sum([rebars.Fs].*[rebars.x])/sum([rebars.Fs]);
            end
        end
        
        function Fs1 = get.Fs1(this)
            %%% GET.FS1 racuna rezultantnu silu u celiku u zategnutoj zoni
            %%% (Fs1) u [N]
            rebars = findobj(this.Rebars, 'zone', 1);
            Fs1 = sum([rebars.Fs]);
        end
        
        function xFs1 = get.xFs1(this)
            %%% XFS1 racuna polozaj sile Fs1 na x osi, mjereno od
            %%% vrha presjeka
            rebars = findobj(this.Rebars, 'zone', 1);
            if isempty(rebars)
                xFs1 = 0;
            else
                xFs1 = sum([rebars.Fs].*[rebars.x])/sum([rebars.Fs]);
            end
        end
        
        
        % get.Mrd funkcija racuna reaktivni moment savijanja presjeka, na
        % osnovu usvojene armature
        function Mrd = get.Mrd(this)
            if isempty(findobj(this.Rebars, 'zone', 1))
                Mrd = 0;
                return;
            end
            dims = this.dims;
            Nsd = this.Nsd;
            delta = 0.01; % zeljena preciznost (dopustena razlika izmedju Msd i Mrd)
            xmin = 0;
            xmax = this.dims.h;
            this.x = (xmax+xmin)/2;
            
            dF = abs(this.Fs1+Nsd-this.Fc-this.Fs2);
            while dF > delta
                % ako je Fc > Fs1, neutralna osa se mora nalaziti iznad
                % trenutnog polozaja (dakle x se smanjuje)
                if (this.Fc+this.Fs2) > (this.Fs1 + Nsd)
                    xmax = this.x;
                else % u suprotnom, x se povecava tako sto povecavamo donju granicu (xmin)
                    xmin = this.x;
                end
                this.x = (xmin+xmax)/2;
                dF = abs(this.Fs1+Nsd-this.Fc-this.Fs2);
            end
            % Moment savijanja u odnosu na teziste zategnute armature
            Mrd = this.Fc*(this.xFs1-this.xFc)+this.Fs2*(this.xFs1-this.xFs2)-Nsd*(this.xFs1-dims.h/2); % [Nmm]
        end
        
        function Ac = get.Ac(this)
            Ac = integral(@this.b, 0, this.dims.h); % [mm2]
        end
    end
    
    methods
        function this = CrossSection(d)
            %%% CROSSSECTION Constructor:
            %%% d - vektor koji sadrzi sledece parametre: 
            %%% bf; sirina flanse
            %%% hf; visina flanse
            %%% hv; visina vute
            %%% bw; sirina rebra
            %%% h; ukupna visina nosaca
            %%% U slucaju kada je nosac bez flanse, dovoljno je
            %%% definisati visinu i sirinu
            if nargin ~= 0
                n = length(d);
                if n >= 5
                    this.dims.bf = d(1);
                    this.dims.hf = d(2);
                    this.dims.hv = d(3);
                    this.dims.bw = d(4);
                    this.dims.h = d(5);
                end
                if n == 2
                    b = d(1);
                    h = d(2);
                    this.dims.bf = b;
                    this.dims.hf = 0;
                    this.dims.hv = 0;
                    this.dims.bw = b;
                    this.dims.h = h;
                end
            else % default dimenzije presjeka
                this.dims.bf = 1350;
                this.dims.hf = 100;
                this.dims.hv = 200;
                this.dims.bw = 400;
                this.dims.h = 1200;
            end
            this.x = this.dims.h*0.9*this.xdRatio;
        end
        
        function plotSection(this,ax)
            %%% PLOTSECTION plota presjek nosaca sa uzengijama
            % uklanjanje postojecih linija sa grafa
            %delete(findobj(ax.Children, 'flat', 'Tag', 'crossSection_line'));
            %delete(findobj(ax.Children, 'flat', 'Tag', 'stirrup_patch'));
            cla(ax);            
            ax.CameraUpVector = [-1 0 0]; % rotira koord. sistem za 90 stepeni u smjeru kazaljke
            ax.DataAspectRatio = [1 1 1];
            ax.XLim(2) = this.dims.h; % podesava opseg na x osi
            
            section_line = line('Color', [0.4 0.4 0.4],...
                'Tag', 'crossSection_line',...
                'LineWidth', 1.3,...
                'Parent',ax);
            % x i y koordinate tacaka koje se spajaju linijom
            x = [this.Points(:,1)' flip(this.Points(:,3)') this.Points(1)];
            y = [this.Points(:,2)' flip(this.Points(:,4)') this.Points(1,2)];
            %ax.YLim = [min(y)-50, max(y)+50];
            section_line.XData = x;
            section_line.YData = y;
            
            % uzengije
            c_nom = this.c_nom;
            dims = this.dims;
            stirrup = this.stirrup;
            xs = [0 dims.h dims.h 0 0];
            y0 = (dims.bf-dims.bw)/2;
            ys = [y0 y0 y0+dims.bw y0+dims.bw y0];
            % vektor transformacije niza da bi se offsetovao za vrrijednost
            % c_nom (gdje se dodaje c_nom, ide 1, a gdje se oduzima, -1)
            transx = [1 -1 -1 1 1];
            transy = [1 1 -1 -1 1];
            xs = xs+c_nom*transx;
            ys = ys+c_nom*transy;
            xs = [xs flip(xs)+stirrup*flip(transx)];
            ys = [ys flip(ys)+stirrup*flip(transy)]; 

            patch('FaceColor', [0 0.6 0.8],...
                'Tag', 'stirrup_patch',...
                'Parent', ax,...
                'LineStyle', 'none',...
                'XData', xs, 'YData', ys);
            
            % horizontalni krak uzengija uz gornju ivicu
            xs = [0 0];
            ys = [0 dims.bf];
            transx = [1 1];
            transy = [1 -1];
            xs = xs+c_nom*transx;
            ys = ys+c_nom*transy;
            xs = [xs flip(xs)+stirrup*flip(transx)];
            ys = [ys flip(ys)]; % +stirrup*flip(transy)
            patch('FaceColor', [0 0.6 0.8],...
                'Tag', 'stirrup_patch',...
                'Parent', ax,...
                'LineStyle', 'none',...
                'XData', xs, 'YData', ys);  
        end %%

        function plotCompression(this, ax)
            %%% PLOTCOMPRESSION Oznacava pritisnutu zonu poprecnog presjeka
            % brise postojece compressedZone graficke objekte (reset grafika)
            cZone = findobj(ax, 'Tag', 'compressedZone');
            delete(cZone);
            compressedZone = patch('FaceColor', 0.4*ones(1,3),...
                'Tag', 'compressedZone',...
                'Parent', ax,...
                'LineStyle', 'none',...
                'FaceAlpha', 0.2);
            % x i y koordinate tacaka koje se spajaju linijom
            x = [this.Points(:,1)' flip(this.Points(:,3)')];
            y = [this.Points(:,2)' flip(this.Points(:,4)')];
            % obiljezavanje pritisnute povrsine betona
            [bx,y_left] = this.b(this.x);
            
            % tacke A i B definisu neutralnu liniju 
            % provjeriti kriticni slucaj kada je neutralna usa u zoni donje
            % flanse (za slucaj da je nosac u obliku I presjeka)
            A = [this.x,y_left];
            B = [this.x,y_left+bx];
            logical_index = x<=this.x;
            xc = x(logical_index);
            yc = y(logical_index);
            middle_index = length(xc)/2;
            xc = [xc(1:middle_index) A(1) B(1) xc(middle_index+1:end)];
            yc = [yc(1:middle_index) A(2) B(2) yc(middle_index+1:end)];

            compressedZone.XData = xc;
            compressedZone.YData = yc;
        end %%
        function plotStrain(this, ax, section_axes)
            %%% PLOTSTRAIN generise grafik dilatacija
            
            % calculate x and y data for patch object
            zone = 1;
            d = this.dims.h - this.c_nom - this.stirrup - this.ds_max(1)/2;
            x = [0 d d 0];
            y = [this.ecu2 this.strain(d) 0 0];
            
            % odnos za uskladjivanje visine grafa dilatacija sa prikazom
            % poprecnog presjeka
            ratio = ax.Position(3)/section_axes.Position(3);
                        
            % setup axes properties
            cla(ax);            
            ax.CameraUpVector = [-1 0 0]; % rotira koord. sistem za 90 stepeni u smjeru kazaljke
            set(ax,'YDir','reverse',...
                'PlotBoxAspectRatio', section_axes.PlotBoxAspectRatio.*[1 ratio 1],...
                'XLim', section_axes.XLim);
            
            % provjeriti da li ima sipki u zategnutoj zoni, 
            % ako ima, provjeriti odnos x/d i promijeniti boju plota
            color = [0.6 0.8 0];
            if this.xFs1 ~= 0
                if this.x/d > this.xdRatio % umjesto xFs1 da ide koordinata zadnjeg reda armature?
                    color = [0.8 0 0];
                end
            end
            patch('FaceColor', color,...
                'Tag', 'strain_patch',...
                'Parent', ax,...
                'LineStyle', 'none',...
                'FaceAlpha', 0.6,...
                'XData', x, 'YData', y);
                
           label = text(0, this.ecu2, num2str(this.ecu2*1000, '%0.1f'),...
               'Parent', ax,...
               'Color', [0.3059 0.3961 0.5804]); % numericki prikaz dilatacije u betonu u promilima
           % mijenjam units u pixels da bih mogao definisati apsolutan polozaj taga
           % nezavisno od razmjere grafika
           label.Units = 'pixels';
           % pomjera tag udesno za 5 px i nadole za 10 px u odnosu na
           % prvobitni polozaj
           label.Position(1:2) = [label.Position(1)+5 label.Position(2)-10];
           
           % dilatacija u celiku u promilima
           text(this.dims.h-50, -1/1000, num2str(this.strain(this.dims.h)*1000, '%0.2f'),...
               'Parent', ax,...
               'Color', [0.3059 0.3961 0.5804]); 
           
           % oznacavanje x koordinata
           ax.XTick = [0 this.x d];
           ax.XTickLabel = sprintf('%.1f\n', ax.XTick);
        end
        
        function plotStress(this, ax, section_axes)
            % odnos za uskladjivanje visine grafa napona sa prikazom
            % poprecnog presjeka
            ratio = ax.Position(3)/section_axes.Position(3);
                        
            % setup axes properties
            cla(ax);            
            ax.CameraUpVector = [-1 0 0]; % rotira koord. sistem za 90 stepeni u smjeru kazaljke
            set(ax,'YDir','normal',...
                'PlotBoxAspectRatio', section_axes.PlotBoxAspectRatio.*[1 ratio 1],...
                'XLim', section_axes.XLim);
            numPoints = 100;
            x = linspace(0, this.x, numPoints);
            y = this.sigmac(x);
            % dodavanje koordinata pocetne tacke na kraju niza da se
            % zatvori kontura
            x = [x 0];
            y = [y 0];
            
            patch('FaceColor', 0.8627*[1 1 1],...
                'Tag', 'stress_patch',...
                'Parent', ax,...
                'LineStyle', 'none',...
                'XData', x, 'YData', y);
            
            % oznacavanje koordinata polozaja sila i n. ose na x osi
            ax.XTick = unique([this.xFc this.xFs1 this.xFs2 this.x]);
            ax.XTickLabel = sprintf('%.1f\n', ax.XTick);
            
            % oznacavanje fcd na y osi
            ax.YTick = [0 this.fcd];
            ax.YTickLabel = sprintf('%.1f\n', ax.YTick);
            
            % ucrtavanje sile u betonu
            fcolor = [0.2 0.2 0.2];
            x = this.xFc;
            drawVector(x, '<', fcolor, this.Fc);
            
            % ucrtavanje sile u zategnutoj armaturi
            if ~isempty(findobj(this.Rebars, 'zone', 1))
                fcolor = [0 0.6 0.8];
                x = this.xFs1;
                drawVector(x, '>', fcolor, this.Fs1);
            end
            
            % ucrtavanje sile u pritisnutoj armaturi
            if ~isempty(findobj(this.Rebars, 'zone', 2))
                fcolor = [0 0.6 0.8];
                x = this.xFs2;
                drawVector(x, '<', fcolor, this.Fs2);
            end
            
            % pomocna funkcija za ucrtavanje vektora sila
            function drawVector(x, triangleStyle, color, force)
                % crta silu
                fx = x;
                fy = 0.7*ax.YLim(2);
                f = line(fx*[1 1], [0 fy], 'Parent', ax, ...
                    'LineWidth', 1.5, 'Color', color);
                % dodaje strelicu na vrh linije
                ax.Units = 'points';
                px = fx;
                py = 0;
                p = line(px*[1 1], [py py], 'Parent', ax,...
                'Marker', triangleStyle, 'MarkerFaceColor', color,...
                'MarkerEdgeColor', 'none');
                dy = 0;
                if strcmp(triangleStyle, '>')
                    dy = 0.7*ax.YLim(2);
                end
                p.YData = dy + p.MarkerSize/ax.OuterPosition(3)*ax.YLim(2)*[1 1];
                ax.Units = 'pixels';
                
                % velicina sile [kN] - tekst
                label = text(x, fy, [num2str(force/1000, '%0.2f') ' kN'],...
                    'Parent', ax,...
                    'Color', color); % numericki prikaz dilatacije u betonu u promilima
                % mijenjam units u pixels da bih mogao definisati apsolutan polozaj taga
                % nezavisno od razmjere grafika
                label.Units = 'pixels';
                % pomjera tag udesno za 5 px i nadole za 10 px u odnosu na
                % prvobitni polozaj
                label.Position(1:2) = [label.Position(1)+10 label.Position(2)];
            end
            
            
        end
        
        
        % funkcija za crtanje sipki armature u poprecnom presjeku
        function drawRebar(this, ax)
            % obrisi posotjece sipke
            rebars = findobj(ax, 'Type', 'rectangle');
            delete(rebars);
            % iscrtaj ponovo sve sipke
            for i = 1:numel(this.Rebars)
                this.Rebars(i).draw(ax);
            end
        end
    end
    
    
    methods (Access = 'public')
        %% daje sirinu nosaca b na odstojanju x od vrha u mm
        function [b,y] = b(this,x) 
            % broj koraka za for petlju jednak je broju parova tacaka
            % umanjenom za jedan (broj redova -1)
            steps = size(this.Points,1)-1;
            b = zeros(1,length(x));
            y = b;
            for i = 1:steps
                % definise konturne tacke za segment "i" nosaca
                A = this.Points(i:i+1,:);
                % odredjuje koje tacke ulaznog niza pripadaju tom segmentu
                logical_index = x>=A(1) & x<=A(2);
                if sum(logical_index)==0
                    continue
                end
                % sirina nosaca izmedju prvog i drugog para tacaka za
                % segment "i"
                dy1 = A(1,4)-A(1,2);
                dy2 = A(2,4)-A(2,2);
                
                % ako sirine nisu jednake, sirina na odstojanju "x" se
                % racuna preko proporcije
                if dy1~=dy2
                    xi = x(logical_index);
                    xi = xi - A(1);
                    dx = A(2)-A(1);
                    b(logical_index) = dy1-(dy1-dy2).*xi./dx;
                    middle = (A(1,4)+A(1,2))/2;
                    % "y" koordinata lijeve ivice presjeka na odstojanju "x"
                    y(logical_index) = middle - b(logical_index)/2; 
                else % u suprotnom, sirina je konstantna na tom segmentu
                    b(logical_index) = dy1;
                    % "y" koordinata lijeve ivice presjeka na odstojanju "x"
                    y(logical_index) = A(1,2);
                end
            end
        end %%
        function strain = strain(this,x)
            %%% STRAIN iz proporcije racuna dilataciju na odstojanju x
            %%% [mm]
            % x - zadana koordinata u mm
            % this.x - x koordinata neutralne ose u mm
            ecu2=this.ecu2;
            strain = ecu2-ecu2*x/this.x;
        end         
        %%
        function [sigma_c] = sigmac(this,x)
            %%% returns stress value in MPa for x in [mm]
            % ec se mnozi sa -1 jer je ec negativno za pritisak, 
            % a trebaju nam pozitivni brojevi za formulu za napon
            % sigmac(epsilon)
            ec = -this.strain(x); 
            sigma_c = this.sigmac_of_eps(ec); % [N/mm^2]
            %fc = sigma_c.*(this.b(x)); % [N]
        end
        function fc = sigmacb(this, x)
            %%% SIGMACB racuna segment sile Fc u liniji na koord. x,
            %%% koristi se u get.Fc i sigmac_moment funkcijama za proracun
            %%% sile Fc i xFc
            sigmac = this.sigmac(x);
            fc = sigmac.*(this.b(x));
        end
        %%
        function sigma_c = sigmac_of_eps(this,ec)
        %%% Returns sigma_c stress [MPa] in concrete for given dillatation ec
            % prelocate memory for sigma_c
            sigma_c = zeros(1,length(ec));
            % concrete classes [12/15 16/20 20/25 25/30 30/37 35/45 40/50 45/55
            % 50/60] N/mm^2 (MPa)
            fck = this.fck; % C30/37
            if fck < 50
                ec2 = 2/1000;
                n = 2;
            else
                ec2 = (2+0.085*(fck-50)^0.53)/1000;
                n = 1.4+23.4*((90-fck)/100)^4;
            end
            fcd = this.fcd; % design concrete strength
            
            % if ec>=2 promils, set sigma_c equal to alpha*fcd constant
            sigma_c(ec>=ec2) = fcd; 
            % calculate sigma_c for the remaining elements based on the parabolic part of
            % the concrete stress-strain graph
            logical_index = (sigma_c == 0);
            sigma_c(logical_index) = fcd*(1-(1-ec(logical_index)/ec2).^n); %[N/mm^2]
        end
        
        function out = sigmac_moment(this,x)
        %%% racuna moment sile dijagrama napona oko gornje ivice grede
        %%% (x=0), na osnovu cega se odredjuje polozaj sile Fc (xFc)
            fc = this.sigmacb(x); % [N]
            out = fc.*x;
        end
    end
    
end

