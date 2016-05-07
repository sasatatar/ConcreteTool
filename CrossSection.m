classdef CrossSection < handle
    %CROSSSECTION klasa predstavlja karakteristike poprecnog presjeka
    %nosaca
    %   Sve duzine su u [mm]
    properties 
        dims = struct('bf', 0, 'hf', 0,... % dimenzije nosaca (struct) 
            'hv', 0, 'bw', 0, 'h', 0);                 
        %   .bf;           % sirina flanse
        %   .hf;           % visina flanse
        %   .hv;           % visina vute (prelaz imzedju flanse i rebra)
        %   .bw;            % sirina rebra
        %   .h;             % ukupna visina nosaca
        fck = 30;           % karakt. cvrtstoca betona [MPa] (C30/37)
        fctm = 2.9;         % cvrstoca betona na zatezanje (potrebno kod proracuna As1,min)
        alpha_cc = 0.85;       % koef. koji uzima u obzir dugorocne negativne 
                            % faktore na cvrstocu betona
        gammac = 1.5;       % koeficijent sigurnosti za beton
        delta = 1;          % faktor redistribucije momenatas
        x;                  % polozaj neutralne ose od vrha presjeka [mm]
        
        % podaci vezano za armaturu
        c_nom = 35;             % zastitni sloj betona [mm]
        Rebars = Rebar.empty;   % matrica sipki armature
        dg = 32;                % najvece zrno agregata  
        ds_max = [19 19];       % maksimalan precnik armature [mm] posebno za zategnutu (1) 
                                % posebno za pritisnutu zonu (2)
        As1_req = 0;            % potrebna povrsina armature [mm^2]
        As2_req = 0;
        fyk = 500;              % karakteristicna cvsrtoca celika [MPa]
        Es = 200000;            % modul elasticnosti celika [MPa]
        gamma_s = 1.15;         % koef. sigurnosti za fyd = fyk / gamma_s
        fyd;                    % racunski dopusteni napon u poduznoj armaturi 
        minRebarSpacing;        % minimalni svijetli razmak izmedju sipki
        xdRatio = 0.45;                % odnos x/d = 0.45 za C<=50/60, 0.35 za C>=55/67
        % uticaji
        Nsd = 0;                % normalna sila koja djeluje u presjeku (pritisak je pozitivan)
        Msd = 0;                    % moment savijanja koji djeluje na presjek [Nmm]
        Vsd = 0;                % transverzalna sila
        Ted = 0;                % moment torzije [Nmm]
        
        % podaci vezano za uzengije - stirrup
        stirrup;           % uzengije [mm]
        fywk = 500;             % karakteristicna cvrstoca za uzengije
        fywd;                   % racunska granica tecenja u uzengijama fywd = fywk / gamma_s
        sl = 0;                      % razmak izmedju uzengija [mm]
        m = 2;                      % sjecnost        
        alpha = 90;                  % ugao uzengija sa horizontalom [deg]
        Asw_req = 0;            % ukupna potrebna kol. poprecne armature [mm^2/mm]
        
        ecu2 = -3.5/1000;       % strain in concrete / dilatacija u betonu
        strainHardening = 0;    % set to 1 to use strain hardening
    end
    %% Dependent
    properties (Dependent)
        % tacke koje definisu poprecni presjek u sledecem formatu:
        % [x11 y11 x12 y12 % x11,y11 koordinate prve konturne tacke
        %  x21 y21 x22 y22...]
        Points;
        fcd; % design concrete strength fcd = alfacc * fck / gammac
        
        centroid;         % center of gravity / teziste presjeka
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
        z;                  % krak unutrasnjih sila [mm]
    end%%
    
    %% Metode za armaturu
    methods       
        %% funkcija za plotanje dijagrama M-fi (moment-zakrivljenost)
        function plotMfi(this)
            ecu2 = this.ecu2;
            ec = [0:-0.25:-3.5];
            ec_yield_zone = [];
            
            % utvrdjivanje dilatacija u celiku pri stanju tecenja
            r = Rebar(this);
            eyd = r.fyd/this.Es;
            
            % formiranje praznih vektora za pohranjivanje vrijednosti M, fi
            % i es
            fi = zeros(1,length(ec));
            Mrd = zeros(1,length(ec));
            es = Mrd;
            % yield flag
            yielded = false; 
            
            wbar = waitbar(0, 'Generisanje dijagrama...');
            % za svaku predefinisanu vrijednost ec, izracunaj moment
            % savijanja i zakrivljenost presjeka
            es1 = 0;
            total = length(ec);
            for i = 1:length(ec)
                this.ecu2 = ec(i)/1000;
                Mrd(i) = this.Mrd/10^6; %[kNm]

                es1 = this.strain(this.xFs1);
                es(i)=es1;
                % predjena je granica tecenja
                if es1 > eyd && yielded == false
                    yielded = true;
                    % za zonu tecenja progusti mrezu tacaka na svakih
                    % 0.0625
                    ec_yield_zone = [ec(i-1):-0.0625:ec(i)];
                end
                fi(i) = -this.ecu2/this.x;
                % update waitbar
                waitbar(i/total, wbar);
            end

            % alociraj prazne vektore za pohranjivanje dodatnih vrijednosti
            % M i fi za guscu mrezu tacaka na potezu gdje se desilo tecenje
            % armature
            Mrd_yield_zone = zeros(1, length(ec_yield_zone));
            fi_yield_zone = Mrd_yield_zone;
            % proracunaj M i fi za dodatnu mrezu tacaka
            waitbar(0, wbar, 'Generisanje dodatnih tacaka u zoni plastifikacije...');
            total = length(ec_yield_zone);
            for i = 1:length(ec_yield_zone)
                this.ecu2 = ec_yield_zone(i)/1000;
                Mrd_yield_zone(i) = this.Mrd/10^6;
                fi_yield_zone(i) = -this.ecu2/this.x;
                waitbar(i/total, wbar);
            end
            close(wbar);
            % reset ecu2
            this.ecu2 = ecu2;
            % dilatacije koje treba oznaciti
            % pronalazi koordinate prije nego sto se niz dopuni
            points = -[3.5]; %2 3 3.5
            x = [];
            y = [];
            for i = 1:length(points)
                point = points(i);
                x = [x fi(ec==point)];
                y = [y Mrd(ec==point)];
            end
            % spoji grube i detaljne vektore vrijednosti i poredaj po
            % velicini
            Mrd_refined = unique([Mrd Mrd_yield_zone]);
            fi_refined = unique([fi fi_yield_zone]);
            % get figure handle if it allready exists
            fig = findobj(0, 'Tag', 'mifi');
            
            % if not, create a new one
            if isempty(fig)
                fig = figure;
                set(fig, 'Tag', 'mifi');
                hold on;
            end
            % set figure background color
            fig.Color = [1 1 1];
            % set that figure as active
            figure(fig);
            
%             subplot(2,1,1);
            % plota M-Fi dijagram
            if this.strainHardening == 0
                line_style = '-';
            else
                line_style = '--';
            end
            line1 = plot(fi_refined,Mrd_refined);
            ax1 = gca;
            set(ax1, 'XGrid', 'on', 'YGrid', 'on');
            % podesavanje oznaka osa
            ax1.XLabel.String = 'Zakrivljenost, \Phi (1/mm)';
            ax1.YLabel.String = 'Moment savijanja, M (kNm)';
            set(line1, 'LineWidth', 1.5, 'LineStyle', line_style, 'Color', 'k'); %, 'Color', 'k'
            
            %return; % izkomentarisati ovu liniju da bi se na dijagramu oznacile tacke za dilatacije u betonu 
            % oznaci pojedine tacke
            markers = plot(x,y);
            set(markers, 'LineStyle', 'none',...
                'MarkerEdgeColor', 'none',...
                'Marker', 'o',...
                'MarkerSize', 6,...
                'MarkerFaceColor', line1.Color);
%             
%             subplot(2,1,2);
%             p2 = plot(es,Mrd, '-sr');
%             p2.MarkerFaceColor = 'red';
%             this.fck = this.fck;
        end
        
        
        %% x koordinata tezista poprecnog presjeka
        function x = get.centroid(this)
            x = integral(@this.Sy, 0, this.dims.h)/this.Ac;
        end
        
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
        
        function fywd = get.fywd(this)
            fywd = this.fywk / this.gamma_s; % [MPa]
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
                rebar = Rebar(this, ds, x, y, row, column, zone);
                this.Rebars(end+1) = rebar;
                % pozvati funkciju za crtanje sipke - dok sredim
                rect = rebar.draw(ax);
                % kontrola kolicine armature
                As_max = 0.04*this.Ac;
                As_tot = this.As1 + this.As2 + rebar.Area;
                if As_tot > As_max
                    msgbox('Kolicina armature je veca od maksimalne dopustene prema EC2.','Presjek je prearmiran', 'help');
                end
            else
                disp('mjesto je zauzeto');
                rect = [];
            end
        end
        
        function rearrangeRebars(this, dh, dy)
            %%% REARRANGEREBARS pravi novi raspored postojece armature u
            %%% zategnutoj zoni usljed promjene visine poprecnog presjeka,
            %%% precnika uzengija, ili zastitnog sloja c_nom
            
            if isempty(this.Rebars)
                return;
            end         
            
            if nargin == 3
                dh = dy;
            else
                dy = 0;
            end
            
            for zone = 1:(nargin-1)
                rebars = findobj(this.Rebars, 'zone', zone);
                % osa simetrije
                ys = this.dims.bf/2;
                % koordinata prve kolone
                y0 = (this.dims.bf-this.dims.bw)/2+this.c_nom+this.stirrup+this.ds_max(zone)/2;
                % provjerava se sipka po sipka
                for i = 1:numel(rebars)
                    r = rebars(i);
                    r.x = r.x - (-1)^zone*dh; % zona 1: +dh; zona 2: -dh
                    if dy == 0 
                        continue;
                    end
                    % usvojeno je da je dy > 0 ako se sirina reda uvecala
                    % ako se sipka nalazi s lijeve strane ose simetrije oduzima se dy
                    if r.y ~= ys
                        r.y = r.y - dy*(ys-r.y)/(ys-y0);
                    end
                end
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
        
        function [As1_pot, As2_pot] = calculateAs(this, d)
            %%% CALCULATEAS racuna potrebnu kolicinu armature 
            %%% Msd - moment savijanja u presjeku u Nmm
            % ako d nije definisano, usvaja se 0.9h
            dims = this.dims;
            if nargin == 1
                d = 0.9*dims.h;
            end
            % povecanje momenta usljed djelovanja normalne sile
            Nsd = this.Nsd;
            Msd = this.Msd + Nsd*(d-this.centroid);
            
            % pretpostavlja se jednostruko armiranje
            As1_pot = 0;
            As2_pot = 0;
            delta = 0.01; % zeljena preciznost (dopustena razlika izmedju Msd i Mrd)
            xmin = 0;
            % gornja granica za x je staticka visina presjeka (od
            % pritisnute ivice do najudaljenijeg reda armature
            xmax = this.xdRatio*(dims.h-this.c_nom-this.stirrup-this.ds_max(1)/2); 
%             esy = this.fyd/this.Es;
%             xmax = abs(d*this.ecu2/(this.ecu2-esy));
            
            % maksimalni moment savijanja sa jednostrukim armiranjem (Mrd
            % za x = xmax)
            this.x = xmax; % maksimalna vrijednost x dozvoljena prema EC2
            %Fc = this.Fc;
            Mrd = getMrd();
            
            % Minimalna i maksimalna kolicina armature:
            % minimalna kolicina armature prema EC2:
            As1_min = max([0.26*this.fctm/this.fyk*dims.bw*d...
                0.0013*dims.bw*d]) % ili d umjesto dims.h?
            % ukupna maskimalna kolicina armature (As1+As2)
            As_max = 0.04*this.Ac;

            % provjera da li je potrebno dvostruko armiranje  
            % pomocni Rebar objekat za ocitavanje napona po visini
            r = Rebar(this);
            if Mrd < Msd
                % krak sila
                z = d-this.xFc;
                Fs1 = Mrd / z - Nsd; % [N]
                dM = Msd - Mrd;
                
                % pretpostavka o polozaju tezista pritisnute armature
                % manja vrijednost od 0.1h i polozaja tezista jednog reda
                % armature precnika 36 mm (najveca moguca armatura)
                d2 = min([0.1*dims.h this.c_nom+this.stirrup+this.ds_max(2)/2]);
                zd = d - d2;
                Fs2 = dM / zd; % [N]
                As1_pot = (Fs1+Fs2) / r.sigmas(d); % mm2
                As2_pot = Fs2 / abs(r.sigmas(d2)); % mm2
                if As2_pot > As1_pot/2
                    message = ['As2 > As1/2, potrebna je veca visina presjeka.'];
                    msgbox(message, 'Upozorenje', 'help');
                    return;
                end
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
            As1_pot = (this.Fc-Nsd)/r.sigmas(d); % [mm2]
            if As1_pot > As_max || As2_pot > As_max
                disp('As > As_max');
                msgbox('As,uk > As,max. Potrebno je usvojiti veci presjek.',...
                    'Prearmiran presjek', 'help');
            elseif As1_pot < As1_min
                As1_pot = As1_min;
            end
            
            % pomocna funkcija koja racuna pribliznu vrijednost reaktivnog
            % momenta Mrd na osnovu poznate sile u betonu (Fc) i kraka
            % unutrasnjih sila (z)
            function Mrd = getMrd()
                xFc = this.xFc;
                z = d-xFc;
                Mrd = this.Fc*z; % [Nmm]
            end
        end
        
        function [out1 out2] = calculateAsw(this, ax)
            %%% CALCULATEASW racuna potreban razmak izmedju uzengija
            % razmak izmedju uzengija [mm]
            sl = this.sl;
            out1 = zeros(1,5);
            out2 = out1;
            Vsd = this.Vsd;
            fcd = this.fcd; % racunska cvrstoca betona
            bw = this.dims.bw; % sirina rebra
            d = this.xFs1; % staticka visina           
            z = 0.9*d; %this.z;
            if z <= 0
                msgbox('Nije definisan unutrasnji krak sila (z), presjek ne sadrzi poduznu armaturu u zoni zatezanja.',...
                    'Presjek nije armiran', 'help');
                return;
            end
            dsw = this.stirrup; % usvojeni precnik uzengija
            m = this.m; % zadana sjecnost
            Asw = m*dsw^2*pi/4; % efektivna povrsina jedne uzengije [mm2]
            alpha = this.alpha; % ugao izmedju uzengija i horizontale [deg]
            fywk = this.fywk;   % karakteristicna cvrstoca poprecne armature
            fywd = fywk/this.gamma_s;   % usvojena proracunska cvrstoca pop. armature
            
            
            % koef. redukcije cvrstoce na pritisak u betonu
            v1 = 0.6*(1-this.fck/250);
            Vrd_max = bw*z*v1*fcd*(1+cotd(alpha))/2; % za theta = 45 stepeni

            % broj tacaka koje ce se generisati
            % sto je veci, veca je preciznost
            num = 500;
            if Vsd <= Vrd_max
                Vrd = unique([linspace(0, Vrd_max, num) Vsd]);
            else
                Vrd = unique(linspace(0, Vrd_max, num));
            end
            Vrd = Vrd(Vrd>0);
            % proracun ugla teta
            b = bw*z*v1*fcd./Vrd;
            cotTheta = (b+sqrt(b.^2-4*(1-b*cotd(alpha))))/2;
            
            % korekcija vrijednosti prema EC2
            cotTheta(cotTheta > 2.5) = 2.5;
            cotTheta(cotTheta < 1) = 1;
            s = Asw*z*fywd*(cotTheta+cotd(alpha))*sind(alpha)./Vrd;
            
            % minimalni i maksimalni dopusteni razmak uzengija
            Rho_w_min = 0.08*sqrt(this.fck)/this.fywk;
            s_max = min(0.75*d*(1+cotd(alpha)), Asw/(Rho_w_min*bw*sind(alpha)));
            s_min = min(s);
            % korekcija vrijednosti s prema EC2
            s(s>s_max) = s_max;
            % al - translacija dijagrama, tj. produzenje zat. arm.
            al = z*(cotTheta-cotd(alpha))./2; % [mm]
            % dodatna sila zatezanja
            dFtd = Vsd.*al/z;
            
            % ukoliko je zadano s, odredjuje se Vrd interpolacijom
            % za interpolaciju je potreban monotono rastuci niz x
            % koordinata, tj. vektor s, tako da za proracun thete
            % koristimo samo taj dio vektora, s_red, cotTheta_red, Vrd_red
            
            % broj jedinstenih elemenata u vektoru s
            elements = numel(unique(s));
            % redukovani vektori sa jedinstvenim vrijednostima
            al_red = al(end-elements+1:end);
            dFtd_red = dFtd(end-elements+1:end);
            s_red =  s(end-elements+1:end);
            cotTheta_red = cotTheta(end-elements+1:end);
            Vrd_red = Vrd(end-elements+1:end);
            Vrd_min = Vrd_red(1);
            
            % output 1 [s cotTheta Vrd dFtd al]
            if Vsd <= 0 || isnan(Vsd)
                out1 = zeros(1, 5);
            elseif Vsd <= Vrd_max
                % ako je zadana sila manja od minimalne, obaviti proracun
                % za Vrd_min
                logical_index = (Vrd==max(Vsd, Vrd_min));
                cotTheta_d = cotTheta(logical_index);
                al_out = z*(cotTheta_d-cotd(alpha))/2;
                dFtd_out = al_out/z*Vsd/1000;
                out1 = [s(logical_index) cotTheta_d max(Vrd_min,Vsd)/1000 dFtd_out al_out];
            else
                out1 = zeros(1,5);
                message = ['Presjek nije u stanju da prenese zadanu transverzalnu silu. '...
                    'Potrebno je usvojiti manji ugao alfa, veci presjek ili visu klasu betona.'];
                msgbox(message, 'Presjek preslab', 'help');
            end
            % output 2 - odredjivanje Vrd i cotTheta za zadanu vrijednost s interpolacijom
            % out2: [s cotTheta Vrd dFtd al]
            if ~isnan(sl)
                if sl == 0
                    out2 = zeros(1,5);
                elseif sl >= s_min && sl <= s_max
                    cotTheta_d = interp1(s_red, cotTheta_red, sl);
                    al_out = z*(cotTheta_d-cotd(alpha))/2;
                    dFtd_out = al_out/z*Vsd/1000;
                    out2 = [sl  cotTheta_d interp1(s_red, Vrd_red, sl)/1000 dFtd_out al_out]; %[mm - kN mm]
                else
                    out2 = zeros(1,5);
                    if sl < s_min
                        message = ['Usvojeni razmak uzengija (s) je manji od minimalno dopustenog pri cemu dolazi do'...
                            ' nenajavljenog loma pritisnute dijagonale u betonu prema EC2 (izraz 6.15) (s(min) = ' sprintf('%.2f',s_min) ' mm). '...
                            'Unesite vrijednost u odgovarajucem opsegu.'];
                        title = 'Prearmiran presjek';
                    end
                    if sl > s_max
                        message = ['Usvojeni razmak uzengija (s) je veci od maksimalno dopustenog prema EC2 (izrazi 9.5N i 9.6N) standardu'...
                            ' (s(max) = ' sprintf('%.2f',s_max) ' mm). '...
                            'Unesite vrijednost u odgovarajucem opsegu.'];
                        title = 'Nedovoljno armiran presjek';
                    end
                    msgbox(message, title, 'help');
                end
            end
            % plotanje Vrd grafa
            % ako nije proslijedjen dijagram, prekini izvrsenje
            if nargin < 2
                return;
            end
            lineColor = 1*[0 0.4471 0.7412];
            cla(ax(1));
            Vrd_line = line('Color', lineColor,... %0.4*[1 1 1] [0.1529 0.2275 0.3725]
                'Tag', 'Vrd_line',...
                'LineWidth', 2,...
                'Parent',ax(1));
            ax(1).XGrid = 'on';
            ax(1).YGrid = 'on';
            ax(1).XLabel.String = 'A_{sw}/s [mm]';
            ax(1).YLabel.String = 'V_{Rd} [kN]';
            XData = Asw./s;
            Vrd_line.XData = XData;
            Vrd_line.YData = Vrd/1000; % [kN]
            space = 0.1*(XData(end)-XData(1));
            ax(1).XLim = [XData(1)-space XData(end)+space];
            % ispisivanje s na x osi umjesto Asw/s
            ax(1).XTickLabel = sprintf('%.1f\n', ax(1).XTick); %1./ax(1).XTick.*Asw
            % plotanje horizontalne crvene linije koja oznacava Vrd_max
            line('XData', ax(1).XLim, 'YData', Vrd_max*[1 1]./1000,...
                'Color', 'red', 'Parent', ax(1));
            
            % plotanje al dijagrama
            cla(ax(2));
            al_line = line('Color', lineColor,...
                'Tag', 'al_line',...
                'LineWidth', 2,...
                'Parent',ax(2));
            ax(2).XGrid = 'on';
            ax(2).YGrid = 'on';
            ax(2).XLabel.String = 'A_{sw}/s [mm]';
            ax(2).YLabel.String = '\DeltaF_{td} [kN]';
            al_line.YData = dFtd_red/1000;
            XData = Asw./s_red;
            al_line.XData = XData;
            space = 0.1*(XData(end)-XData(1));
            ax(2).XLim = [XData(1)-space XData(end)+space];
            % ispisivanje s na x osi umjesto Asw/s
            ax(2).XTickLabel = sprintf('%.1f\n', ax(2).XTick); %1./ax(2).XTick.*Asw
            
            % oznacavanje tacke na dijagramu koja odgovara Vsd (crvena
            % tacka)
            if Vsd ~= 0 
                % Vrd dijagram
                line(Vrd_line.XData(Vrd==min(Vrd_max, Vsd))*[1 1], Vsd/1000*[1 1], 'Parent', ax(1),...
                    'Marker', 'o', 'MarkerFaceColor', 'red',...
                    'MarkerEdgeColor', 'none');
                % dFtd dijagram
                line(Asw/out1(1)*[1 1], out1(4)*[1 1], 'Parent', ax(2),...
                    'Marker', 'o', 'MarkerFaceColor', 'red',...
                    'MarkerEdgeColor', 'none');
            end
            % oznacavanje tacke na dijagramu koja odgovara usvojenom
            % rastojanju sl (zelena tacka)
            if ~isnan(sl) && sl ~= 0
                % Vrd dijagram
                line(Asw/out2(1)*[1 1], out2(3)*[1 1], 'Parent', ax(1),...
                    'Marker', 'o', 'MarkerFaceColor', [0 0.8 0],...
                    'MarkerEdgeColor', 'none');
                % dFtd dijagram
                line(Asw/out2(1)*[1 1], out2(4)*[1 1], 'Parent', ax(2),...
                    'Marker', 'o', 'MarkerFaceColor', [0 0.8 0],...
                    'MarkerEdgeColor', 'none');
            end
        end
        
        % proracun torzione armature
        function [out, s_min] = calculateTrd(this, ax)
            out = zeros(1, 5);
            % racunski moment torzije
            Ted = this.Ted;
            % racunska transverzalna sila
            Ved = this.Vsd;
            dsw = this.stirrup;
            % cvrstoca poprecne armature na zatezanje
            fywd = this.fywd;
            % cvrstoca poduzne armature na zatezanje
            fyd = this.fyd;
            fcd = this.fcd;
            d = this.xFs1; % staticka visina           
            z = 0.9*d; %this.z;

            % povrsina jednog kraka uzengija
            Asw = dsw^2*pi/4; % evektivna povrsina jedne uzengije [mm2]
            % dimenzije rebra
            bw = this.dims.bw;
            h = this.dims.h;
            % povrsina rebra:
            A = bw*h;
            % obim rebra
            u = 2*bw + 2*h;
            % efektivna debljina zida min 2*c_nom+dsw, 
            tef = max([A/u, 2*this.c_nom+dsw]);
            % povrsina Ak
            Ak = (bw-tef)*(h-tef);
            % obim povrsine Ak
            uk = 2*(bw-tef+h-tef);
            % tok smicanja
            q = Ted/2/Ak;
            % generise 100 vrijednosti cotTheta izmedju 1 i 2.5
            cotThetaArray = linspace(2.5, 1, 100);
            % koef. redukcije cvrstoce na pritisak u betonu
            v1 = 0.6*(1-this.fck/250);
            Vrd_max_array = bw*z*v1*fcd*1./(cotThetaArray + 1./cotThetaArray);
            Trd_max_array = v1*fcd*2*Ak*tef./(cotThetaArray + 1./cotThetaArray);
            VTcombined = Ved./Vrd_max_array + Ted./Trd_max_array;
            % plot function
            % clear axes
            cla(ax);
            line(cotThetaArray, VTcombined, 'Parent', ax);
            
            % interpolira vrijednost cotTheta za koju je zbir interakcije 1
            cotTheta = interp1(VTcombined, cotThetaArray, 1);
            
            % provjeriti da li je postojeci presjek dovoljan
            if VTcombined(end) > 1
                message = ['Potreban je jaci presjek. Ved/Vrd + Ted/Trd > 1'];
                msgbox(message, 'Upozorenje', 'help');
                return;
            end
            
            % u slucaju da je presjek prejak pa se ni jednu vrijednost ne
            % prekoraci vrijednost 1 interakcione funkcije, treba usvojiti
            % cotTheta = 2.5 (minimalna vrijednost ugla)
            if isnan(cotTheta)
                cotTheta = 2.5;
            end
            % PRORACUN POPRECNE ARMATURE
            dsw = this.stirrup; % usvojeni precnik uzengija
            m = this.m; % zadana sjecnost poprecne armature za Ved
            Asw = dsw^2*pi/4; %  povrsina jedne uzengije [mm2]
            
            function [Asw_Ved Asw_Ted] = calcAsw(cotTheta)
                % Asw_sum = Asw_Ved + 2Asw_Ted
                % potrebna kolicina poprecne armature za Ved po metru duznom
                % grede
                Asw_Ved = Ved/(fywd*z*cotTheta);
                % potrebna kolicina poprecne armature za Ted po metru duznom
                % grede 
                Asw_Ted = Ted/(fywd*2*Ak*cotTheta);
            end
            [Asw_Ved Asw_Ted] = calcAsw(cotTheta); 

            % Sabiranje poprecne armature
            % potrebna kolicina konturnih uzengija
            % Asw_Ted se mnozi sa 2 jer je formula izvedena za jednu stranicu, 
            % tj. jedan krak uzengija
            Asw_out = Asw_Ved/m*2 + 2*Asw_Ted;
            % potrebna kolicina unutrasnjih uzengija
            Asw_in = Asw_Ved/m*(m-2);
            % ukupno
            Asw_sum = Asw_in + Asw_out;
            this.Asw_req = Asw_sum;
            
            % potreban razmak izmedju konturnih uzengija
            s_out = 2*Asw/Asw_out;
            % potreban razmak izmedju unutrasnjih uzengija
            if Asw_in > 0
                s_in = (m-2)*Asw/Asw_in;
            else
                s_in = 0;
            end
            
            % minimalni i maksimalni dopusteni razmak uzengija
            % alpha = 90, theta = 45 -> cotTheta = 1
            [maxAsw_Ved maxAsw_Ted] = calcAsw(1);
            maxAsw_out = maxAsw_Ved/m*2 + 2*maxAsw_Ted;
            Rho_w_min = 0.08*sqrt(this.fck)/this.fywk;
            s_max = min([0.75*d m*Asw/(Rho_w_min*bw)]); % u/8 bw
            s_min = ceil(2*Asw/maxAsw_out/5)*5;
            
            if s_out > s_max
                s_out = s_max;
            end
            
            out = [s_out s_in Asw_out Asw_in Asw_sum];
            
            % dodatna poduzna armatura za torziju
            %Asl = q*uk*cotTheta/fywd;
            % potrebna poprecna armatura
        end
        
        function out = recalculateAsw(this, s)
            out = zeros(1,5);
            Asw_sum = this.Asw_req;
            m = this.m;
            Asw = this.stirrup^2*pi/4;
            Asw_out = 2*Asw/s;
            n = 0;
            if Asw_out >= Asw_sum
                Asw_in = 0;
                nmax = 0;
            else 
                Asw_in = Asw_sum - Asw_out;
                % broj razmaka s izmedju dopunskih unutrasnjih uzengija
                n = ((m-2)*Asw)/(Asw_in*s);
                nmax = floor(n);
                % korigovano Asw_in
                Asw_in = (m-2)*Asw/(nmax*s);
            end
            
            out = [s nmax*s Asw_out Asw_in (Asw_out+Asw_in)];
        end
        
    end 
    
    
    %% SET i GET metode
    methods   
        %% set Msd
        function set.Msd(this, Msd)
            if ~isnan(Msd) && Msd > 0
                this.Msd = Msd;
            else
                this.Msd = 0;
            end
        end
        
        %% set Ted
        function set.Ted(this, Ted)
            if ~isnan(Ted) && Ted > 0
                this.Ted = Ted;
            else
                this.Ted = 0;
            end
        end
        
        function set.fck(this, fck)
            this.fck = fck;
            % maksimalna dopustena dilatacija u betonu
            % usvaja se da je negativna
            if this.fck <= 50
                this.ecu2 = -3.5/1000;
            else
                % formula iz EC2
                this.ecu2 = -(2.6+35*((90-this.fck)/100)^4)/1000;
            end
        end
        
        function set.stirrup(this, stirrup)
            if ~isempty(this.stirrup)
                ds = this.stirrup - stirrup;
                this.rearrangeRebars(0, ds); 
            end
            this.stirrup = stirrup;
        end
        
        %% c_nom setter function
        function set.c_nom(this, c_nom)
            if ~isempty(this.c_nom)
                dc = this.c_nom - c_nom;
                this.rearrangeRebars(0, dc);
            end
            this.c_nom = c_nom;
        end
        
        function set.dims(this, dims)
            fields = {'bf','hf','hv','bw','h'};
            changes = zeros(1,5);
            for i = 1:numel(fields)
                field = fields{i};
                changes(i) = eval(sprintf(['dims.' field '- this.dims.' field]));
            end
            this.dims = dims;
            % ukoliko je promijenjena samo visina presjeka, sipke se
            % rearanziraju, inace se brisu
            dh = changes(5); % promjena visine presjeka
            if dh ~= 0
                this.rearrangeRebars(dh);
            elseif sum(changes) ~= 0
                this.Rebars = Rebar.empty();
            end
        end
        
        
        
        function set.Nsd(this, Nsd)
            if isnan(Nsd)
                this.Nsd = 0;
            else
                this.Nsd = Nsd;
            end
        end 
        
        function set.xdRatio(this, xdRatio)
            %%% odnos x/d = 0.45 za C<=50/60, 0.35 za C>=55/67
            this.xdRatio = xdRatio;
            if this.fck <= 50 && xdRatio > 0.45
                msgbox('Odnos x / d za ovu klasu betona ne bi trebao biti veci od 0.45.','Paznja', 'help');
            elseif this.fck > 50 && xdRatio > 0.35
                msgbox('Odnos x / d za ovu klasu betona ne bi trebao biti veci od 0.35.','Paznja', 'help');
            end
        end
        
        function set.x(this, x)
            this.x = x;
            % FC racuna intenzitet rezultantne sile u betonu Fc [N]
            Fc = abs(integral(@this.sigmacb, 0, x)); % N
            this.Fc = Fc;
            % XFC racuna polozaj tezista sile u betonu Fc racunajuci od
            % vrha presjeka [mm]
            this.xFc = abs(integral(@this.sigmac_moment, 0, x)/Fc);  %[mm]
        end
        
        function x = get.x(this)
            x = this.x;
        end
        
        function p = get.Points(this)
            %%% POINTS generise matricu tacaka koje definisu poprecni
            %%% presjeka nosaca na osnovu zadanih dimenzija (tacke su
            %%% potrebne za graficki prikaz nosaca). 
            [x,y] = this.generatePoints(this.dims);
            p.x = x;
            p.y = y;
        end
        
        % pomocna funkcija koja generise tacke presjeka
        function [x,y] = generatePoints(this,dim)
            % x koordinate
            x = [0 dim.hf dim.hf+dim.hv dim.h];
            x = [x flip(x) x(1)];
            % y koordinate
            y = [0 0 (dim.bf-dim.bw)/2*[1 1]];
            y = [y flip([y(1:2)+dim.bf y(3:4)+dim.bw]) y(1)];
        end
        
        function [x,y] = offsetPolygon(this, dims, c_nom)
            %%% offsetPolygon funkcija racuna tacke presjeka definisanog
            %%% dimenzijama dims transliranim unutra za vrijednost c_nom,
            %%% koristi se kod iscrtavanja uzengija i zastitnog sloja
            %%% betona
            % kopija dimenzija
            newDims = dims;
            % podesavanje novih dimenzija da se dobije offsetovan T presjek
            newDims.h = dims.h - 2*c_nom;
            newDims.bf = dims.bf - 2*c_nom;
            newDims.bw = dims.bw - 2*c_nom;
            % sirina vute
            c = (dims.bf-dims.bw)/2;
            % alpha: ugao izmedju vute i horizontale
            alpha = atan(dims.hv/c);
            % ugao izmedju vute i vertikalne ivice flanse
            beta = pi/2 + alpha;
            % promjena visine flanse hf u zavisnosti od ugla vute
            dhf = c_nom + c_nom/tan(beta/2);
            newDims.hf = dims.hf - dhf;
            % nove tacke presjeka
            [x,y] = this.generatePoints(newDims);
            % translacija tacaka za vrijednost c_nom
            x = x+c_nom;
            y = y+c_nom;
        end
        
        function fcd = get.fcd(this)
            fcd = this.alpha_cc*this.fck / this.gammac;
        end
        
        function fctm = get.fctm(this)
            % cvrstoca na zatezanje
            if this.fck <= 50
                fctm = 0.3*this.fck^(2/3);
            else
                fctm = 2.12*log(1+(this.fck+8)/10);
            end
        end
        
        function Fc = get.Fc(this)
            Fc = this.Fc;
        end
        
        function xFc = get.xFc(this)
            xFc = this.xFc;
        end
        
        function Fs2 = get.Fs2(this)
            %%% FS2 racuna rezultantnu silu u celiku u pritisnutoj zoni
            rebars = findobj(this.Rebars, 'zone', 2);
            Fs2 = sum([rebars.Fs]); 
        end
        
        function xFs2 = get.xFs2(this)
            %%% XFS2 racuna polozaj sile Fs2 na x osi mjereno od
            %%% vrha presjeka, 
            %%% if no compression reinforcement, returns 0
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
            if isnan(xFs1)
                xFs1 = 0;
            end
        end
        
        function z = get.z(this)
            xFs1 = this.xFs1;
            if xFs1 == 0
                z = 0;
                return;
            end
            Fc = this.Fc;
            Fs2 = this.Fs2;
            x = (Fc*this.xFc + Fs2*this.xFs2)/(Fc + Fs2);
            z = xFs1 - x;
        end
        
        % get.Mrd funkcija racuna reaktivni moment savijanja presjeka, na
        % osnovu usvojene armature u [Nmm]
        function Mrd = get.Mrd(this)
            if isempty(findobj(this.Rebars, 'zone', 1))
                Mrd = 0;
                return;
            elseif this.ecu2 == 0
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
            Mrd = this.Fc*(this.xFs1-this.xFc)+this.Fs2*(this.xFs1-this.xFs2)-Nsd*(this.xFs1-this.centroid); % [Nmm]
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
                    this.stirrup = 8;
                end
                if n == 2
                    b = d(1);
                    h = d(2);
                    this.dims.bf = b;
                    this.dims.hf = 0;
                    this.dims.hv = 0;
                    this.dims.bw = b;
                    this.dims.h = h;
                    this.stirrup = 8;
                end
            else % default dimenzije presjeka
                this.dims.bf = 1350;
                this.dims.hf = 100;
                this.dims.hv = 200;
                this.dims.bw = 400;
                this.dims.h = 1200;
                this.stirrup = 10;
            end
            % pretpostavka polozaja neutralne ose
            this.x = this.dims.h*0.9*this.xdRatio;
        end
        
        function plotSection(this,ax)
            %%% PLOTSECTION plota presjek nosaca sa uzengijama
            % uklanjanje postojecih linija sa grafa
            %delete(findobj(ax.Children, 'flat', 'Tag', 'crossSection_line'));
            %delete(findobj(ax.Children, 'flat', 'Tag', 'stirrup_patch'));
            cla(ax);            
            ax.CameraUpVector = [-1 0 0]; % rotira koord. sistem za 90 stepeni u smjeru kazaljke
            ax.DataAspectRatio = [1 1 1]; % podesava odnos dimenzija x,y i z osa da je 1:1:1
            ax.XLim(2) = this.dims.h; % podesava opseg na x osi
            
            section_line = line('Color', [0.4 0.4 0.4],...
                'Tag', 'crossSection_line',...
                'LineWidth', 1.3,...
                'Parent',ax);
            % x i y koordinate tacaka koje se spajaju linijom
            x = this.Points.x;
            y = this.Points.y;
            %ax.YLim = [min(y)-50, max(y)+50];
            section_line.XData = x;
            section_line.YData = y;
            
            % oznake dimenzija presjeka na x osi
            ticksX = unique(x);
            ax.XTick = ticksX(ticksX~=0);
            ax.XTickLabel = sprintf('%.1f\n', ax.XTick);
            
            % oznake dimenzija presjeka na y osi
            ticksY = unique(y);
            ax.YTick = ticksY;
            ax.YTickLabel = sprintf('%.1f\n', ax.YTick);
            
            % uzengije
            c_nom = this.c_nom;
            dims = this.dims;
            stirrup = this.stirrup;
            % koordinate spoljasnje konture uzengija na odstojanju c_nom od
            % ivice presjeka
            [xo, yo] = this.offsetPolygon(dims, c_nom);
            % koordinate unutrasnje konture uzengija na odstojanju
            % c_nom+stirrup od ivice presjeka
            [xi, yi] = this.offsetPolygon(dims, c_nom+stirrup);
            
            stirrup_x = [xo flip(xi)];
            stirrup_y = [yo flip(yi)];
            
            % konturne uzengije
            patch('FaceColor', 0.4*[1 1 1],...
                'Tag', 'stirrup_patch',...
                'Parent', ax,...
                'LineStyle', 'none',...
                'XData', stirrup_x, 'YData', stirrup_y); 
            
            % uzengije u rebru presjeka
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

            patch('FaceColor', 0.4*[1 1 1],...
                'Tag', 'stirrup_patch',...
                'Parent', ax,...
                'LineStyle', 'none',...
                'XData', xs, 'YData', ys);
            
            % horizontalni krak uzengija uz gornju ivicu
%             xs = [0 0];
%             ys = [0 dims.bf];
%             transx = [1 1];
%             transy = [1 -1];
%             xs = xs+c_nom*transx;
%             ys = ys+c_nom*transy;
%             xs = [xs flip(xs)+stirrup*flip(transx)];
%             ys = [ys flip(ys)]; % +stirrup*flip(transy)
%             patch('FaceColor', 0.4*[1 1 1],...
%                 'Tag', 'stirrup_patch',...
%                 'Parent', ax,...
%                 'LineStyle', 'none',...
%                 'XData', xs, 'YData', ys);  
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
                'FaceAlpha', 0.25);
            % x i y koordinate tacaka koje se spajaju linijom
            x = this.Points.x(1:end-1);
            y = this.Points.y(1:end-1);
            % obiljezavanje pritisnute povrsine betona
            [bx,y_left] = this.b(this.x);
            
            % tacke A i B definisu neutralnu liniju 
            % provjeriti kriticni slucaj kada je neutralna usa u zoni donje
            % flanse (za slucaj da je nosac u obliku I presjeka)
            % neutralna linija
            nline_x = this.x*[1 1];
            nline_y = [y_left y_left+bx];
            A = [this.x,y_left];
            B = [this.x,y_left+bx];
            logical_index = x<=this.x;
            xc = x(logical_index);
            yc = y(logical_index);
            middle_index = length(xc)/2;
            xc = [xc(1:middle_index) nline_x xc(middle_index+1:end)];
            yc = [yc(1:middle_index) nline_y yc(middle_index+1:end)];

            compressedZone.XData = xc;
            compressedZone.YData = yc;
        end %%
        
        function plotStrain(this, ax, section_axes)
            %%% PLOTSTRAIN generise grafik dilatacija
            
            % calculate x and y data for patch object
            zone = 1;
            % x koordinata prvog reda sipki zategnute armature (zona 1)
            row1 = this.dims.h - this.c_nom - this.stirrup - this.ds_max(1)/2;
            x = [0 row1 row1 0];
            y = [this.ecu2 this.strain(row1) 0 0];
            
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
            % privremeno: pravim pomocni Rebar objekat da ocitam eud
            r = Rebar(this);
            % kasnije: euk se cuva unutar crossSection objekta, i sipka ga
            % cita odatle (trenutno je pohranjen unutar Rebar klase
            if this.xFs1 ~= 0 
                if this.x/row1 > this.xdRatio || ...
                        (this.strainHardening == true && this.strain(row1) > r.eud)
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
           x = this.dims.h-this.c_nom-this.stirrup-this.ds_max(1)/2;
           text(x, -1/1000,...
               num2str(this.strain(x)*1000, '%0.2f'),...
               'Parent', ax,...
               'Color', [0.3059 0.3961 0.5804]); 
           
           % oznacavanje x koordinata
           ax.XTick = [0 this.x row1];
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
                'XLim', section_axes.XLim, 'YLim', [0 this.fcd]);
            
            % plotanje dijagrama napona u betonu
            numPoints = 100;
            x = linspace(0, this.x, numPoints);
            y = this.sigmac(x);
            % dodavanje koordinata pocetne tacke na kraju niza da se
            % zatvori kontura
            x = [x 0];
            y = [y 0];
            
            patch('FaceColor', 0.4*[1 1 1],...%0.8627*[1 1 1],...
                'Tag', 'stress_patch',...
                'Parent', ax,...
                'FaceAlpha', 0.3,...
                'LineStyle', 'none',...
                'XData', x, 'YData', y);
            
            % oznacavanje koordinata polozaja sila i n. ose na x osi
            ticks = unique([this.xFc this.xFs1 this.xFs2]);
            % ukoliko je normalna sila Nsd razlicita od nule, dodati i njen
            % polozaj
            if this.Nsd ~= 0
                ticks = sort([ticks this.centroid]);
            end
            ax.XTick = ticks(ticks~=0);
            ax.XTickLabel = sprintf('%.1f\n', ax.XTick);
            
            % oznacavanje fcd na y osi
            ax.YTick = [0 this.fcd];
            ax.YTickLabel = sprintf('%.1f\n', ax.YTick);
            
            
            % enum forceKind: 1 - beton, 2 - celik, 3 normalna sila
            % ucrtavanje sile u betonu
            forceKind = 1;
            x = this.xFc;
            drawVector(x, '<', forceKind, this.Fc);
            
            % ucrtavanje sile u zategnutoj armaturi
            if ~isempty(findobj(this.Rebars, 'zone', 1))
                forceKind = 2;
                x = this.xFs1;
                drawVector(x, '>', forceKind, this.Fs1);
            end
            
            % ucrtavanje sile u pritisnutoj armaturi
            if ~isempty(findobj(this.Rebars, 'zone', 2))
                %fcolor = [0 0.4 0.6];
                forceKind = 2;
                x = this.xFs2;
                drawVector(x, '<', forceKind, this.Fs2);
            end
            
            % ucrtavanje normalne sile
            if this.Nsd ~= 0
                forceKind = 3;
                x = this.centroid;
                drawVector(x, '>', forceKind, this.Nsd);
            end
            
            % pomocna funkcija za ucrtavanje vektora sila
            function drawVector(x, triangleStyle, forceKind, force)
                % crta silu
                % definisanje boje
                if forceKind == 2
                    color = [0 0.6 0.6]; % boja za Fs
                else
                    color = 0.2*[1 1 1]; % boja za Fc i N
                end
                fx = x;
                k = 0.6;
                fy = k*ax.YLim(2);
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
                    dy = k*ax.YLim(2);
                end
                p.YData = dy*[1 1]; % + p.MarkerSize/ax.OuterPosition(3)*ax.YLim(2)*[1 1];
                ax.Units = 'pixels';
                
                % vrijednost sile [kN] - tekst
                label = text(x, 0, [num2str(force/1000, '%0.2f') ' kN'],...
                    'Parent', ax,...
                    'Color', color); % numericki prikaz dilatacije u betonu u promilima
                % mijenjam units u pixels da bih mogao definisati apsolutan polozaj taga
                % nezavisno od razmjere grafika
                label.Units = 'pixels';
                % pomjera oznaku udesno za 5 px u odnosu na
                % prvobitni polozaj
                if forceKind == 1 || forceKind == 3
                    dx = -10; % ako je u pitanju Fc ili N, pomjeri oznaku nadole
                else
                    dx = 10; % ako je u pitanju Fs, pomjeri oznaku nagore
                end
                label.Position(1:2) = [label.Position(1)+5 label.Position(2)+dx];
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
        %% pomocna funkcija koja racuna staticki moment povrsine oko y ose
        % koristi se u proracunu x koordinate tezista poprecnog presjeka
        function Sy = Sy(this, x)
            b = this.b(x);
            Sy = b.*x;
        end
        
        %% daje sirinu nosaca b na odstojanju x od vrha u mm
        %% output: sirina b i lijeva y koordinata u odgovarajucem presjeku
        function [b,y] = b(this,x) 
            % x koordinate presjeka
            x1 = this.Points.x;
            % y koordinate presjeka
            y1 = this.Points.y;
            % horizontalna prava na zadanoj visini x
            % generise po 2 x i dvije y koordinate za svaku vrijednost x, koje definisu 
            % horizontalne linije kroz konturu i razdvaja ih s NaN, radi vektorizacije koda
            x2 = reshape([x; x; NaN*ones(1,length(x))], 1, 3*length(x));
            % svakoj vrijednosti x pridruzuje y koordinate krajnje lijeve (0) i
            % krajnje desne (bf) tacke presjeka umanjeno, tj. uvecano za 10
            % radi malih odstupanja u dimenzijama grede usljed gresaka
            % zaokruzivanja
            y2 = repmat([-10 this.dims.bf+10 NaN], 1, length(x));
            
            % presjecne tacke horizontalne prave i konture presjeka u
            % obliku vektor kolona x i y
            [xo,yo] = polyxpoly(x1, y1, x2, y2);
            % trazene sirine presjeka
            b = [yo(2:2:end)-yo(1:2:end)]';
            % lijeve y koordinate
            y = yo(1:2:end);
        end %%
        function strain = strain(this,x)
            %%% STRAIN iz proporcije racuna dilataciju na odstojanju x
            %%% [mm]
            % x - zadana koordinata u mm
            % this.x - x koordinata neutralne ose u mm
            ecu2=this.ecu2;
            strain = ecu2-ecu2*x./this.x;
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
            
            % if ec>=2 promils, set sigma_c equal to fcd constant
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
    
    methods (Access = 'private')
        function [x0 y0] = origin(this, zone)
            %%% ORIGIN(zone) racuna koordinate prvog reda i prve kolone
            %%% sipki armature
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
        end
    end
    
end

