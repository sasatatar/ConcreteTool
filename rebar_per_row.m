function out = rebar_per_row(ds) % ds precnik sipki armature u [mm]
    %%% Odredjuje maksimalan broj sipki u jednom redu
    global bw; % sirina rebra: bw = 0.4m
    bw_mm = bw * 10; % sirina rebra u mm
    dg = 32; % najvece zrno agregata
    
    global c_nom; % ukupan zastitni sloj betona
    global stirrup; % uzengije
    
    spacing = zeros(3, length(ds)); 
    % uslovi za min. razmak izmedju sipki (prema EC2)
    spacing(1,:) = ds;   % precnik armature
    spacing(2,:) = dg+5; % precnik najveceg zrna agregata + 5mm
    spacing(3,:) = 20;   % 20mm
    global rebar_spacing;
    rebar_spacing = max(spacing); % minimalan razmak izmedju sipki
    
    % maksimalan broj sipki u jednom redu
    out = floor((bw_mm - 2*(c_nom + stirrup) - ds)./(rebar_spacing+ds))+1;
    
end