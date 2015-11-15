function [x_Fs As pref] = pack_rebar(As_req)
    %%% PCAK_REBAR v2 chooses the best combination of rebars and returns their center
    %%% of gravity, cross section and best choice (closest to the assumed
    %%% value of d
    global h;      % height of cross section in [mm]
    h_mm = h*10;
    global ds = [12 14 16 19 22 25 28 32 36] % rebar diameters [mm]
    c_nom = 37;  % concrete cover  [mm]
    global stirrup % uzengije     12 [mm]
    global rebar_spacing;
    global d;      % x_Fs assumption
    global fyd;    % dopusteni racunski napon u celiku
   
    
    rebar_csection = ds.^2*pi/4/100; % [cm^2]
    
    % get the maximum number of rebars per row for each rebar diameter
    rebars_per_row = rebar_per_row(ds);
    
    rebar_req = ceil(As_req./rebar_csection); % rebars required
    
    % number of rows required
    rebar_rows_req = ceil(rebar_req./rebars_per_row); 
    
    % number of rebars in top row (as they are mostly left unfilled)
    last_row = rebars_per_row - (rebar_rows_req .* rebars_per_row - rebar_req);
    
    % number of blank rebar spaces in last row
    row_difference = rebars_per_row - last_row;
    
    % calculate the center of gravity of rebars (x_Fs) and their total area
    % As
    x_Fs = zeros(1, length(ds));
    As = rebar_req.*rebar_csection;
    for i=1:length(rebar_rows_req)
        x_row = zeros(1, rebar_rows_req(i)); % x coords for each row
        rebars_in_row = x_row;
        x_row(1) = h_mm-c_nom(i)-stirrup-ds(i)/2;
        if rebar_rows_req(i) == 1
            rebars_in_row(1) = rebar_req(i);
        else
            rebars_in_row(1) = rebars_per_row(i);
        end
        for j=2:rebar_rows_req(i)
            % x coord. of the next row obtained by subtracting rebar
            % spacing from the x coord of the row underneath
            x_row(j) = x_row(j-1) - (rebar_spacing(i) + ds(i));
            if x_row(j) < 0
                x_Fs(i) = -1;
                break;
            end
            % set nuber of bars in last row equal to last_row(i)
            if j == rebar_rows_req(i)
                rebars_in_row(j) = last_row(i);
            else % otherwise, maximal number of rebars per row
                rebars_in_row(j) = rebars_per_row(i);
            end
        end
        x_Fs(i) = x_row * rebars_in_row'/sum(rebars_in_row);
        %es = -ec_of_x(x_row/10);
        %Fs = rebars_in_row.*rebar_csection(i).*es.*
    end
    
    % preffered choice of rebar
    assumption_diff = abs(x_Fs-d_mm);
    [~, pref] = min(assumption_diff);
    while last_row(pref)==0
        assumption_diff(pref) = [];
        [~, pref] = min(assumption_diff);
    end
    if isempty(pref)
        error('Ni jedna sema armiranja ne odgovara, usvojiti armaturu rucno');
    end
    
    x_Fs = x_Fs/10;
    
end