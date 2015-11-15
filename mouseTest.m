function window_motion_test
figure('WindowButtonDownFcn','')
% SortMethod = 'childorder' ubrzava render
ah = axes('SortMethod','childorder', 'ButtonDownFcn', @wbdcb);
ah.DataAspectRatio = [1 1 1];

axis ([1 10 1 10])
title('Click and drag')
   function wbdcb(src,callbackdata)   
       % src - window handle
     % vrsta klika
     f = gcf;
       seltype = f.SelectionType;
       
     % ako je korisnik napravio lijevi klik 
     % SelectionType = normal
     if strcmp(seltype,'normal')
        % strelicu za pointer zamijeni krugom
        %src.Pointer = 'circle';
        % trenutna koordinata misa u axes objektu
        cp = ah.CurrentPoint;
        x = cp(1,1);
        y = cp(1,2);
        w = 0.3;
        h = 0.3;
        rectangle('Position', [x-w/2 y-h/2 w h], 'Curvature', [1 1],...
            'FaceColor', 'red',...
            'ButtonDownFcn', @onClick);
        % cuvanje pocetnog polozaja 
%         xinit = cp(1,1);
%         yinit = cp(1,2);
%         hl = line('XData',xinit,'YData',yinit,...
%         'Marker','p','color','b');
        % funkcija koja se poziva pri svakom pomjeranju misa
%         src.WindowButtonMotionFcn = @wbmcb;
%         % button release funkcija
%         src.WindowButtonUpFcn = @wbucb;
     end    
       function onClick(this, ~)
           seltype = f.SelectionType;
           if strcmp(seltype,'alt')
               delete(this);
               return;
           end
       end
        % funkcija koja se poziva pri svakom pomjeranju misa
        function wbmcb(src,callbackdata)
           cp = ah.CurrentPoint;
           xdat = [xinit,cp(1,1)];
           ydat = [yinit,cp(1,2)];
           hl.XData = xdat;
           hl.YData = ydat;
           % Update figure
           drawnow
        end
        % button release funkcija
        function wbucb(src,callbackdata)
           % cuva koji je taster misa pritisnut 
           % normal - lijevi
           % alt - desni
           last_seltype = src.SelectionType;
            
           % ako je zadnji klik bio desni
           if strcmp(last_seltype,'alt')
              % vrati strelicu za pointer misa
              src.Pointer = 'arrow';
              % obrisi callback funkcije za trenutni prozor
              src.WindowButtonMotionFcn = '';
              src.WindowButtonUpFcn = '';
           else
               % ako je bio lijevi klik, ne radi nista, izadji iz funkcije
              return
           end
        end
  end
end