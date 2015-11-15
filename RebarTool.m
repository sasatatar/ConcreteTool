function varargout = RebarTool(varargin)
% REBARTOOL MATLAB code for RebarTool.fig
%      REBARTOOL, by itself, creates a new REBARTOOL or raises the existing
%      singleton*.
%
%      H = REBARTOOL returns the handle to a new REBARTOOL or the handle to
%      the existing singleton*.
%
%      REBARTOOL('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in REBARTOOL.M with the given input arguments.
%
%      REBARTOOL('Property','Value',...) creates a new REBARTOOL or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before RebarTool_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to RebarTool_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help RebarTool

% Last Modified by GUIDE v2.5 03-Nov-2015 17:03:17

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @RebarTool_OpeningFcn, ...
                   'gui_OutputFcn',  @RebarTool_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before RebarTool is made visible.
function RebarTool_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to RebarTool (see VARARGIN)

% Choose default command line output for RebarTool
handles.output = hObject;
hObject.CloseRequestFcn = @close;
handles.ds = [12 14 16 19 22 25 28 32 36];
% Update handles structure
guidata(hObject, handles);
% provjera da li je proslijedjen CrossSection objekat zajedno sa ostalim 
% handles-ima
if nargin>=4
    % preuzimanje potrebnih referenci za objekte
    ConcreteToolHandles = varargin{1};
    % crossSection handle
    handles.crossSection = ConcreteToolHandles.crossSection;
    % As_uitable handle iz ConcreteTool prozora
    handles.As_uitable = ConcreteToolHandles.As_uitable;
    % section_axes iz ConcreteTool prozora
    handles.section1_axes = ConcreteToolHandles.section_axes;
    % Update handles structure
    guidata(hObject, handles);
    
    % popunjavanje tabele rebar_uitable
    % table handle
    cs = handles.crossSection;
    table = handles.rebar_uitable;
    table.Position = [297 422 383 186];
    % podrzane dimenzije sipki armature
    ds = handles.ds; % [12 14 16 19 22 25 28 32 36]
    % podesavanje da su nazivi po redovima precnici armature
    set(table, 'RowName', ds);
    tableData = zeros(length(ds), 5); 
    % definise N/red i As kolone
    zone = 1; % prvo otvaramo zategnutu zonu
    tableData(:,1) = [cs.RPR(ds, zone)]'; % N/red 
    tableData(:,2) = [ds.^2*pi/4]'; % As
    tableData(:,3) = tableData(:,1).*tableData(:,2); % As/red
    % Potrebno
    As1_req = cs.As1_req; % str2double(handles.As1_edit.String)
    % potreban broj sipki
    tableData(:,4) = ceil(As1_req./tableData(:,2));
    % As_uk
    tableData(:,5) = tableData(:,2).*tableData(:,4);
    table.Data = tableData;
    % formatiranje kolona (numeric za cijele brojeve, bank za 2 decimale)
    set(table, 'ColumnFormat',{'numeric', 'bank', 'bank', 'numeric', 'bank'});
    
    % unos podataka u As1_uitable
    table = handles.As1_uitable;
    tableData = zeros(3, 2);
    tableData(1, :) = [cs.As1_req cs.As2_req];
    tableData(2, :) = [cs.As1 cs.As2];
    tableData(3, :) = tableData(1, :) - tableData(2, :);
    table.Data = tableData;
    table.RowName = {'Potrebno', 'Ugradjeno', 'Nedostaje'};
    table.ColumnName = {'As1 [mm2]', 'As2 [mm2]'};
    table.ColumnFormat = {'bank', 'bank'};
    % fino podesavanje dimenzija tabele za As
    table.Position(3:4) = [253 76];
    %handles.As1_req_edit.String = num2str(cs.As1_req,'%.2f');
    
    % setup dsmax_popup 
    zone = 1;
    index = find(ds==cs.ds_max(zone),1);
    handles.dsmax_popup.String = ds;
    handles.dsmax_popup.Value = index;
    
    % setup ds_popup
    setDsPopup(handles);
    
    % prikaz presjeka
    ax = handles.section_axes;
    cs.plotSection(ax);
    % Ogranicava prikaz na zategnutu zonu presjeka
    ax.XLim = [0.65*cs.dims.h cs.dims.h];
    ax.YLim = [(cs.dims.bf-cs.dims.bw)/2 (cs.dims.bf-cs.dims.bw)/2+cs.dims.bw];
    % prikaz pomocnih linija
    cs.plotGrid(ax, zone);
    
    % prikaz sipki armature
    cs.drawRebar(ax);
    % dodaje right click funkcionalnost nacrtanim objektima, da bi se
    % ranije nacrtane sipke mogle brisati (nakon ponovnog otvaranja
    % RebarTool prozora)
    rebars = findobj(ax, 'Type', 'rectangle');
    for i = 1:numel(rebars)
        set(rebars(i), 'ButtonDownFcn', {@RebarOnClick, handles});
    end
else
    disp('Nije definisan poprecni presjek');  
end

function close(obj, eventdata, handles)
delete(obj);


% podesava popup meni za odabir sipke za ugradnju
function setDsPopup(handles)
% setup ds_popup
cs = handles.crossSection;
zone = handles.zone_popup.Value;
index = find(handles.ds==cs.ds_max(zone), 1);
handles.dsmax_popup.Value = index;
ds = handles.dsmax_popup.String;
handles.ds_popup.String = handles.ds(1:index);
handles.ds_popup.Value = index;


% --- Outputs from this function are returned to the command line.
function varargout = RebarTool_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
keyboard


% --- Executes on selection change in dsmax_popup.
function dsmax_popup_Callback(hObject, eventdata, handles)
% hObject    handle to dsmax_popup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns dsmax_popup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from dsmax_popup
cs = handles.crossSection;
zone = handles.zone_popup.Value;
index = handles.dsmax_popup.Value;
cs.ds_max(zone) = str2double(handles.dsmax_popup.String(index,:));
% podesava popup meni za odabir sipke za ugradnju
setDsPopup(handles);
reset_pushbutton_Callback(0, 0, handles);



% --- Executes during object creation, after setting all properties.
function dsmax_popup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dsmax_popup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function As1_req_edit_Callback(hObject, eventdata, handles)
% hObject    handle to As1_req_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of As1_req_edit as text
%        str2double(get(hObject,'String')) returns contents of As1_req_edit as a double


% --- Executes during object creation, after setting all properties.
function As1_req_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to As1_req_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in ds_popup.
function ds_popup_Callback(hObject, eventdata, handles)
% hObject    handle to ds_popup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns ds_popup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from ds_popup


% --- Executes during object creation, after setting all properties.
function ds_popup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ds_popup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on mouse press over axes background.
function section_axes_ButtonDownFcn(ax, eventdata, handles)
% hObject    handle to section_axes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
fig = handles.RebarTool_figure;

% left click action
if strcmp(fig.SelectionType, 'normal')
    cs = handles.crossSection;
    % ocitava index odabrane vrijednosti ds
    index = handles.ds_popup.Value;
    % na osnovu ocitanog index-a, cita vrijednost i pretvara je u broj
    ds = str2double(handles.ds_popup.String(index,:));
    % ocitava koordinate misa u ax koord. sistemu
    cp = ax.CurrentPoint;
    mouseX = cp(1,1);
    mouseY = cp(1,2);
    zone = handles.zone_popup.Value;
    % dodaje sipku armature na odgovarajuci polozaj u presjeku
    % rectangle object koji predstavlja sipku
    rebar = cs.addRebar(ax, ds, mouseX, mouseY, zone);
    set(rebar, 'ButtonDownFcn', {@RebarOnClick, handles});
    updateView(handles);
end

function updateView(handles)
% As_uitable iz ConcreteTool prozora
As_uitable = handles.As_uitable;
% As_uitable iz ovog prozora
As1_uitable = handles.As1_uitable;
% CrossSection objekat
cs = handles.crossSection;
% azurira As_uitable tabelu u ConcreteTool prozoru (As1 usvojeno)
As_uitable.Data(2) = cs.As1;
As_uitable.Data(4) = cs.As2;
% azurira As1_uitable u ovom prozoru
As1_uitable.Data(2:3,1) = [cs.As1 As1_uitable.Data(1)-cs.As1]';
As1_uitable.Data(2:3,2) = [cs.As2 As1_uitable.Data(4)-cs.As2]';
% azuriranje crteza armature u ConcreteTool prozoru
% ponovo crta sipke armature u PRVOM ConcreteTool prozoru
cs.drawRebar(handles.section1_axes);

% ponovo crta sipke armature u RebarTool prozoru
ax = handles.section_axes;
zone = handles.zone_popup.Value;
cs.plotGrid(ax, zone);
cs.drawRebar(ax)
% dodaje right click funkcionalnost nacrtanim objektima, da bi se
% ranije nacrtane sipke mogle brisati (nakon ponovnog otvaranja
% RebarTool prozora)
rebars = findobj(ax, 'Type', 'rectangle');
for i = 1:numel(rebars)
    set(rebars(i), 'ButtonDownFcn', {@RebarOnClick, handles});
end

% azuriranje rebar tabele
table = handles.rebar_uitable;
tableData = table.Data;
% Nedostaje As
if zone == 1
    As1_req = As1_uitable.Data(3);
else
    As1_req = As1_uitable.Data(6);
end
if As1_req < 0
    As1_req = 0;
end
% potreban broj sipki
tableData(:,4) = ceil(As1_req./tableData(:,2));
% As_uk
tableData(:,5) = tableData(:,2).*tableData(:,4);
table.Data = tableData;

function RebarOnClick(rect, ~, handles)
fig = handles.RebarTool_figure;
cs = handles.crossSection;
seltype = fig.SelectionType;
% provjera da li je napravljen desni klik
if strcmp(seltype,'alt')
    row = rect.UserData(1);
    column = rect.UserData(2);
    zone = rect.UserData(3);
    % ukloni Rebar objekat iz section.Rebars niza
    rebar = findobj(cs.Rebars, 'row', row, 'column', column, 'zone', zone);
    cs.Rebars(cs.Rebars == rebar) = [];
    % ukoliko je prikazan, obrisi rebar tag (balon koji prikazuje precnik
    % armature)
    tag = findobj(handles.section_axes, 'Tag', 'rebarTag');
    delete(tag);
    % ukloni rectangle i Rebar objekte koji su predstavljali tu sipku
    delete(rect);
    delete(rebar);
    % azuriraj prikaz
    updateView(handles);    
end


% --- Executes on button press in reset_pushbutton.
function reset_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to reset_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
cs = handles.crossSection;
zone = handles.zone_popup.Value;
rebars = findobj(cs.Rebars, 'zone', zone);
% uklanja sipke iz odgovarajuce zone
for i = 1:numel(rebars)
    cs.Rebars(cs.Rebars == rebars(i)) = [];
end
% ponovo iscrtava sve sipke
updateView(handles);



% --- Executes on selection change in zone_popup.
function zone_popup_Callback(hObject, eventdata, handles)
% hObject    handle to zone_popup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns zone_popup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from zone_popup
updateView(handles);
cs = handles.crossSection;
zone = handles.zone_popup.Value;
setDsPopup(handles);
ax = handles.section_axes;
if zone == 1
    ax.XLim = [0.65*cs.dims.h cs.dims.h];
    ax.YLim = [(cs.dims.bf-cs.dims.bw)/2 (cs.dims.bf-cs.dims.bw)/2+cs.dims.bw];
elseif zone == 2
    ax.XLim = [0 0.5*cs.dims.h];
    ax.YLim = [0 cs.dims.bf];
end


% --- Executes during object creation, after setting all properties.
function zone_popup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to zone_popup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function rebar_uitable_CreateFcn(hObject, eventdata, handles)
% hObject    handle to rebar_uitable (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on scroll wheel click while the figure is in focus.
function RebarTool_figure_WindowScrollWheelFcn(hObject, eventdata, handles)
% hObject    handle to RebarTool_figure (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.FIGURE)
%	VerticalScrollCount: signed integer indicating direction and number of clicks
%	VerticalScrollAmount: number of lines scrolled for each click
% handles    structure with handles and user data (see GUIDATA)

% podesavanje precnika armature pomocu tockica na misu
ds_popup = handles.ds_popup;
if eventdata.VerticalScrollCount > 0
    if ds_popup.Value < size(ds_popup.String, 1)
        ds_popup.Value = ds_popup.Value + 1;
    end
else
    if ds_popup.Value > 1
        ds_popup.Value = ds_popup.Value - 1;
    end
end


% --- Executes on mouse motion over figure - except title and menu.
function RebarTool_figure_WindowButtonMotionFcn(hObject, eventdata, handles)
% hObject    handle to RebarTool_figure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ax = handles.section_axes;
cs = handles.crossSection;
cp = ax.CurrentPoint;
tags = findobj(ax, 'Tag', 'rebarTag');
if ~isempty(tags)
   delete(tags);
end
for i = 1:numel(cs.Rebars)
    rebar = cs.Rebars(i);
    if rebar.mouseOver(cp)
        rebar.showTag(ax, cp);
    end
end

