function varargout = ConcreteTool(varargin)
% CONCRETETOOL MATLAB code for ConcreteTool.fig
%      CONCRETETOOL, by itself, creates a new CONCRETETOOL or raises the existing
%      singleton*.
%
%      H = CONCRETETOOL returns the handle to a new CONCRETETOOL or the handle to
%      the existing singleton*.
%
%      CONCRETETOOL('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CONCRETETOOL.M with the given input arguments.
%
%      CONCRETETOOL('Property','Value',...) creates a new CONCRETETOOL or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ConcreteTool_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ConcreteTool_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ConcreteTool

% Last Modified by GUIDE v2.5 04-Feb-2016 20:45:25

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ConcreteTool_OpeningFcn, ...
                   'gui_OutputFcn',  @ConcreteTool_OutputFcn, ...
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


% --- Executes just before ConcreteTool is made visible.
function ConcreteTool_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ConcreteTool (see VARARGIN)

% Choose default command line output for ConcreteTool
handles.output = hObject;
% custom close function
hObject.CloseRequestFcn = @close;
% enable axes toolbar
set(hObject,'toolbar','figure');
% podesavanje dimenzija prozora
hObject.Position(3:4) = [707 550];

% generise oznake od grckih slova u interfejsu
drawLabels(handles);

% load the default CrossSection object from a template and store it into handles structure
section_template = fullfile(pwd,'section_templates','default_template.mat');
if exist(section_template, 'file')==2
    load(section_template);
else
    % if not available, create the default one
    section = CrossSection;
end

handles.crossSection = section;
% set tabs and save their handles
tabGroup = uitabgroup(hObject,'Position',[0 0 1 1], 'Tag', 'tabGroup');
handles.tabGroup = tabGroup;
handles.section_uitab = uitab(tabGroup,'Title','Presjek', 'Tag', 'section_uitab');
handles.rebar_uitab = uitab(tabGroup,'Title','Podužna armatura', 'Tag', 'rebar_uitab');
handles.stirrup_uitab = uitab(tabGroup, 'Title', 'Poprecna armatura', 'Tag', 'stirrup_uitab');
handles.torsion_uitab = uitab(tabGroup, 'Title', 'Proracun torzije', 'Tag', 'torsion_uitab');

handles.rebar_panel.Parent = handles.rebar_uitab;
handles.rebar_panel.Position(1:2) = [15 15];
handles.section_panel.Parent = handles.section_uitab;
handles.section_panel.Position(1:2) = [15 15];

%% ucitavanje programa za uzengije i torziju
handles = loadStirrup(hObject, handles);
handles = loadTorsionTool(hObject, handles);

% Update handles structure
guidata(hObject, handles);
handles.crossSection.plotSection(handles.crossSection_axes);

% set callback function for tab click
tabGroup.SelectionChangedFcn = {@TabCallback};
% load data to text fields
readData(handles);

function handles = loadStirrup(hObject, handles)
    %% ucitavanje programa za uzengije
    %% funkcija vraca azuriran handles
    StirrupTool;
    % get stirrup_figure handle
    % HandleVisibility mora da bude 'on', da bi findobj funkcija pronasla objekat !!!
    stirrup_figure = findobj('Type', 'figure', 'Tag', 'stirrup_figure');
    % get stirrup_panel handle
    Vsd_panel = findobj(stirrup_figure(1), 'Tag', 'Vsd_panel');
    % clear children in stirrup_uitab (just in case)
    delete(handles.stirrup_uitab.Children);
    % transfer stirrup_panel to stirrup_uitab
    Vsd_panel(1).Parent = handles.stirrup_uitab;
    % save handles from stirrup_figure
    handles.stirrup_tab = guidata(stirrup_figure(1));
    % brisanje prozora stirrup_figure / delete stirrup_figure
    delete(stirrup_figure);
    
    % podesavanje tabele za prikaz podataka
    table = handles.stirrup_tab.Vrd_table;
    %table.ColumnFormat = {'char', 'bank', 'bank', 'bank', 'bank'};
    tableData = {'Potr.' [] [] [] [] []; 'Usvo.' [] [] [] [] []};    
    table.Data = tableData;
    table.ColumnWidth = {50 62 62 62 62 64};
    
%% funkcija za ucitavanje programa za torziju
function handles = loadTorsionTool(hObject, handles)
    %% ucitavanje programa za torziju
    %% funkcija vraca azuriran handles
    TorsionTool;
    % get torsion_figure handle
    % HandleVisibility mora da bude 'on', da bi findobj funkcija pronasla objekat !!!
    torsion_figure = findobj('Type', 'figure', 'Tag', 'torsion_figure');
    % get Ted_panel handle
    Ted_panel = findobj(torsion_figure(1), 'Tag', 'Ted_panel');
    % clear children in torsion_uitab (just in case)
    delete(handles.torsion_uitab.Children);
    % transfer torsion_panel to torsion_uitab
    Ted_panel(1).Parent = handles.torsion_uitab;
    % save handles from torsion_figure
    handles.torsion_tab = guidata(torsion_figure(1));
    % brisanje prozora torsion_figure / delete torsion_figure
    delete(torsion_figure);
    
    % podesavanje tabele za prikaz podataka
    table = handles.torsion_tab.Trd_table;
    %table.ColumnFormat = {'char', 'bank', 'bank', 'bank', 'bank'};
    tableData = {'Potr.' [] [] [] [] [] []; 'Usvo.' [] [] [] [] [] []};    
    table.Data = tableData;
    table.ColumnWidth = {35 54 54 54 53 60 52};
    

    

function drawLabels(handles)
% oznaka za ecu2 unutar axes objekta:
axis(handles.ecu2Label_axes, 'off');
text(0, 0.5, ['$\varepsilon_{cu2}$'],... % generise grcko epsilon sa cu2 u indeksu
    'HorizontalAlignment', 'Left',...
    'Interpreter', 'latex', 'FontSize', 12,...
    'Parent', handles.ecu2Label_axes);
text(0.5, 0.5, ['[' char(8240) ']:'],... % char(8240) generise simbol za promil
    'HorizontalAlignment', 'Left',...
    'FontSize', 11,...
    'Parent', handles.ecu2Label_axes);
% oznaka za alfa i gama_c unutar axes objekta:
axis(handles.alphaLabel_axes, 'off'); % skriva ose
text(0, 0.5, ['$\alpha_{cc}=$'],... % generise grcko alfa sa cc u indeksu
    'HorizontalAlignment', 'Left',...
    'Interpreter', 'latex', 'FontSize', 11,...
    'Parent', handles.alphaLabel_axes);
text(0.7, 0.5, ['$\gamma_c=$'],... 
    'HorizontalAlignment', 'Left',...
    'Interpreter', 'latex', 'FontSize', 11,...
    'Parent', handles.alphaLabel_axes);

% oznaka za gamma_s unutar axes objekta:
axis(handles.gammas_axes, 'off');
text(0, 0.5, ['$\gamma_s=$'],... % generise grcko gama_s
    'HorizontalAlignment', 'Left',...
    'Interpreter', 'latex', 'FontSize', 12,...
    'Parent', handles.gammas_axes);


function close(obj, data, handles)
delete(obj);

function TabCallback(hObject, callbackdata)
handles = guidata(hObject);
cs = handles.crossSection;
% ucitava podatke iz CrossSection objekta i azurira prikaz u poljima
% poziva i updateAxes funkciju
readData(handles); 
if hObject.SelectedTab == handles.section_uitab
    
elseif hObject.SelectedTab == handles.rebar_uitab

elseif hObject.SelectedTab == handles.stirrup_uitab
    cs.calculateAsw([handles.stirrup_tab.Vsd_axes handles.stirrup_tab.al_axes]);
elseif hObject.SelectedTab == handles.torsion_uitab
    % provjera sjecnosti i ukljucivanje/iskljucivanje dugmica za s2
    cs = handles.crossSection;
    if cs.m == 2
        handles.torsion_tab.s2plus_button.Enable = 'off';
        handles.torsion_tab.s2minus_button.Enable = 'off';
    else
        handles.torsion_tab.s2plus_button.Enable = 'on';
        handles.torsion_tab.s2minus_button.Enable = 'on';
    end
end



% --- Outputs from this function are returned to the command line.
function varargout = ConcreteTool_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes during object creation, after setting all properties.
function concreteClass_popup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to concreteClass_popup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function c_nom_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to c_nom_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function dg_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dg_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function ecu2_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ecu2_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function bf_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to bf_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function hf_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to hf_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function hv_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to hv_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function bw_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to bw_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function h_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to h_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% ucitava podatke iz polja i pohranjuje ih u crossSection objekat
function updateData(hObject, eventdata, handles)
cs = handles.crossSection;
if handles.tabGroup.SelectedTab == handles.section_uitab
    % zastitni sloj c_nom
    cs.c_nom = str2double(handles.c_nom_edit.String);
    % maks. zrno agregata
    cs.dg = str2double(handles.dg_edit.String);
    % dimenzije
    cs.dims.bf = str2double(handles.bf_edit.String);
    cs.dims.hf = str2double(handles.hf_edit.String);
    cs.dims.hv = str2double(handles.hv_edit.String);
    cs.dims.bw = str2double(handles.bw_edit.String);
    cs.dims.h = str2double(handles.h_edit.String);
    cs.stirrup = str2double(handles.stirrup_edit.String);
    % provjera ako je bw >= bf => u pitanju je pravougaoni presjek
    % pa pomazemo korisniku tako sto automatski prepravljamo bw i ostale dimenzije
    if (cs.dims.bw >= cs.dims.bf && hObject==handles.bf_edit)
        cs.dims.bw = cs.dims.bf;
        handles.bw_edit.String = num2str(cs.dims.bw);
        cs.dims.hv = 0;
        cs.dims.hf = 0;
        handles.hf_edit.String = '0';
        handles.hv_edit.String = '0';
    end
    cs.x = cs.xdRatio*0.85*cs.dims.h;
    % reset d/h ratio:
    handles.dhRatio_edit.String = '0.9';
    handles.d_edit.String = num2str(0.9*cs.dims.h, '%.1f');
    
    cs.As1_req = 0;
    cs.As2_req = 0;
    
elseif handles.tabGroup.SelectedTab == handles.rebar_uitab
    % update rebar
    cs.fyk = str2double(handles.fyk_edit.String);
    cs.Es = str2double(handles.Es_edit.String);
    cs.gamma_s = str2double(handles.gammas_edit.String);
    cs.delta = str2double(handles.xdRatio_edit.String);
    cs.strainHardening = handles.strainHardening_checkbox.Value;
    % update concrete
    concreteClass = [12 16 20:5:60 70:10:90;... % fck [MPa]
        3.5*ones(1,9) 3.1 2.9 2.7 2.6 2.6]; % ecu2 [promili]
    cs.fck = concreteClass(1,handles.concreteClass_popup.Value);
    cs.alpha_cc = str2double(handles.alpha_cc_edit.String);
    cs.gammac = str2double(handles.gammac_edit.String);
    handles.ecu2_edit.String = num2str(concreteClass(2,handles.concreteClass_popup.Value));
    if cs.fck <= 50
        xdRatio = 0.45;
    else
        xdRatio = 0.35;
    end
    if hObject == handles.xdRatio_edit
        xdRatio = str2double(handles.xdRatio_edit.String);
    else
        handles.xdRatio_edit.String = num2str(xdRatio);
    end
    cs.xdRatio = xdRatio;
    % axial force
    cs.Nsd = str2double(handles.Nsd_edit.String)*1000; % [N]
    % staticka visina d - sinhronizacija polja
    if hObject == handles.dhRatio_edit
        dhRatio = str2double(handles.dhRatio_edit.String);
        handles.d_edit.String = num2str(dhRatio*cs.dims.h, '%.1f');
    elseif hObject == handles.d_edit
        d = str2double(handles.d_edit.String);
        handles.dhRatio_edit.String = num2str(d/cs.dims.h, '%.2f');
    end
elseif handles.tabGroup.SelectedTab == handles.stirrup_uitab
    fields = {'fywk', 'stirrup', 'm', 'alpha', 'sl'};
    for i = 1:numel(fields)
        field = fields{i};
        eval(['cs.' field ' = str2double(handles.stirrup_tab.' field '_edit.String);']);
    end
    % Vsd ne moze u for petlju jer se konvertuje u N
    cs.Ved = str2double(handles.stirrup_tab.Vsd_edit.String)*1000;
elseif handles.tabGroup.SelectedTab == handles.torsion_uitab
    fields = {'fywk', 'stirrup', 'm'};
    for i = 1:numel(fields)
        field = fields{i};
        eval(['cs.' field ' = str2double(handles.torsion_tab.' field '_edit.String);']);
    end
    % Ted ne moze u for petlju jer se konvertuje u Nmm
    cs.Ted = str2double(handles.torsion_tab.Ted_edit.String)*10^6; % Nmm
    cs.Ved = str2double(handles.torsion_tab.Ved_edit.String)*10^3; % N
    
    % provjera sjecnosti i ukljucivanje/iskljucivanje dugmica za s2
    cs = handles.crossSection;
    if cs.m == 2
        handles.torsion_tab.s2plus_button.Enable = 'off';
        handles.torsion_tab.s2minus_button.Enable = 'off';
    else
        handles.torsion_tab.s2plus_button.Enable = 'on';
        handles.torsion_tab.s2minus_button.Enable = 'on';
    end
end
updateAxes(handles);

% azurira graficke prikaze i tabelu armature
function updateAxes(handles)
cs = handles.crossSection;
if handles.tabGroup.SelectedTab == handles.section_uitab
    cs.plotSection(handles.crossSection_axes);
elseif handles.tabGroup.SelectedTab == handles.rebar_uitab
    cs.plotSection(handles.section_axes);
    cs.drawRebar(handles.section_axes);
    handles.As_uitable.Data = [cs.As1_req cs.As2_req; cs.As1 cs.As2];
    
    % plotStrain funkciji se prosljedjuje i section_axes referenca da bi
    % izracunala velicine dijagrama (da budu u ravnini)
    cs.plotStrain(handles.strain_axes, handles.section_axes);
    cs.plotCompression(handles.section_axes);
    cs.plotStress(handles.stress_axes, handles.section_axes);
end


function readData(handles)
cs = handles.crossSection;
if handles.tabGroup.SelectedTab == handles.rebar_uitab
    handles.fyk_edit.String = num2str(cs.fyk);
    handles.Es_edit.String = num2str(cs.Es);
    handles.gammas_edit.String = num2str(cs.gamma_s);
    handles.xdRatio_edit.String = num2str(cs.xdRatio);
    
    handles.Nsd_edit.String = num2str(cs.Nsd/1000); %[kN]
    handles.Msd_edit.String = num2str(cs.Msd*10^-6); %[kNm]
    
    handles.d_edit.String = num2str(0.9*cs.dims.h, '%.2f');
    
    % azurira karakteristicnu cvrstocu betona
    concreteClass = [12 16 20:5:60 70:10:90]; % fck [MPa]
    handles.concreteClass_popup.Value = find(concreteClass==cs.fck);
    handles.alpha_cc_edit.String = num2str(cs.alpha_cc);
    handles.gammac_edit.String = num2str(cs.gammac);
    % load section data
    handles.ecu2_edit.String = num2str(cs.ecu2*-1000);
elseif handles.tabGroup.SelectedTab == handles.section_uitab
    % zastitni sloj c_nom
    handles.c_nom_edit.String = num2str(cs.c_nom);
    % maks. zrno agregata
    handles.dg_edit.String = num2str(cs.dg);
    % stirrups / uzengije
    handles.stirrup_edit.String = num2str(cs.stirrup);
    % dimenzije
    handles.bf_edit.String = num2str(cs.dims.bf);
    handles.hf_edit.String = num2str(cs.dims.hf);
    handles.hv_edit.String = num2str(cs.dims.hv);
    handles.bw_edit.String = num2str(cs.dims.bw);
    handles.h_edit.String = num2str(cs.dims.h);
elseif handles.tabGroup.SelectedTab == handles.stirrup_uitab
    % ucitavanje podataka za poprecnu armaturu
    fields = {'fywk', 'stirrup', 'm', 'alpha', 'sl'};
    for i = 1:numel(fields)
        field = fields{i};
        eval(['handles.stirrup_tab.' field '_edit.String = num2str(cs.' field ');']);
    end
    handles.stirrup_tab.Vsd_edit.String = num2str(cs.Ved/1000);
elseif handles.tabGroup.SelectedTab == handles.torsion_uitab
    % ucitavanje podataka za proracun torzije
    fields = {'fywk', 'stirrup', 'm'};
    for i = 1:numel(fields)
        field = fields{i};
        eval(['handles.torsion_tab.' field '_edit.String = num2str(cs.' field ');']);
    end
    handles.torsion_tab.Ted_edit.String = num2str(cs.Ted*10^-6);
    handles.torsion_tab.Ved_edit.String = num2str(cs.Ved*10^-3);
end
updateAxes(handles);



% --- Executes during object creation, after setting all properties.
function section_axes_CreateFcn(hObject, eventdata, handles)
% hObject    handle to section_axes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called



function stirrup_edit_Callback(hObject, eventdata, handles)
% hObject    handle to stirrup_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of stirrup_edit as text
%        str2double(get(hObject,'String')) returns contents of stirrup_edit as a double


% --- Executes during object creation, after setting all properties.
function stirrup_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to stirrup_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function alpha_cc_edit_Callback(hObject, eventdata, handles)
% hObject    handle to alpha_cc_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of alpha_cc_edit as text
%        str2double(get(hObject,'String')) returns contents of alpha_cc_edit as a double


% --- Executes during object creation, after setting all properties.
function alpha_cc_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to alpha_cc_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function gammac_edit_Callback(hObject, eventdata, handles)
% hObject    handle to gammac_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of gammac_edit as text
%        str2double(get(hObject,'String')) returns contents of gammac_edit as a double


% --- Executes during object creation, after setting all properties.
function gammac_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to gammac_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function fyk_edit_Callback(hObject, eventdata, handles)
% hObject    handle to fyk_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fyk_edit as text
%        str2double(get(hObject,'String')) returns contents of fyk_edit as a double


% --- Executes during object creation, after setting all properties.
function fyk_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fyk_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Es_edit_Callback(hObject, eventdata, handles)
% hObject    handle to Es_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Es_edit as text
%        str2double(get(hObject,'String')) returns contents of Es_edit as a double


% --- Executes during object creation, after setting all properties.
function Es_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Es_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function As1_edit_Callback(hObject, eventdata, handles)
% hObject    handle to As1_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of As1_edit as text
%        str2double(get(hObject,'String')) returns contents of As1_edit as a double


% --- Executes during object creation, after setting all properties.
function As1_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to As1_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function As2_edit_Callback(hObject, eventdata, handles)
% hObject    handle to As2_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of As2_edit as text
%        str2double(get(hObject,'String')) returns contents of As2_edit as a double

% --- Executes during object creation, after setting all properties.
function As2_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to As2_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in debug_pushbutton.
function debug_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to debug_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
keyboard

% --- Executes on key press with focus on debug_pushbutton and none of its controls.
function debug_pushbutton_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to debug_pushbutton (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in addAs_pushbutton.
function addAs_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to addAs_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
RebarTool(handles);

function dhRatio_edit_Callback(hObject, eventdata, handles)
% hObject    handle to dhRatio_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dhRatio_edit as text
%        str2double(get(hObject,'String')) returns contents of dhRatio_edit as a double


% --- Executes during object creation, after setting all properties.
function dhRatio_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dhRatio_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Msd_edit_Callback(hObject, eventdata, handles)
% hObject    handle to Msd_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Msd_edit as text
%        str2double(get(hObject,'String')) returns contents of Msd_edit as a double


% --- Executes during object creation, after setting all properties.
function Msd_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Msd_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Mrd_edit_Callback(hObject, eventdata, handles)
% hObject    handle to Mrd_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Mrd_edit as text
%        str2double(get(hObject,'String')) returns contents of Mrd_edit as a double


% --- Executes during object creation, after setting all properties.
function Mrd_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Mrd_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function As_uitable_CreateFcn(hObject, eventdata, handles)
% hObject    handle to As_uitable (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
hObject.Position(3) = 220;
hObject.Position(4) = 64;
hObject.Data = zeros(2,2);
hObject.RowName = {'Potr.','Usvo.'};
hObject.ColumnName = {'As1 [mm2]', 'As2 [mm2]'};
hObject.ColumnFormat = {'bank', 'bank'};


% --- Executes on button press in addAs2_pushbutton.
function addAs2_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to addAs2_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in calculateMrd_pushbutton.
function calculateMrd_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to calculateMrd_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
cs = handles.crossSection;
handles.Mrd_edit.String = num2str(cs.Mrd*10^-6,'%.2f'); % [kNm]
% plotStrain funkciji se prosljedjuje i section_axes referenca da bi
% izracunala velicine dijagrama (da budu u ravnini)
% azuriraj prikaz
updateAxes(handles);


function d_edit_Callback(hObject, eventdata, handles)
% hObject    handle to d_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of d_edit as text
%        str2double(get(hObject,'String')) returns contents of d_edit as a double
d = str2double(handles.d_edit.String);
cs = handles.crossSection;
handles.dhRatio_edit.String = num2str(d/cs.dims.h, '%.2f');

% --- Executes during object creation, after setting all properties.
function d_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to d_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on button press in calculateAs_pushbutton.
function calculateAs_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to calculateAs_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
cs = handles.crossSection;
Msd = str2double(handles.Msd_edit.String)*10^6; % [Nmm]
if ~isnan(Msd) && Msd > 0
    cs.Msd = Msd;
    d = str2double(handles.d_edit.String);
    [As1 As2] = cs.calculateAs(d);
    cs.As1_req = As1;
    cs.As2_req = As2;
else
    As1 = 0;
    As2 = 0;
    cs.As1_req = As1;
    cs.As2_req = As2;
end
% azurira grafikone i tabele
updateAxes(handles);



function Nsd_edit_Callback(hObject, eventdata, handles)
% hObject    handle to Nsd_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Nsd_edit as text
%        str2double(get(hObject,'String')) returns contents of Nsd_edit as a double


% --- Executes during object creation, after setting all properties.
function Nsd_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Nsd_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function alphaLabel_axes_CreateFcn(hObject, eventdata, handles)
% hObject    handle to alphaLabel_axes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate alphaLabel_axes



% --- Executes during object creation, after setting all properties.
function gammas_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to gammas_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function xdRatio_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to xdRatio_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
tooltip = ['Maksimalan dopusten odnos x /d\n'...
    '(x - polozaj neutralne ose, d - staticka visina presjeka).'];
hObject.TooltipString = sprintf(tooltip);


% --- Executes during object creation, after setting all properties.
function gammas_axes_CreateFcn(hObject, eventdata, handles)
% hObject    handle to gammas_axes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate gammas_axes



function gammas_edit_Callback(hObject, eventdata, handles)
% hObject    handle to gammas_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of gammas_edit as text
%        str2double(get(hObject,'String')) returns contents of gammas_edit as a double




% --- Executes on selection change in concreteClass_popup.
function concreteClass_popup_Callback(hObject, eventdata, handles)
% hObject    handle to concreteClass_popup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns concreteClass_popup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from concreteClass_popup



function ecu2_edit_Callback(hObject, eventdata, handles)
% hObject    handle to ecu2_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ecu2_edit as text
%        str2double(get(hObject,'String')) returns contents of ecu2_edit as a double


% --- Executes on button press in open_button.
function open_button_Callback(hObject, eventdata, handles)
% hObject    handle to open_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% load the default CrossSection object from a template and store it into handles structure
path = fullfile(pwd,'section_templates');
file = uigetfile(fullfile(path,'*.mat;'),'Select file');
if file==0
    return;
end
load(fullfile(path, file));
handles.crossSection = section;
% Update handles structure
guidata(hObject, handles);
handles.crossSection.plotSection(handles.crossSection_axes);
% ucitavanje podataka u tekstualna polja u interfejsu
readData(handles);


% --- Executes on button press in save_button.
function save_button_Callback(hObject, eventdata, handles)
% hObject    handle to save_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
section = handles.crossSection;
file = uiputfile(fullfile(pwd,'section_templates','*.mat;'),'Save file');
% check if file is selected
if file ~= 0
    save(fullfile(pwd,'section_templates',file), 'section');
end


% --- Executes during object creation, after setting all properties.
function delta_axes_CreateFcn(hObject, eventdata, handles)
% hObject    handle to delta_axes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate delta_axes



function delta_edit_Callback(hObject, eventdata, handles)
% hObject    handle to delta_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of delta_edit as text
%        str2double(get(hObject,'String')) returns contents of delta_edit as a double


% --- Executes during object creation, after setting all properties.
function delta_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to delta_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in strainHardening_checkbox.
function strainHardening_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to strainHardening_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of strainHardening_checkbox


% --- Executes on button press in Mfi_button.
function Mfi_button_Callback(hObject, eventdata, handles)
% hObject    handle to Mfi_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
cs = handles.crossSection;
cs.plotMfi();
