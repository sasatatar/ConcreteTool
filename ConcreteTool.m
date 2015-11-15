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

% Last Modified by GUIDE v2.5 04-Nov-2015 10:39:04

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
hObject.CloseRequestFcn = @close;
% podesavanje dimenzija prozora
hObject.Position(3:4) = [600 550];

% oznaka za ecu2 unutar axes objekta:
axis(handles.ecu2Label_axes, 'off');
text(0, 0.5, ['$\varepsilon_{cu2}$'],... % generise grcko epsilon sa cu2 u indeksu
    'HorizontalAlignment', 'Left',...
    'Interpreter', 'latex', 'FontSize', 12,...
    'Parent', handles.ecu2Label_axes);
text(0.6, 0.5, ['[' char(8240) ']:'],... % char(8240) generise simbol za promil
    'HorizontalAlignment', 'Left',...
    'FontSize', 11,...
    'Parent', handles.ecu2Label_axes);
% oznaka za alfa i gama_c unutar axes objekta:
axis(handles.alphaLabel_axes, 'off');
text(0, 0.5, ['$\alpha=$'],... % generise grcko epsilon sa cu2 u indeksu
    'HorizontalAlignment', 'Left',...
    'Interpreter', 'latex', 'FontSize', 11,...
    'Parent', handles.alphaLabel_axes);
text(0.8, 0.5, ['$\gamma_c=$'],... 
    'HorizontalAlignment', 'Left',...
    'Interpreter', 'latex', 'FontSize', 11,...
    'Parent', handles.alphaLabel_axes);

% create the CrossSection object and store it into handles structure
handles.crossSection = CrossSection;
handles.crossSection.plotSection(handles.crossSection_axes);


% set tabs
tabGroup = uitabgroup(hObject,'Position',[0 0 1 1], 'Tag', 'tabGroup');
sectionTab = uitab(tabGroup,'Title','Presjek', 'Tag', 'section_uitab');
rebarTab = uitab(tabGroup,'Title','Armatura', 'Tag', 'rebar_uitab');

handles.rebar_panel.Parent = rebarTab;
handles.rebar_panel.Position(1:2) = [15 15];
handles.section_panel.Parent = sectionTab;
handles.section_panel.Position(1:2) = [15 15];

% save tab handles
handles.tabGroup = tabGroup;
handles.section_uitab = sectionTab;
handles.rebar_uitab = rebarTab;
% set callback function
tabGroup.SelectionChangedFcn = {@TabCallback, handles};
% Update handles structure
guidata(hObject, handles);


function close(obj, data, handles)
delete(obj);

function TabCallback(hObject, callbackdata, handles)
% izvrsava se ako je odabran drugi tab (armatura)
if hObject.SelectedTab == handles.rebar_uitab
    % plot poprecnog presjeka i dilatacija
    updateRebarData(hObject, callbackdata, handles);
    dhRatio = str2double(handles.dhRatio_edit.String);
    cs = handles.crossSection;
    handles.d_edit.String = num2str(dhRatio*cs.dims.h, '%.1f');
    handles.As_uitable.Data = [cs.As1_req cs.As2_req; cs.As1 cs.As2];
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
function updateData(hobject, eventdata, handles)
cs = handles.crossSection;
% azurira karakteristicnu cvrstocu betona
concreteClass = [20:5:60 70:10:90;... % fck [MPa]
    3.5*ones(1,7) 3.1 2.9 2.7 2.6 2.6;... % ecu2 [promili]
    2.2 2.6 2.9 3.2 3.5 3.8 4.1 4.2 4.4 4.6 4.8 5]; % fctm {MPa] - cvrstoca na zatezanje
cs.fck = concreteClass(1,handles.concreteClass_popup.Value);
% azurira dopustenu dilataciju ecu2 u zavisnosti od klase, prema EC2
cs.ecu2 = -concreteClass(2,handles.concreteClass_popup.Value)/1000;
% azurira dopustenu cvrstocu na zatezanje
cs.fctm = concreteClass(3,handles.concreteClass_popup.Value);
handles.ecu2_edit.String = num2str(concreteClass(2,handles.concreteClass_popup.Value));
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
cs.alpha = str2double(handles.alpha_edit.String);
cs.gammac = str2double(handles.gammac_edit.String);
cs.Rebars = Rebar.empty;
cs.x = cs.xdRatio*0.9*cs.dims.h;
cs.As1_req = 0;
cs.As2_req = 0;
% provjera ako je bw = bf => u pitanju je pravougaoni presjek
if cs.dims.bw == cs.dims.bf
    cs.dims.hv = 0;
    cs.dims.hf = 0;
    handles.hf_edit.String = '0';
    handles.hv_edit.String = '0';
end
cs.plotSection(handles.crossSection_axes);

% ucitava podatke iz polja pod tabom Armatura i pohranjuje ih u crossSection objekat
function updateRebarData(hobject, eventdata, handles)
cs = handles.crossSection;
cs.plotSection(handles.section_axes);
cs.drawRebar(handles.section_axes);
cs.fyk = str2double(handles.fyk_edit.String);
cs.Es = str2double(handles.Es_edit.String);
% plotStrain funkciji se prosljedjuje i section_axes referenca da bi
% izracunala velicine dijagrama (da budu u ravnini)
cs.plotStrain(handles.strain_axes, handles.section_axes);
cs.plotCompression(handles.section_axes);


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


% --- Executes on button press in togglebutton1.
function togglebutton1_Callback(hObject, eventdata, handles)
% hObject    handle to togglebutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of togglebutton1


% --- Executes on button press in togglebutton2.
function togglebutton2_Callback(hObject, eventdata, handles)
% hObject    handle to togglebutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of togglebutton2



function alpha_edit_Callback(hObject, eventdata, handles)
% hObject    handle to alpha_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of alpha_edit as text
%        str2double(get(hObject,'String')) returns contents of alpha_edit as a double


% --- Executes during object creation, after setting all properties.
function alpha_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to alpha_edit (see GCBO)
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



function edit12_Callback(hObject, eventdata, handles)
% hObject    handle to alpha_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of alpha_edit as text
%        str2double(get(hObject,'String')) returns contents of alpha_edit as a double


% --- Executes during object creation, after setting all properties.
function edit12_CreateFcn(hObject, eventdata, handles)
% hObject    handle to alpha_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function gama_c_edit_Callback(hObject, eventdata, handles)
% hObject    handle to gammac_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of gammac_edit as text
%        str2double(get(hObject,'String')) returns contents of gammac_edit as a double


% --- Executes during object creation, after setting all properties.
function gama_c_edit_CreateFcn(hObject, eventdata, handles)
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
dhRatio = str2double(handles.dhRatio_edit.String);
cs = handles.crossSection;
handles.d_edit.String = num2str(dhRatio*cs.dims.h, '%.1f');

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
cs.plotStrain(handles.strain_axes, handles.section_axes);
cs.plotCompression(handles.section_axes);



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
if Msd~=0
    d = str2double(handles.d_edit.String);
    [As1 As2] = cs.calculateAs(Msd,d);
    cs.As1_req = As1;
    cs.As2_req = As2;
    handles.As_uitable.Data(1,1) = As1;
    handles.As_uitable.Data(1,2) = As2;
end
updateRebarData(hObject, 0, handles);



function N_edit_Callback(hObject, eventdata, handles)
% hObject    handle to N_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of N_edit as text
%        str2double(get(hObject,'String')) returns contents of N_edit as a double


% --- Executes during object creation, after setting all properties.
function N_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to N_edit (see GCBO)
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
