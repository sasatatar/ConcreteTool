function varargout = StirrupTool(varargin)
% STIRRUPTOOL MATLAB code for StirrupTool.fig
%      STIRRUPTOOL, by itself, creates a new STIRRUPTOOL or raises the existing
%      singleton*.
%
%      H = STIRRUPTOOL returns the handle to a new STIRRUPTOOL or the handle to
%      the existing singleton*.
%
%      STIRRUPTOOL('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in STIRRUPTOOL.M with the given input arguments.
%
%      STIRRUPTOOL('Property','Value',...) creates a new STIRRUPTOOL or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before StirrupTool_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to StirrupTool_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help StirrupTool

% Last Modified by GUIDE v2.5 01-May-2016 17:50:03

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @StirrupTool_OpeningFcn, ...
                   'gui_OutputFcn',  @StirrupTool_OutputFcn, ...
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


% --- Executes just before StirrupTool is made visible.
function StirrupTool_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to StirrupTool (see VARARGIN)

% Choose default command line output for StirrupTool
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes StirrupTool wait for user response (see UIRESUME)
% uiwait(handles.stirrup_figure);


% --- Outputs from this function are returned to the command line.
function varargout = StirrupTool_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in debug_button.
function debug_button_Callback(hObject, eventdata, handles)
% hObject    handle to debug_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%fprintf('Es: %s\n', handles.Es_edit.String);
keyboard



function Vsd_edit_Callback(hObject, eventdata, handles)
% hObject    handle to Vsd_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Vsd_edit as text
%        str2double(get(hObject,'String')) returns contents of Vsd_edit as a double


% --- Executes during object creation, after setting all properties.
function Vsd_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Vsd_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



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



function m_edit_Callback(hObject, eventdata, handles)
% hObject    handle to m_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of m_edit as text
%        str2double(get(hObject,'String')) returns contents of m_edit as a double


% --- Executes during object creation, after setting all properties.
function m_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to m_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function sl_edit_Callback(hObject, eventdata, handles)
% hObject    handle to sl_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of sl_edit as text
%        str2double(get(hObject,'String')) returns contents of sl_edit as a double


% --- Executes during object creation, after setting all properties.
function sl_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sl_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function fywk_edit_Callback(hObject, eventdata, handles)
% hObject    handle to fywk_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fywk_edit as text
%        str2double(get(hObject,'String')) returns contents of fywk_edit as a double


% --- Executes during object creation, after setting all properties.
function fywk_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fywk_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function gammas2_edit_Callback(hObject, eventdata, handles)
% hObject    handle to gammas2_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of gammas2_edit as text
%        str2double(get(hObject,'String')) returns contents of gammas2_edit as a double


% --- Executes during object creation, after setting all properties.
function gammas2_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to gammas2_edit (see GCBO)
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


% --- Executes on button press in calculateAsw.
function calculateAsw_Callback(hObject, eventdata, handles)
% hObject    handle to calculateAsw (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
cs = handles.crossSection;
[out, ~] = cs.calculateAsw([handles.stirrup_tab.Vsd_axes handles.stirrup_tab.al_axes]);
% convert output to cell array and round to two decimal places
out = arrayfun(@(x)(sprintf('%.2f',x)), out, 'unif', 0);
table = handles.stirrup_tab.Vrd_table;
table.Data(1,2:6) = out(1:end);


% --- Executes on button press in calculateVrd.
function calculateVrd_Callback(hObject, eventdata, handles)
% hObject    handle to calculateVrd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
cs = handles.crossSection;
[~, out] = cs.calculateAsw([handles.stirrup_tab.Vsd_axes handles.stirrup_tab.al_axes]);
% convert output to cell array and round to two decimal places
out = arrayfun(@(x)(sprintf('%.2f',x)), out, 'unif', 0);
table = handles.stirrup_tab.Vrd_table;
table.Data(2,2:6) = out(1:end);
