function varargout = GUI(varargin)
% GUI MATLAB code for GUI.fig
%      GUI, by itself, creates a new GUI or raises the existing
%      singleton*.
%
%      H = GUI returns the handle to a new GUI or the handle to
%      the existing singleton*.
%
%      GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI.M with the given input arguments.
%
%      GUI('Property','Value',...) creates a new GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUI

% Last Modified by GUIDE v2.5 24-Jun-2014 15:33:42

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI_OutputFcn, ...
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


% --- Executes just before GUI is made visible.
function GUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUI (see VARARGIN)

% Choose default command line output for GUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);
global equation tspan beta gamma corr year;
%equation = 'relativistic';
set(handles.text2,'String','PAT^2');
% t0 = 0;
% tf = 1;
% steps = 365;
% tspan = linspace(t0,tf,steps);
% beta = 1;
% gamma = 1;
% corr = 1e5;
% year = 1*365*(23*60*60+56*60+4);
 
% UIWAIT makes GUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = GUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in body.
function body_Callback(hObject, eventdata, handles)
% hObject    handle to body (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns body contents as cell array
%        contents{get(hObject,'Value')} returns selected item from body


% --- Executes during object creation, after setting all properties.
function body_CreateFcn(hObject, eventdata, handles)
% hObject    handle to body (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in trajectory.
function trajectory_Callback(hObject, eventdata, handles)
% hObject    handle to trajectory (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in plot_oss.
%function plot_oss_Callback(hObject, eventdata, handles)
% hObject    handle to plot_oss (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

 

% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in equ_rel.
function equ_rel_Callback(hObject, eventdata, handles)
% hObject    handle to equ_rel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global equation;
% Hint: get(hObject,'Value') returns toggle state of equ_rel
%if get(hObject,'Value')
    equation = 'relativistic';
%end;


% --- Executes on button press in equ_newt.
function equ_newt_Callback(hObject, eventdata, handles)
% hObject    handle to equ_newt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global equation;
% Hint: get(hObject,'Value') returns toggle state of equ_newt
%if get(hObject,'Value')
    equation = 'newtonian';
%end;

% --- Executes on button press in run.
function run_Callback(hObject, eventdata, handles)
% hObject    handle to run (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Y;
global  beta gamma tol corr year runtime n pos vel actual_position mu;
t0=str2double(get(handles.init_time,'String'));
tf=str2double(get(handles.final_time,'String'));
steps=str2double(get(handles.num_steps,'String'));
%tspan = linspace(t0,tf,steps);
beta=str2double(get(handles.beta_value,'String'));
gamma=str2double(get(handles.gamma_value,'String'));
tol = str2double(get(handles.tolerance,'String'));
corr = str2double(get(handles.dist_corr,'String'));
year = 365*(23*60*60+56*60+4)*str2double(get(handles.time_corr,'String'));
n = 15;

[Y,runtime] = PAT2;
set(handles.out_runtime,'String',runtime);

% --- Executes on button press in pushbutton10.
function pushbutton10_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton11.
function pushbutton11_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in plot_iss.
function plot_iss_Callback(hObject, eventdata, handles)
% hObject    handle to plot_iss (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Y;
cla reset
set(handles.text2,'String','Displaying the outputed trajectory of the inner solar system');
axes(handles.traj_plot);
 Y_sun = Y(:,4:6);
 Y_mer = Y(:,7:9);
 Y_ven = Y(:,10:12);
 Y_ear = Y(:,13:15);
 Y_mar = Y(:,16:18);
 hold on;
 orange = [1 0.5 0];
 plot3(Y_sun(:,1),Y_sun(:,2),Y_sun(:,3),'LineStyle','o','Color',orange,'LineWidth',5);
 plot3(Y_mer(:,1),Y_mer(:,2),Y_mer(:,3),'g');
 plot3(Y_ven(:,1),Y_ven(:,2),Y_ven(:,3),'k');
 plot3(Y_ear(:,1),Y_ear(:,2),Y_ear(:,3),'b');
 plot3(Y_mar(:,1),Y_mar(:,2),Y_mar(:,3),'r');
 grid on;
 axis([-3e8 3e8 -3e8 3e8 -3e8 3e8]);
 xlabel('x [km]')
 ylabel('y [km]')
 zlabel('z [km]')
 hold off;
set(handles.text2,'String','Displaying the outputed trajectory of the inner solar system');

% --- Executes on button press in plot_asteroid.
function plot_asteroid_Callback(hObject, eventdata, handles)
% hObject    handle to plot_asteroid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Y;

set(handles.text2,'String','Displaying the outputed trajectory of the asteroid');
axes(handles.traj_plot);
cla reset
 Y_ast = Y(:,1:3);
 %Y_sun = Y(:,4:6);
 hold on;
 plot3(Y_ast(:,1),Y_ast(:,2),Y_ast(:,3),'m');
 %orange = [1 0.5 0];
 %plot3(Y_sun(:,1),Y_sun(:,2),Y_sun(:,3),'LineStyle','o','Color',orange,'LineWidth',4);
 xlabel('x [km]')
 ylabel('y [km]')
 zlabel('z [km]')
 axis([-3e8 3e8 -3e8 3e8 -3e8 3e8]);
 grid on;
 hold off;

% --- Executes on button press in plot_oss.
function plot_oss_Callback(hObject, eventdata, handles)
% hObject    handle to plot_oss (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Y;
axes(handles.traj_plot);
 cla reset
 Y_jup = Y(:,19:21);
 Y_sat = Y(:,22:24);
 Y_ura = Y(:,25:27);
 Y_nep = Y(:,28:30);
 Y_plu = Y(:,31:33);
 hold on;
 plot3(Y_jup(:,1),Y_jup(:,2),Y_jup(:,3),'g');
 plot3(Y_sat(:,1),Y_sat(:,2),Y_sat(:,3),'r');
 plot3(Y_ura(:,1),Y_ura(:,2),Y_ura(:,3),'k');
 plot3(Y_nep(:,1),Y_nep(:,2),Y_nep(:,3),'m');
 plot3(Y_plu(:,1),Y_plu(:,2),Y_plu(:,3),'b');
 grid on;
 xlabel('x [km]')
 ylabel('y [km]')
 zlabel('z [km]')
 hold off;
 axis([-2e9 2e9 -2e9 2e9 -2e9 2e9]);
set(handles.text2,'String','Displaying the outputed trajectory of the outer solar system');

% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function init_time_Callback(hObject, eventdata, handles)
% hObject    handle to init_time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of init_time as text
%        str2double(get(hObject,'String')) returns contents of init_time as a double


% --- Executes during object creation, after setting all properties.
function init_time_CreateFcn(hObject, eventdata, handles)
% hObject    handle to init_time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function final_time_Callback(hObject, eventdata, handles)
% hObject    handle to final_time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of final_time as text
%        str2double(get(hObject,'String')) returns contents of final_time as a double


% --- Executes during object creation, after setting all properties.
function final_time_CreateFcn(hObject, eventdata, handles)
% hObject    handle to final_time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function num_steps_Callback(hObject, eventdata, handles)
% hObject    handle to num_steps (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of num_steps as text
%        str2double(get(hObject,'String')) returns contents of num_steps as a double


% --- Executes during object creation, after setting all properties.
function num_steps_CreateFcn(hObject, eventdata, handles)
% hObject    handle to num_steps (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function beta_value_Callback(hObject, eventdata, handles)
% hObject    handle to beta_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of beta_value as text
%        str2double(get(hObject,'String')) returns contents of beta_value as a double


% --- Executes during object creation, after setting all properties.
function beta_value_CreateFcn(hObject, eventdata, handles)
% hObject    handle to beta_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function gamma_value_Callback(hObject, eventdata, handles)
% hObject    handle to gamma_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of gamma_value as text
%        str2double(get(hObject,'String')) returns contents of gamma_value as a double


% --- Executes during object creation, after setting all properties.
function gamma_value_CreateFcn(hObject, eventdata, handles)
% hObject    handle to gamma_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function time_corr_Callback(hObject, eventdata, handles)
% hObject    handle to time_corr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of time_corr as text
%        str2double(get(hObject,'String')) returns contents of time_corr as a double


% --- Executes during object creation, after setting all properties.
function time_corr_CreateFcn(hObject, eventdata, handles)
% hObject    handle to time_corr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function dist_corr_Callback(hObject, eventdata, handles)
% hObject    handle to dist_corr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dist_corr as text
%        str2double(get(hObject,'String')) returns contents of dist_corr as a double


% --- Executes during object creation, after setting all properties.
function dist_corr_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dist_corr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in energy_ratio.
function energy_ratio_Callback(hObject, eventdata, handles)
% hObject    handle to energy_ratio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global mu mass Y;
axes(handles.other_plot);
cla reset
n = 15; 
pos = Y(:,1:3*n);
vel = Y(:,3*n+1:6*n);

apo_pos = Y(:,1:3);

kinetic = zeros(length(pos),n);
potential = zeros(length(pos),n);
ratio = zeros(length(pos),n);

mu_sun = mu(2,1);

for i = 1:n
    for j = 1:length(apo_pos)
    if i == 2
        kinetic(:,i) = 0;
        potential(:,i) = 0;
        ratio(:,i) = 0;
        
    else
        mass_i = mass(i,1);
        kinetic(j,i) = 1/2*mass_i*norm(vel(j,3*(i-1)+1:3*(i-1)+3))^2;
        potential(j,i) = -mass_i*mu_sun/norm(pos(j,3*(i-1)+1:3*(i-1)+3));
        ratio(j,i) = -potential(j,i)/kinetic(j,i);

    end
    end
end
time = linspace(1,length(pos),length(pos));

plot(time,ratio(:,1))
xlabel('time [days]')
ylabel('Energy ration (V/T)')
set(gca,'FontSize',14)
axis([time(1) time(end) min(ratio(:,1)) max(ratio(:,1))])
set(handles.text2,'String','Displaying the ratio of potential over kinetic energy of the asteroid');



% --- Executes on button press in norm_error.
function norm_error_Callback(hObject, eventdata, handles)
% hObject    handle to norm_error (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Error vector is the same size as the output
global Y
n = 15;
pos = Y(:,1:3);
load apo;
actual_position = apo(:,1:3);
m = length(pos);
% Initialize error vector
error = zeros(m,1);

% Compute the error. The sign is changed when 
% position(i)<actual_position(i) so that the oscillatory behavior can be
% seen (if there is one)
for i = 1:length(pos)
%     if pos(i,:)>actual_position(i,:)
%         error(i) = norm(pos(i,:)-actual_position(i,:));
%     elseif pos(i)<actual_position(i)
        error(i,1) = -norm(pos(i,:)-actual_position(i,:));
%     end
end     

% Compute maximal error
%max_error_apo = max(abs(error));
% Compute error in the last time step
%last_error = abs(error(end));

% Plot the error (in norm)
axes(handles.other_plot);
plot(error);
set(gca,'FontSize',14)
%title('Numerical error in the position')
xlabel('Time [days]');
set(gca,'FontSize',14)
ylabel('Error [km]')

% global Y;
% axes(handles.other_plot);
% Y_ast = Y(:,1:3);
% load apo
% apo = apo(:,1:3);
% error = zeros(length(Y_ast),1);
% for i = 1:length(Y_ast)
% %     if apo-Y_ast<0
% %         error(:,1) = -norm(apo-Y_ast);
% %     else
%     error(i,1) = norm(apo(i,:)-Y_ast(i,:));
% %     end
% end
% plot(error)
% xlabel('Time [days]')
% ylabel('Error')
set(handles.text2,'String','Displaying the numerical error (in [km])');

% --- Executes on button press in error_component.
function error_component_Callback(hObject, eventdata, handles)
% hObject    handle to error_component (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Y;
axes(handles.other_plot);
Y_ast = Y(:,1:3);
load apo
apo = apo(:,1:3);
plot(Y_ast-apo);
xlabel('Time [days]')
ylabel('Error in each component')
legend('x','y','z')

% --- Executes on button press in energy_tot.
function energy_tot_Callback(hObject, eventdata, handles)
% hObject    handle to energy_tot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global mu mass Y;
[mu,mass,actual_position,y0,tspan] ...
          = setup2;  
%cla reset
axes(handles.other_plot);

n = 15; 
pos = Y(:,1:3*n);
vel = Y(:,3*n+1:6*n);

apo_pos = Y(:,1:3);

kinetic = zeros(length(pos),n);
potential = zeros(length(pos),n);
total = zeros(length(pos),n);

mu_sun = mu(2,1);

for i = 1:n
    for j = 1:length(apo_pos)
    if i == 2
        kinetic(:,i) = 0;
        potential(:,i) = 0;
        total(:,i) = kinetic(:,i) + potential(:,i); 
    else
        mass_i = mass(i,1);
        kinetic(j,i) = 1/2*mass_i*norm(vel(j,3*(i-1)+1:3*(i-1)+3))^2;
        potential(j,i) = -mass_i*mu_sun/norm(pos(j,3*(i-1)+1:3*(i-1)+3));
        total(j,i) = kinetic(j,i) + potential(j,i); 
    end
    end
end
time = linspace(1,length(pos),length(pos));

TOTAL = sum(total');
mean_total(:,1) = mean(total(:,1));
mean_TOTAL = mean(TOTAL);

TOTAL = TOTAL/mean_TOTAL;
total(:,1) = total(:,1)/mean_total(:,1);
 
plot(time,total(:,1))
axis([time(1) time(end) min(total(:,1)) max(total(:,1))])
xlabel('Time [days]')
ylabel('Normalized total Energy')
set(handles.text2,'String','Displaying the normalized total energy of the asteroid');


% --- Executes on button press in ast_and_inner.
function ast_and_inner_Callback(hObject, eventdata, handles)
% hObject    handle to ast_and_inner (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Y;

set(handles.text2,'String','Displaying the outputed trajectory of the asteroid');
axes(handles.traj_plot);
cla reset
 Y_ast = Y(:,1:3);
 %Y_sun = Y(:,4:6);
 Y_mer = Y(:,7:9);
 Y_ven = Y(:,10:12);
 Y_ear = Y(:,13:15);
 Y_mar = Y(:,16:18);
 hold on;
 %orange = [1 0.5 0];
 plot3(Y_ast(:,1),Y_ast(:,2),Y_ast(:,3),'m');
 %plot3(Y_sun(:,1),Y_sun(:,2),Y_sun(:,3),'LineStyle','o','Color',orange,'LineWidth',5);
 plot3(Y_mer(:,1),Y_mer(:,2),Y_mer(:,3),'g');
 plot3(Y_ven(:,1),Y_ven(:,2),Y_ven(:,3),'k');
 plot3(Y_ear(:,1),Y_ear(:,2),Y_ear(:,3),'b');
 plot3(Y_mar(:,1),Y_mar(:,2),Y_mar(:,3),'r');
 %orange = [1 0.5 0];
 %plot3(Y_sun(:,1),Y_sun(:,2),Y_sun(:,3),'LineStyle','o','Color',orange,'LineWidth',1);
 xlabel('x [km]')
 ylabel('y [km]')
 zlabel('z [km]')
 axis([-3e8 3e8 -3e8 3e8 -3e8 3e8]);
 grid on;
 hold off;



function tolerance_Callback(hObject, eventdata, handles)
% hObject    handle to tolerance (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tolerance as text
%        str2double(get(hObject,'String')) returns contents of tolerance as a double


% --- Executes during object creation, after setting all properties.
function tolerance_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tolerance (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in save_results.
function save_results_Callback(hObject, eventdata, handles)
% hObject    handle to save_results (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Y;
save Y;
load Y;