%% ENM 056 : Machine Design Assignment
% **************************************************************************
% This is a template file to help the students start with machine design 
% assignment in ENM 056.
% Uncomment the code one section at a time and fill in the relevant details.
% In order to run the code section by section, use "Run Section" option in
% MATLAB or press Ctrl + Enter.
% How to enable code folding to work with long piece of code
% Home -> Preference -> Editor/Debugger -> Code Folding -> Enable if/else blocks and Sections
% How to enable Datatips to view the variable values directly on code
% instead of workspace
% Home -> Preference -> Editor/Debugger -> Display -> Enable datatips in edit mode
% Note: Please be careful when changing variable names in the template
% file. Also, the students may still need to write some sections by
% themselves.
% **************************************************************************

clear
close all
clc

% Parameters of the reference machine
%Please fill the parameters of the machine below in SI units

mm = 1e-3;            % mm to SI unit
OD_stator = 184*mm;     % Outer diameter of stator 
ID_stator = 130*mm;     % Inner diameter of stator
OD_rotor = 128*mm;          % Outer diameter of rotor
ID_rotor = 60*mm;          % Inner diameter of rotor
Hs0 = 0.5*mm;               % Slot opening height
Hs1 = 0.5*mm;               % Slot wedge height
Hs2 = 14*mm;               % Slot body height
w_tooth = 4.4*mm;           % Tooth width
Rs = 0.5*mm;                % Slot bottom radius fillet
Bs0 = 2*mm;               % Slot opening
N_pole = 8;            % Number of poles
N_slot = 48;            % Number of slots
N_parallel = 2;         % Number of parallel branch
r = 2;                 % Number of winding layers
y = 5;                 % Coil pitch
d_strand = 0.813*mm;          % Wire diameter
f_Cu_max = 45;          % Maximum Cu-fill factor
mu_0 = 4 * pi * 1e-7; % Magnetic permeability of vacuum [H/m]
rho_Cu = 1.72e-8;     % Resistivity of Cu [ohm/m]
t_mag = 5*mm;       % Thickness of magnet, uncomment when necessary
w_mag = 15*mm;        % Width of magnet, uncomment when necessary
L_stack = 115*mm;           % Stack length
V_DC = 650;              % DC link voltage [V]
d_sy_1 = ID_stator+2*(Hs0+Hs1+Hs2+Rs);
d_sy_2 = OD_stator;
d_sy_m = (d_sy_1+OD_stator)/2;

%% Load line
% 
% % Start by assuming an air-gap flux density [T]
B_gap = 0.01: 0.001: 1;
%B_gap = 0.7;
% 
% % The effective air-gap cross section perpendicular to flux crossing the air-gap [sq.m]
A_gap = pi/16*(ID_stator+OD_stator)/2*L_stack; 
% 
% % Flux in the air-gap [Wb]
Phi_gap = A_gap*B_gap;
% 
% % Cross-section of flux path in stator tooth. Assuming that the flux passes through 2.5 stator teeth [sq.m]
A_tooth = 2.5*w_tooth*L_stack;
% 
% % Cross-section of stator yoke perpendicular to flux [sq.m]
A_yoke = (d_sy_2 - d_sy_1)/2;
% 
% % Cross section of flux path in rotor assuming that the flux is concentrated in a tube with width equal to magnet width [sq.m]
A_rotor = w_mag*L_stack;
% 
% % Output cross sections of different parts
disp('Cross-section of different parts')
fprintf('Air-gap cross-section = % .2f [mm^2] \n',A_gap * 1e6)
fprintf('Stator tooth cross-section = % .2f [mm^2] \n',A_tooth * 1e6)
fprintf('Stator yoke cross-section = % .2f [mm^2] \n \n',A_yoke * 1e6)
fprintf('Rotor cross-section = % .2f [mm^2] \n \n',A_rotor * 1e6)

% % Flux densities in different parts of the machine [T]
B_tooth = Phi_gap/A_tooth; % Flux density in stator tooth
B_yoke =  Phi_gap/A_yoke; % Flux density in stator yoke 
B_rotor = Phi_gap/A_rotor; % Flux density in rotor 

index = B_gap == 0.7;

% Output corresponding flux densities when air gap flux density is 0.7 [T]
fprintf('Flux densities in the different part of the machine when air-gap flux density is % .2f [T]\n', B_gap(index))
fprintf('Stator tooth flux density = % .2f [T] \n',B_tooth(index))
fprintf('Stator yoke flux density = % .2f [T] \n',B_yoke(index))
fprintf('Rotor yoke flux density = % .2f [T] \n \n',B_rotor(index))

%Import the B-H curve of the M235-50A steel from a TAB file
BH_data = importdata('SURA M250-35A - BH Curve @ 50 Hz.tab'); % import data
H_data = BH_data(:,1);           % copy (Row All , Column One) as H
B_data = BH_data(:,2);   % copy (Row All , Column Two) as B

figure(10)
plot(H_data, B_data)
xlabel('H')
ylabel('B')
grid on

%%
% Calculated magnetic field intensity [A/m]
% Interpolation
method   = 'linear';     % 'linear' or 'spline' can be selected as interpolation method
H_stator_tooth = interp1(B_data,H_data,B_tooth,method);   % Interpolate stator tooth flux density, B_tooth to calculate corresponding H_tooth
H_stator_yoke = interp1(B_data, H_data, B_yoke, method);       % Interpolate stator yoke flux density, B_yoke to calculate corresponding H_yoke
H_rotor = interp1(B_data,H_data, B_rotor, method);              % Interpolate rotor flux density, B_rotor to calculate corresponding H_rotor
H_gap = B_gap/mu_0;                % Calculate magnetic field intensity in the air-gap. Note: Air-gap does not contain iron.

% Output magnetic field intensities when the air gap flux density is 0.7 T
fprintf('Magnetic field intensities in the different part of the machine when air-gap flux density is % .2f [T] \n', B_gap(index))
fprintf('Stator tooth Magnetic field intensity = % .2f [A/m] \n',H_stator_tooth(index))
fprintf('Stator yoke Magnetic field intensity = % .2f [A/m] \n',H_stator_yoke(index))
fprintf('Rotor Magnetic field intensity = % .2f [A/m] \n \n',H_rotor(index))
fprintf('Air-gap Magnetic field intensity = % .2f [A/m] \n',H_gap(index))
% 
% Length of flux path [m]
l_stator_tooth = Hs0+Hs1+Hs2+Rs;  % Length of flux path in stator tooth
l_stator_yoke = 4/48*pi*d_sy_m+(d_sy_2-d_sy_1)/2;   % Length of flux path in stator yoke
l_rotor = pi/2*(4/48*pi*OD_rotor)-2*t_mag;         % Length of flux path in rotor
l_gap = (ID_stator-OD_rotor)/2;           % Length of flux path in air-gap
% 
% Output length of flux path in different parts of the machine
% disp('Length of flux path in different parts of the machine')
fprintf('Length of flux path in stator tooth = % .2f [mm] \n',l_stator_tooth * 1e3)
fprintf('Length of flux path in stator yoke = % .2f [mm] \n',l_stator_yoke * 1e3)
fprintf('Length of flux path in rotor = % .2f [mm] \n',l_rotor * 1e3)
fprintf('Length of flux path in air-gap = % .2f [mm] \n \n',l_gap * 1e3)
% 
% MMF drop in different parts [Aturn]
MMF_stator_tooth = l_stator_tooth*H_stator_tooth;     % MMF drop in stator tooth
MMF_stator_yoke = l_stator_yoke*H_stator_yoke;      % MMF drop in stator yoke
MMF_rotor = l_rotor*H_rotor;            % MMF drop in rotor
MMF_gap = l_gap*H_gap;              % MMF drop in air-gap
% 
% MMF drop in different parts when the air gap flux density is 0.7T 
fprintf('MMF drop in different parts of the machine when air-gap flux density is % .2f [T] \n', B_gap(index))
fprintf('MMF drop in stator tooth = % .2f [A-turn] \n',MMF_stator_tooth(index))
fprintf('MMF drop in stator yoke = % .2f [A-turn] \n',MMF_stator_yoke(index))
fprintf('MMF drop in rotor = % .2f [A-turn] \n',MMF_rotor(index))
fprintf('MMF drop in air-gap = % .2f [A-turn] \n\n',MMF_gap(index))
% 
% Total MMF drop
MMF_total = 2*MMF_gap+MMF_rotor+MMF_stator_yoke+2*MMF_stator_tooth;
% 
% MMF drop in iron and air-part
MMF_iron = MMF_total-2*MMF_gap;           % MMF drop in the iron parts
MMF_air = MMF_gap*2;            % MMF drop in air-gap
% 
% Plot load line
figure(1)               % create Figure 1
clf                     % clear figure
subplot(1,3,1)
plot(MMF_total,Phi_gap*1e3, 'LineWidth', 2)
xlabel('Total MMF drop [A-turn]')
ylabel('Flux in air-gap[mWb]')
legend('Load line')
grid on

subplot(1,3,2)
hold on
plot(MMF_iron,Phi_gap*1e3, 'LineWidth', 2)
plot(MMF_air,Phi_gap*1e3, 'LineWidth', 2)
xlabel('MMF drop [A-turn]')
ylabel('Flux in air-gap[mWb]')
legend('MMF drop in iron','MMF drop in air')
grid on

subplot(1,3,3)
hold on
plot(MMF_air,Phi_gap*1e3, 'LineWidth', 2)
xlabel('MMF drop [A-turn]')
ylabel('Flux in air-gap[mWb]')
legend('MMF drop in air')
grid on

%Source Line

B_mag = [0, 0.5912, 1.1824]; % B data of magnet [T]
H_mag = [-902285, -451142, 0]; % H data of magnet [A/m]

MMF_mag = 2*H_mag*t_mag;
Phi_mag = 2*B_mag*w_mag*L_stack;

%Magnet characteristic
figure(2)
clf
plot(MMF_mag,Phi_mag * 1e3, 'LineWidth', 2)
xlim([min(MMF_mag) 0])
xlabel('MMF produced by magnet [A-turn]')
ylabel('Flux due to magnet [mWb]')
legend('Demagnetization characteristic')
grid on

%Plot load line and sorce line together
figure(3)
clf
plot(MMF_mag,Phi_mag * 1e3, 'LineWidth', 2)
hold on
plot(-MMF_total,Phi_gap * 1e3, 'LineWidth', 2)
hold off
xlim([min(MMF_mag) 0])
xlabel('MMF [A-turn]')
ylabel('Flux [mWb]')
legend('Demagnetization characteristic' , 'Load line')
grid on

%Plot the B-H curve of the core material
figure(4)

%Find the intersection
%To avoid extrapolation of data, we can limit the x-axis of the load line to maximum MMF that can be produced by magnet
index = MMF_total < max(abs(MMF_mag));
MMF_total_con = MMF_total(index);
Phi_gap_con = Phi_gap(index);

%Interpolate magnet characteristic corresponding to total MMF drop to have same x-axis. Don't forget the negative sine to move it to second quadrant
Phi_mag_interp = interp1( abs(MMF_mag), Phi_mag, MMF_total_con, method);
[value, index] = min(abs(Phi_mag_interp - Phi_gap_con));

Phi_gap_no_load = Phi_gap_con(index);          % Flux at the noload operating point
B_gap_no_load = Phi_gap_no_load/A_gap;       % Flux density in airgap
B_tooth_no_load = Phi_gap_no_load/A_tooth;     % Flux density in stator tooth
B_yoke_no_load = Phi_gap_no_load/A_yoke;      % Flux density in stator yoke
B_rotor_no_load = Phi_gap_no_load/A_rotor;     % Flux density in rotor
% 
fprintf('Open cicruit MMF = % .2f [A] \n',min(MMF_mag))
fprintf('Short circuit flux = % .2f [mWb] \n',max(Phi_mag)*1e3)
fprintf('Magnet thickness = % .2f [mm] \n',t_mag * 1e3)
fprintf('Magnet width = % .2f [mm] \n',w_mag * 1e3)
fprintf('No load flux in the air-gap = % .5f [Wb] \n',Phi_gap_no_load)
fprintf('No load flux density in the air-gap = % .2f [T] \n',B_gap_no_load)
fprintf('No load flux density in stator tooth = % .2f [T] \n',B_tooth_no_load)
fprintf('No load flux density in stator yoke = % .2f [T] \n\n',B_yoke_no_load)
fprintf('No load flux density in rotor = % .2f [T] \n\n',B_rotor_no_load)

%Perform sensitivity analysis
%analysis 0: Do nothing
%analysis 1: Change in magnet thickness
%analysis 2: Change in magnet width

analysis = 1;

switch analysis
    
    case 0
    DO NOTHING    
  
    case 1       
    %Change in magnet thickness        
    t_mag_1 = t_mag * 1.1; % 10% increase in magnet thickness
        t_mag_2 = t_mag * 0.9; % 10% decrease in magnet thickness     
        MMF_mag_1 = 2 * H_mag * t_mag_1;
        MMF_mag_2 = 2 * H_mag * t_mag_2;      
        %Magnet characteristic
        figure(5)
        clf
        hold on
        plot(MMF_mag,Phi_mag * 1e3, 'LineWidth', 2)
        plot(MMF_mag_1,Phi_mag * 1e3, 'LineWidth', 2)
        plot(MMF_mag_2,Phi_mag * 1e3, 'LineWidth', 2)
        plot(-MMF_total,Phi_gap * 1e3, 'LineWidth', 2)
        hold off
        xlim([min(MMF_mag_1) 0])
        xlabel('MMF [A-turn]')
        ylabel('Flux [mWb]')
        legend('Magnet thickness 8 mm' , 'Magnet thickness 8.8 mm','Magnet thickness 7.2 mm', 'Load line')
        grid on
        
    case 2
    %Change in magnet width
    w_mag_1 = w_mag * 1.1; % 10% increase in magnet thickness
        w_mag_2 = w_mag * 0.9; % 10% decrease in magnet thickness     
        Phi_mag_1 = B_mag*w_mag_1*L_stack;
        Phi_mag_2 = B_mag*w_mag_2*L_stack;    
        %Magnet characteristic
        figure(5)
        clf
        hold on
        plot(MMF_mag,Phi_mag * 1e3, 'LineWidth', 2)
        plot(MMF_mag,Phi_mag_1 * 1e3, 'LineWidth', 2)
        plot(MMF_mag,Phi_mag_2 * 1e3, 'LineWidth', 2)
        plot(-MMF_total,Phi_gap * 1e3, 'LineWidth', 2)
        hold off
        xlim([min(MMF_mag) 0])
        xlabel('MMF [A-turn]')
        ylabel('Flux [mWb]')
        legend('Magnet width 22 mm' , 'Magnet width 24.2 mm','Magnet width 19.8 mm', 'Load line')
        grid on
    
end
% % 
% % Stator winding design 
% % 
% % Slot area
% % A_slot = ;          % Slot area
% % 
% % Tooth area
% % A_tooth_axial = ;   % Tooth area
% % 
% % Number of turns per coil. Use floor to get integer value.
% % N_turn = ;
% % 
% % Number of strands per turn. Use floor to get integer value.
% % N_strand = ; 
% % 
% % Fill factor
% % f_Cu = ;
% % 
% % fprintf('Tooth area = % .1f [mm^2] \n',A_tooth_axial * 1e6)
% % fprintf('Slot area = % .1f [mm^2] \n',A_slot * 1e6)
% % fprintf('Number of turns per coil = % .0f \n',N_turn)
% % fprintf('Number of strands per turn = % .0f \n',N_strand)
% % fprintf('Fill factor = % .2f \n\n',f_Cu)
% % 
% % Resistance and inductance calculation
% % 
% % Resistance per phase [ohm]
% % Rs = ;  
% % 
% % Permeability of magnet
% % mu_mag = ;
% % 
% % Relative permeability of magnet
% % mu_r_mag = ; 
% % 
% % Magnet reluctance
% % Re_mag = ;
% % 
% % d- and q- axis reluctance
% % Rq = ;
% % Rd = ;
% % 
% % d- and q- axis inductance
% % Nd = ;
% % Nq = ;
% % Ld = ;
% % Lq = ;
% % 
% % fprintf('Resistance per phase = % .1f [mOhm] \n',Rs*1e3)
% % fprintf('Relative permeability of magnet = % .1f \n',mu_r_mag)
% % fprintf('Reluctance of d-flux path = % .1f [H^-1] \n',Rd)
% % fprintf('Reluctance of q-flux path = % .1f [H^-1] \n',Rq)
% % fprintf('d-axis inductance = % .1f [mH] \n',Ld*1e3)
% % fprintf('q-axis inductance = % .1f [mH] \n',Lq*1e3)
% % 
% % Load calculation
% % 
% % Maximum current calculation
% % J_max = ;   % Maximum current density [A/m^2]
% % Isamp = ;   % Maximum current of the machine [A]
% % 
% % d- and q- currents calculation
% % Id = ;      % d-axis current [A]
% % Iq = ;      % q-axis current [A]
% % 
% % d- and q- flux linkage calculation
% % Psid = ;    % d-axis flux linkage [Wb]
% % Psiq = ;    % q-axis flux linkage [Wb]
% % 
% % Electromagnetic torque calculation
% % Te = ;      % Electromagnetic torque [Nm]
% % 
% % Terminal voltage calculation
% % Ud = ;      % d-axis voltage [V]
% % Uq = ;      % q-axis voltage [V]
% % Uamp = ;    % Voltage amplitude [V]
% % 
% % Electrical frequency calculation
% % f = ;       % Electrical frequency [Hz]
% % 
% % Cu losses calculation
% % Pcu = ;     % CU losses [W]
% % 
% % Iron losses calculation
% % Kh = ;      % Hysteresis coefficient 
% % Kc = ;      % Eddy current coefficient
% % Pfe = ;     % Iron losses [W]
% % 
% % Electromagnetic power is same as shaft power by nelecting mechanical losses
% % Pe = ;      % Shaft power [W]
% % 
% % Input power calculation
% % Pi = ;      % Input power [W]
% % 
% % Efficiency calculation 
% % eff = ;     % Efficiency
% % 
% % Power factor calculation
% % S = ;       % Apparent power [W]
% % PF = ;      % Power factor
% % 
% % Output values
% % fprintf('Maximum current % .2f [A] \n', Isamp)
% % fprintf('D-axis current % .2f [A] \n', Id)
% % fprintf('Q-axis current % .2f [A] \n', Iq)
% % fprintf('D-axis flux linkage % .2f [Wb] \n', Psid)
% % fprintf('Q-axis flux linkage % .2f [Wb] \n', Psiq)
% % fprintf('Electromagnetic torque % .2f [Nm] \n', Te)
% % fprintf('Amplitude of terminal voltage % .2f [V] \n', Uamp)
% % fprintf('Electrical frequency % .2f [Hz] \n', f)
% % fprintf('Cu loss % .2f [W] \n', Pcu)
% % fprintf('Iron loss % .2f [W] \n', Pfe)
% % fprintf('Shaft power % .2f [W] \n', Pe)
% % fprintf('Input power % .2f [W] \n', Pi)
% % fprintf('Efficiency % .2f \n', eff )
% % fprintf('Power factor % .2f \n\n', PF)
% % 
% % Current angle sweep
