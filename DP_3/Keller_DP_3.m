%% Michael James Keller
clc
clear
% ¯\_(?)_/¯
%% Varible setup 

% shafts
diameter_motor_shaft = 0.75; %mm
diameter_intermediate_shaft = 0.75; %19.05 mm to 0.75 inches  
diameter_drive_shaft = 1.25; %31.75 mm to 1.25 inches


% Other
step = 2;

%% stepping through 14-19 design layout, page 767

% function
    % Load 
    torque_out = 12*47.13022; %69.3Nm, convert to inch-lbs
    % Speed - RPM 
    rpm_in = 360; %360 rpm to omega
    rpm_out = 90; % 90 rpm to omega
    % Reliability 
    safety_factor = 1.5;
    % life 
    life_hours = 24*365*5;
    life_minutes = life_hours*60;
    % K_o
    
% Unquantifiable risk 

% Tooth System 
    % Pressure Angle 
        phi = 14.5; % pressure angle
    % Addendum 
        % a pitch of 4 seems to give standardized numbers for diameter, 8.5
        % and 4.25 
        %pitch = [2 2.25 2.5 3 4 6 8 10 12 16 20 24 32 40 48]; % Set pitch:  p = N/d, 16 is availible at McMaster
        pitch = 10;
        addendum = 1./pitch;
    % Dedendum 
        deddendum = 1.25./pitch;
    % Root fillet radius 

% Gear Ratio
    % Gear Ratio
        % Page 692, example 13-3
        gear_ratio = rpm_out/rpm_in;
        ratio_per_step = gear_ratio^(1/step);
        fprintf("Question 1: \n")
        fprintf("Reduction per step = %.2f or 2:1\n",ratio_per_step)
    % N_p_theoretical (number of teeth on pinion)
        % Page 678, Equation 13-11
        m = 1/ratio_per_step; % gear teeth ration, per step 
        k = .8; % 1 for full depth teeth, 0.8 for stub teeth
        % need to use stub teeth
        N_p_theory = ceil(((2*k)/((1+2*m)*sind(phi)^2))*(m+sqrt(m^2+(1+2*m)*sin(phi)^2)));
    % N_g_theoretical (number of teeth on gear)
        % Page 679, Equation 13-12
        N_g_max = (N_p_theory^2 * sind(phi)^2 - 4*k^2)/(4*k - 2*N_p_theory*sind(phi)^2); 
        % Find required mating gear for pinion 
        N_g_theory = N_p_theory * (1/ratio_per_step);
        
% Choseing gears 
    %6867K79 for pinion 1 and 3, 
    %6325K9 for gear 2 and 4 
    gear1.name = '6867K79';
    gear1.pitch = 10;
    gear1.pitch_diameter = 3;
    gear1.face_width = 1;
    gear1.teeth = 30;

    gear2.name = '6325K9';
    gear2.pitch = 10;
    gear2.pitch_diameter = 6;
    gear2.face_width = 1;
    gear2.teeth = 60;
        fprintf("Number of teeth on pinion 1 = %2.0f \n",gear1.teeth)
        fprintf("Number of teeth on pinion 3 = %2.0f \n",gear1.teeth)
        fprintf("Number of teeth on gear 2 = %2.0f \n",gear2.teeth)
        fprintf("Number of teeth on gear 4 = %2.0f \n",gear2.teeth)
% Quality Number
% Diametral Pitch 
    diameter_pinion = gear1.pitch_diameter; %N_p_theory ./ pitch;
    diameter_gear = gear2.pitch_diameter; %N_g_theory ./ pitch;
    % gears 1 and 3 are pinions, 2 and 4 are gears 
    % 
    diameter_1 = diameter_pinion;
    diameter_2 = diameter_gear;
    diameter_3 = diameter_pinion;
    diameter_4 = diameter_gear; 
% Face Width 
    face_width = gear1.face_width; % from Mcmaster 
% Pinion Material 
    s_t_steel = 25000; %psi
    s_t_brass = 10000; %psi, page 742, hardness 350
    s_c_steel = 141000; %psi, change if you harden it 
    s_c_brass = 40000; %psi
% Gear Material 

%% Solving Example 
% Pinion Bending 
    F = 1;% 4*pi/pitch; % face width theoretical 
    
%% Finding omega for each gear 
rpm_2 = rpm_in*(N_p_theory/N_g_theory);
rpm_3 = rpm_2; % on same shaft 
fprintf("Question 4: \n")
fprintf("Intermediate shaft rotational speed = %3.0f rpm \n",rpm_2)
%% Finding Torque on each gear 
torque_3 = torque_out*(rpm_out/rpm_3); %pound feet
torque_2 = torque_3*(rpm_3/rpm_2);
torque_in = torque_2*(rpm_2/rpm_in);
fprintf("Torque on intermediate shaft = %f inch-Lbs \n",torque_2)

%% w_t, transmitted tangental load, on each gear - page 697, 13-14 
w_t_1 = 2*torque_in/diameter_1; %pounds 
w_t_2 = 2*torque_2/diameter_2;
w_t_3 = 2*torque_3/diameter_3;
w_t_4 = 2*torque_out/diameter_4; 
fprintf("Question 5: \n")
fprintf("Transmitted tangental load on drive shaft = %f Lbs \n",w_t_4)
%% w_r, transmitted radial load, on each gear - page 701, 13-15
w_r_1 = w_t_1*tand(phi); %pounds 
w_r_2 = w_t_2*tand(phi);
w_r_3 = w_t_3*tand(phi);
w_r_4 = w_t_4*tand(phi);
% Use this for bending in shafts 
fprintf("Transmitted radial load on drive shaft = %f Lbs \n",w_r_4)
% if the actual radial load is less than what we used to calculate shaft
% diameter, then we should be good to use a previously designed shaft
% diameter 
fprintf("The forces are just barely greater than what was used for the \ncalculations to design the original shaft but the calculations should still be reran.\nPlugging in the numbers a shaft of 1.25 inches should still work.\n")
%% Lewis form factors  - 730
Y_1 = 0.318; % N = 30, used chart 14-2 on page 730
Y_2 = 0.355; % N = 60, 
Y_3 = 0.318; % N = 30,
Y_4 = 0.355; % N = 60,

%% K_v for machined gears - 731
v_1 = (diameter_1/12)*pi*rpm_in;
K_v_1 = (1200+v_1)/1200; 

v_2 = (diameter_2/12)*pi*rpm_2;
K_v_2 = (1200+v_2)/1200; 

v_3 = (diameter_3/12)*pi*rpm_3;
K_v_3 = (1200+v_3)/1200;

v_4 = (diameter_4/12)*pi*rpm_out;
K_v_4 = (1200+v_4)/1200; 

%% Bending stress in tooth
% from figure 14-6, page 745
J_1 = 0.385;
J_2 = 0.425;
J_3 = 0.385;
J_4 = 0.425;
% from page 752 
c_mc = 1; %uncrowned 
c_pf_1 = 0.025;
c_pf_2 = 0.025;
c_pm = 1.1;
c_ma = 0.127 + 0.0158*F;
c_e = 1;
k_m_1 = 1+c_mc*(c_pf_1*c_pm + c_ma*c_e);
k_m_2 = 1+c_mc*(c_pf_2*c_pm + c_ma*c_e);
k_b = 1;
k_o = 1.25;
k_s = 1;
sigma_1 = K_v_1*w_t_1*k_s*k_o*((pitch/F)*k_m_1*k_b/J_1);
sigma_2 = K_v_2*w_t_2*k_s*k_o*((pitch/F)*k_m_2*k_b/J_2);
sigma_3 = K_v_3*w_t_3*k_s*k_o*(pitch*k_m_1*k_b/(J_3*F));
sigma_4 = K_v_4*w_t_4*k_s*k_o*(pitch*k_m_2*k_b/(J_4*F));
% compare these to s_t, yeild in tention 

%% Calculating safety factor, must be more than 1.5 
load_cycles = rpm_out*life_minutes;
y_n = .9;
k_t = 1;
k_r = 1.25;
safety_factor_theoretical = (s_t_steel * y_n/(k_t*k_r))/sigma_4 ;
% if this is greater than 1.5 we are good. 
if (safety_factor_theoretical > 1.5)
    fprintf("Passes bending \n")
else
    fprintf("Fails bending \n")
end

%% Finding minimum gear face width

sigma_theoretical = (s_t_steel*y_n/(k_t*k_r))/safety_factor;
f_min_pinion = w_t_3*k_o*K_v_4*k_s*(pitch/sigma_theoretical)*((k_m_1*k_b)/J_3);
f_min_gear = w_t_4*k_o*K_v_4*k_s*(pitch/sigma_theoretical)*((k_m_2*k_b)/J_4);
fprintf("Question 6: \n")
fprintf("Bending Stress: %f \n",sigma_4)
fprintf("Material Properties: \n")
fprintf("sigma_theoretical: %f \n",sigma_theoretical)
fprintf("S_c steel: %f \n",s_t_steel)
fprintf("Smallest face for gear 3 = %f \n",f_min_pinion)
fprintf("Smallest face for gear 4 = %f \n",f_min_gear)
%% Finding Sigma c 
c_f = 1;
I = (cosd(phi)*sind(phi)/2)*(ratio_per_step)/(ratio_per_step+1);
%cp find on page 736 and 749
c_p = ((pi*(((1-.3^2)/(28*10^6))+((1-.3^2)/(28*10^6)))))^(-1/2);
sigma_c_1 = c_p*(w_t_1*k_o*K_v_1*k_s*k_m_1*c_f/(diameter_1*F*I))^(1/2);
sigma_c_2 = c_p*(w_t_1*k_o*K_v_1*k_s*k_m_1*c_f/(diameter_2*F*I))^(1/2);
sigma_c_3 = c_p*(w_t_3*k_o*K_v_3*k_s*k_m_1*c_f/(diameter_3*F*I))^(1/2);
sigma_c_4 = c_p*(w_t_4*k_o*K_v_4*k_s*k_m_2*c_f/(diameter_4*F*I))^(1/2);

fprintf("Contact Stress: %f \n",sigma_c_4)
%% Finding wear factor of safety
z_n = 0.9; % page 755
c_h = 1; %page 753, same material, same hardness
s_h = (s_c_steel*z_n*c_h/(k_t*k_r))/(sigma_c_3); % use 3 because it is biger 
% if greater than 1.5 we are good 
if (s_h > 1.5)
    fprintf("Passes wear \n")
else
    fprintf("Fails wear \n")
end

%% Finding minimum gear face width wear 
sigma_c_theoretical = s_c_steel*z_n*c_h/(k_t*k_r)/sqrt(safety_factor);
f_pinion_theoretical = c_p^2*(w_t_3*k_o*K_v_3*k_s*k_m_1*c_f/(diameter_3*sigma_c_theoretical^2*I));
f_gear_theoretical = c_p^2*(w_t_4*k_o*K_v_4*k_s*k_m_2*c_f/(diameter_4*sigma_c_theoretical^2*I));
fprintf("Smallest face for gear 3 = %f \n",f_pinion_theoretical)
fprintf("Smallest face for gear 4 = %f \n",f_gear_theoretical)
%% Estimating Horse power rating page 732, 14-2 

%% Gears
% Helpful example at 760, 14-4 and till 14-8
% 14-19 goes over how to design a gear set 

fprintf("Question 7: \n")
fprintf("Pinion number %s \n",gear1.name)
fprintf("Gear number %s \n",gear2.name)
fprintf("Question 8: \n")
fprintf("The gears would need to be hardened to a brinell hardness of about 400 to work for my calculations.\n")