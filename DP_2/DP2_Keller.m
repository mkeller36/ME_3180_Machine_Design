%% Michael James Keller
clc
clear

%% Varible setup 

% forces 
drf_N_x = -2013; %drive_roller_force_N_x
belt_force_N_x = -1131; %loaded_belt_side_N_x
belt_force_N_y = -595; %loaded_belt_side_N_y

% Distances 
length_between_supports_mm = 600; %length_between_supports_mm
drive_roller_length_mm = 580; %
drive_roller_diameter_mm = 63.5; %
bearing_width_mm = 12.7; %bearing_width_mm
support_to_mid_pulley_mm = 70; %support_to_mid_pulley_mm
support_to_left_drive_roller_mm = 16.35; %
support_to_right_drive_roller_mm = 16.35; %
right_pin_mm = bearing_width_mm/2 + length_between_supports_mm - 75;
pulley_width_mm = 25; %
pulley_radius_mm = 50; %
key_length_mm = 38; % should not excede 1.5 * shaft diameter =<39; 1.5 inch
key_height_mm = 6.35; % set later - 3/16 - made square with width
key_width_mm = 6.35; % 1/4 inch
keyway_depth_mm = 2.38; % 3/32

% Speeds
conveyor_speed_mm_per_s = 300; %

% Force Application Distances 
rx1_mm = 0; %reaction_1_application_distance_mm
rx2_mm = length_between_supports_mm + bearing_width_mm; %reaction_2_application_distance_mm
drapp_mm = (bearing_width_mm/2) + (length_between_supports_mm/2); %drive_roller_application_distance_mm
pullapp_mm = (bearing_width_mm/2) + length_between_supports_mm + bearing_width_mm + support_to_mid_pulley_mm; %pulley_application_distance_mm

% Prices 
price = 0; 

%% Part 1: Shaft Design Calculations

% Finding reactions 
drive_roller_distributed_load_N_per_mm = drf_N_x /drive_roller_length_mm;

% Sum of Moments 
reaction_2_N_x = -1*((( drf_N_x * (drapp_mm)) + ( belt_force_N_x * (pullapp_mm) ))/rx2_mm);
reaction_2_N_y = -1*(((pullapp_mm) * belt_force_N_y) / rx2_mm);
reaction_2_N_magnitude = sqrt(reaction_2_N_y^2 + reaction_2_N_x^2);
fprintf('Right bearing reaction = %.2f N \n',reaction_2_N_magnitude);

% Sum of forces 
reaction_1_N_x = -1*(drf_N_x + belt_force_N_x + reaction_2_N_x);
reaction_1_N_y = -1*(belt_force_N_y + reaction_2_N_y); 
reaction_1_N_magnitude = sqrt(reaction_1_N_y^2 + reaction_1_N_x^2);
fprintf('Left bearing reaction = %.2f N \n',reaction_1_N_magnitude);

% Finding Torque 
belt_force_N_magnitude = sqrt(belt_force_N_x^2 + belt_force_N_y^2);
belt_torque_N_mm = belt_force_N_magnitude * pulley_radius_mm;
belt_torque_N_m = belt_torque_N_mm/1000;

% Where would we expect the highest stress in the shaft?
shaft_length_mm = bearing_width_mm + length_between_supports_mm + bearing_width_mm + support_to_mid_pulley_mm +pulley_width_mm/2;

% Moment Diagram: X - Z plane 
x = linspace(0,(shaft_length_mm),1000); 
V_xz = (reaction_1_N_x .* (x>0)) + (drive_roller_distributed_load_N_per_mm.*(x-support_to_left_drive_roller_mm).*((x<=(rx2_mm - support_to_right_drive_roller_mm))&(x>=rx1_mm + support_to_left_drive_roller_mm))) + drf_N_x.*(x>rx2_mm - support_to_right_drive_roller_mm) + ((reaction_2_N_x).*(x>(rx2_mm))) + (belt_force_N_x.*(x>pullapp_mm));
M_xz = cumtrapz(x,V_xz);
figure (1)
subplot(211)
xlabel('Distance (mm)')
ylabel('Shear (N)')
plot(x,V_xz)
title('X - Z plane Shear')
grid on
subplot(212)
xlabel('Distance (mm)')
ylabel('Moment (N*mm)')
plot(x,M_xz)
title('X - Z plane Moment')
grid on

% Moment Diagram: Y - Z plane 
V_yz = (reaction_1_N_y .* (x>0)) + ((reaction_2_N_y).*(x>(rx2_mm))) + (belt_force_N_y.*(x>pullapp_mm));
M_yz = cumtrapz(x,V_yz);
figure (2)
subplot(211)
xlabel('Distance (mm)')
ylabel('Shear (N)')
plot(x,V_yz)
title('Y - Z plane Shear')
grid on
subplot(212)
xlabel('Distance (mm)')
ylabel('Moment (N*mm)')
plot(x,M_yz)
title('Y - Z plane Moment')
grid on

% sum of squares moments 
V_prime = sqrt(V_yz.^2 + V_xz.^2);
M_prime = sqrt(M_yz.^2 + M_xz.^2);
figure (3)
subplot(211)
xlabel('Distance (mm)')
ylabel('Shear (N)')
plot(x,V_prime)
title('resultant plane Shear')
grid on
subplot(212)
xlabel('Distance (mm)')
ylabel('Moment (N*mm)')
plot(x,M_prime)
title('resultant plane Moment')
grid on

% torque graph 
torque = belt_torque_N_mm .*((x> (bearing_width_mm/2)+ 600 - 75)&(x<shaft_length_mm));
figure (4)
plot(x,torque)
xlabel('Distance (mm)')
ylabel('Torque (N*mm)')
title('Torque')
grid on

% I am expecting to find the highest stress in that shaft at the right most
% retaining ring as the stress concentrations at this point will be great. 

%% Part 1: Finding Diameter 

close all

% Finding k_a: surface finish factor (table 6-2), McMaster says ground          1346K39	
s_ut_MPA = 2206; %http://www.carbidedepot.com/formulas-hardness.htm rockwell c25 (855)  use c60 (2206) instead 
k_a = 1.58*s_ut_MPA^(-0.085); % should be .8901

% stress concentrations 
% shaft_length_mm

% Finding Kt (bending) and Kts (tortional): End mill keyseat table 7-1
k_t = ones(1,1000);
k_t(x>((12.7/2)+600+12.7)&x<(12.7/2)+600+12.7+4) = 4.8; % was 5, changed to 4.8 retaining ring 
k_t(x>((12.7/2)+600-75)&x<(12.7/2)+600+1-75) = 2; % 2. for pins 
k_t(x>((12.7/2)+75)&x<(12.7/2)+1+75) = 2; % 2. for pins 
k_t(x>(shaft_length_mm-key_length_mm)) = 2.14; % 45mm key length
k_t(x>(shaft_length_mm-support_to_left_drive_roller_mm)) = 00; 

% going back to find real stress concentration using table A-16-18
% retaining ring 97633A300	
% used
% https://www.fastenermart.com/files/external-retaining-ring-specifications.pdf
% for values, edge margin is 0.09 inches 
ring_t = (1 - 0.94)/2;
ring_r = 0.005;
ring_a = 0.046;
ring_r_per_t = ring_r/ring_t;
ring_a_per_t = ring_a/ring_t;
% K_t is 4.8;
% k_ts is 2.75;

%k_ts = 3.0;
k_ts = ones(1,1000);
k_ts(x>((12.7/2)+600+12.7)&x<(12.7/2)+600+12.7+4) = 2.75; % initial guess 5, changed to 2.75
% input 2 
k_ts(x>((12.7/2)+600-75)&x<(12.7/2)+600+1-75) = 3.; % 3.6 for pins 
k_ts(x>((12.7/2)+75)&x<(12.7/2)+1+75) = 3.; % 3.6 for pins
k_ts(x>(shaft_length_mm-key_length_mm)) = 3.0; % 45mm key length
k_ts(x>(shaft_length_mm-support_to_left_drive_roller_mm)) = 00; 

% Finding k_b: size factor page 296 (6-20)
d = 31.75; % experimental_diamter_mm
if (d < 51 )
    k_b = 1.24*d^(-0.107);
elseif (d > 51 && d < 254)
    k_b = 1.51*d^(-0.157);
end
 % k_b should be .8899 - maybe 

% Finding Se value 
k_c = 1; % Loading Type, combined 
k_d = 1; % temp factor, normal temp 
k_e = .62; % reliability factor 99.9999, It will pretty much never fail and that is how I want it
% table 6-5, page 301, argue certain reliability factor

% Going back to find stress concentrations for pins, table A-15_10 
% pins used 98381A552
% .25 inch pin 
pin_diameter = .25;
pin_d_per_D = pin_diameter / 1.25;
% x value 0.2
% table A-15-10 & A-15-11
% gives moment stress concentration of 2.2
% gives tortion stress concentration of 3.2 

% finding k_f - page 303 - Have no idea, k_f = k_t,
% page 304, notch sensitivity Neuber equation to find k_f
% !!!!!!!!!!!!!!!!!!!!! use page 304 to find q based off r and a 
% q_r is 3.17 and q is 1.787
q = 0.85; % from table 6-20 with Sut 0.855 and .5mm radius, was 0.8 
q_shear = .85; % was 0.8
% k_f = 1 + q*(k_t - 1); %estimate k_f for first itteration 
 k_fs = 1 + q_shear*(k_ts-1);
% used k_f and k_fs in goodman equation 
k_f = 1;

% finding s_y, material prop
s_y = 517; % 517 for 75,000 | 310 for 45,000

% Finding s_e
s_e_prime = s_ut_MPA / 2;

% finding s_e_prime - Merin factors 
s_e = k_a * k_b * k_c * k_d * k_e * k_f * s_e_prime; 

% find safety factor n in shaft
n = 1.5;

% find mean torque 
t_m = torque; %N*mm 

% find bending moment amplitude 
m_a = M_prime; 

% DE - Goodman Equation 
d = ((16*n/pi).*(((2.*k_t.*m_a)./s_e)+(sqrt(3).*k_ts.*t_m)./s_ut_MPA)).^(1/3);

figure (5)
plot(x,d)
xlabel('Distance (mm)')
ylabel('Minimum required diameter')
title('Diameter required vs distance')
grid on

% print varibles - location of greatest necessary diameter should be at
% location of greatest stress
[val, loc] = max(d);
fprintf('Required Diameter = %.2f mm \n',val);

fprintf('Required Diameter = %.2f inches \n',(val*0.0393701)); 
max_moment = M_prime(loc);
max_torque = torque(loc);

fprintf('Required Length = %.2f inches \n',(shaft_length_mm*0.0393701)); 

% maybe put equations into van mises and tresca stresses
fprintf('Max moment = %.2f Nmm \n',max_moment)
fprintf('Max torque = %.2f Nmm \n',max_torque)

%   Move shaft Diamter up one standard size to 1.25 inches, because even
%   though one inch may just barely work, the stress concentrations would  
%   increase and push it over one inch, potentially breaking the shaft,
%   	5947K55

%   retaining ring part number      97633A320	
% add 0.09 inches to shaft to have between retaining ring and end of shaft 
%% Verifying with correct shaft diameter yeild check

diameter = 32.75;% was 25.5
c = diameter/2;
i = (pi*diameter^4)/64;
% sigma_prime_m = sqrt(((k_f.*max_moment.*c)./i).^2+(3*max_torque).^2); % set equal to torque 
% sigma_prime_a = sqrt(((k_fs.*t_m.*c)./(i/2)).^2+(3.*max_torque).^2); % set equal to bending moment 
sigma_prime_m = max_torque;
sigma_prime_a = max_moment;

simple_shaft_equation = n*max((((k_fs.*torque.*c)./(i/2)))./s_y + (((k_f.*M_prime.*c)./i))./s_y);
% if this is less then 1, the yeil check succeeds 
if (simple_shaft_equation < 1)
    fprintf('Shaft passes yeild check \n')
else
    fprintf('Shaft fails yeild check \n')
end
% Do von mises check for shaft with given moments and torque. 

%% Part 2: Pulley Key design 

%   Key width and height from table 7-6, page 383
%   Diameter of shaft being 1", the key dimentions should be as put in
%   above, or use stock for cheaper option, 98510A125
%   Key selected    Part#: 	98870A215	    Sy = 1606 MPA	 
key_shear_N_per_mm2 = 4*belt_torque_N_mm/(diameter * key_height_mm * key_length_mm);
key_torque_N_per_mm2 = 2*belt_torque_N_mm/(diameter * key_width_mm * key_length_mm);
% find Von Mises of stress to make sure it doesnt break
key_Von_mises = 1.2*sqrt(key_shear_N_per_mm2^2+(3*key_torque_N_per_mm2^2)); %add safey factor 
fprintf('Key Stress = %.2f MPA \n',key_Von_mises)
% This is well below the yeild stress 

%% Part 2: Set screw design 

% Page 381 - table 7-4
% No axial load so what set screw we chose doesnt matter that much 

%   Screw selected      Part #: 94105A142

%% Part 2: Connection Pin design 

% One pin must resist all torque - 
force_shaft_Nmm = 1.2 * belt_torque_N_mm /(diameter/2); %add saftey factor 
pin_area = 2*force_shaft_Nmm/V_prime(751);
pin_radius = sqrt(pin_area/pi);
fprintf('Minimum pin radius = %.2f mm \n',pin_radius)
% rounded up to .25 inch pin because it was cheaper than a 3/16
% 98381A552
%% Part 3: Roller Bearing design 

% part number:  5709K85
% table 7-9 for fit, page 389, H7, g6
% The shaft directly from McMaster is toleranced to LC 1 - h5
% The bearing directly from McMaster comes toleranced to LC 1 - H6
% This is a tight locational clearance 
% from machinery handbook, page 32
% + 1 thou
% - .6 thou 
% table A-12 for shaft fits and limits 
% need to make shaft H,7 tolerance, 
% bearing has to be thinner than 12.7 mm
f_dynamic_n = 9341.265; %newtons 
f_static_n = 44482.22; % newtons 
a = 10/3; %10/3 for thrust bearings 
roller_circumference = 199.49;
l_catalog = 1*10^6;
l_design = conveyor_speed_mm_per_s / roller_circumference * 3.154*10^7; % seconds per year 
f_catalog = reaction_2_N_magnitude * (l_design/l_catalog)^(1/a);
if((f_catalog < f_dynamic_n) && (reaction_2_N_magnitude < f_static_n))
    fprintf('Bearing passes \n')
else
    fprintf('Bearing fails \n')
end

% The outside of the bearing should have an interference fit with the
% journal of the rails. The tolerancing for the out side of the bearing is 
% +0.001 / -0.0, similar to a n5 fit. We want an interference fit, H7, p6
% on Machine design handbook table 12 ANSI B4.1-1967 (R1999), giving us 0.0016 inches of
% interference worse case. Because of this, the journal the bearing sits
% inside should be should have a nominal measurement of 2.328125 inches
% from the drawing of the bearing with tolerances +0.001 / -0.0016 inches. 
% Because both the shaft and the bearing come within these tolerances, we
% will not need to machine them, saving on cost. 