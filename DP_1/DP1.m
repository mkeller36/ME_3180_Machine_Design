%% Michael Keller Design Problem 1
% part 3
FyLoad = 5525/2;%N
FxLoad = 2011/2;%N
Fxpulley = 1276*cosd(27.7);%N
Fypulley = 1276*sind(27.7);%N
Lengthac = 2.9;%m
Lengthab = 1.950;%m 

%Sum of moments
Fr2y = ((Fypulley*Lengthac)+(FyLoad*Lengthac/2))/Lengthab;
%Sum of y forces 
Fr1y = - Fr2y + Fypulley + FyLoad;
%Sum of x forces 
Fr1x = Fxpulley+FxLoad;

%%

FD = FyLoad;
Fpmm = FD/2900;
F_ay = Fr1y;
F_py = Fypulley
F_by = Fr2y;


x = linspace(0,2905);
V = (F_ay - Fpmm*x).*(x>0)+F_by*(x>1950)-F_py*(x>2900);
plot(x,V)
M = cumtrapz(x,V);
M = M./1000;% mm to m
zero = V * 0;
figure (1)
subplot(211)
plot(x,V)
hold on
plot(x, zero)
xlabel('Distance (mm)')
ylabel('Shear Force (N)')
title('Shear Diagram')

x2 = [x, fliplr(x)];
inBetween = [V, fliplr(zero)];
fill(x2, inBetween, 'cyan');

hold off
subplot(212)
hold on 
plot(x,M)
plot(x,zero)
xlabel('Distance (mm)')
ylabel('Moment (Nmm)')
title('Moment Diagram')


x2 = [x, fliplr(x)];
inBetween = [M, fliplr(zero)];
fill(x2, inBetween, 'cyan');


hold off

%% Calculations
Vmax = max(V);
Mmax = max(M);
Vmin = min(V);
Mmax = abs(min(M));
%height = 3.5 inches
%width = 0.25 inches 
%Max moment at 2900mm
%Max shear at 1966mm
fprintf('Max shear = %f \n',Vmax)
fprintf('Max Moment = %f \n', Mmax)
prompt = 'What is the height in inches? ';
height = 3.5*.0254;
c = height/2;
prompt = 'What is the width in inches? ';
width = .25*.0254;% between 1/8 and 1.5 inches
I = (1/12)*width*height^3;
bendingStress = (Mmax*c)/I;%N/mm^2
force = Fr1x; %pully in x direction
axialStress = force/(height*width);%N/mm^2
totalStress = axialStress + bendingStress;


totalStressI = totalStress*0.000145038;
fprintf('Total Stress = %f PSI \n', totalStressI)
fprintf('Total Stress < 60%% of 50,000 PSI, Sy of steel bar \n')
%% stress concentrations 

x_axis = 5/40;
line = height/40;
%look at chart for Kt
Kt = 1.75;

%% stress in end of beam 

endHeight = .40;
endWidth = width/1000;
endAxial = force/(endHeight*endWidth);
endc = 20;
endI = (1/12)*endWidth*endHeight^3;
endbendingStress = (abs(M(94))*endc)/endI;%N/mm^2
endtotalStress = endAxial + endbendingStress;
endtotalStressI = endtotalStress * .000000000145038;

% multiply max end total stress to get stress at stress concentration 
endtotalStressI = endtotalStressI * Kt;

fprintf('Total Stress in end = %f PSI \n', endtotalStressI)

%% deflection 

w = .9525;%N/mm

%Cantilever - point force - thick 
ymax1 = (F_py*.800^3)/(3*215*10^9*I*2);% Table A-9 in mm
%Cantilever - point force - thin 
ymax1end = (F_py*.150^3)/(3*215*10^9*endI*2);% Table A-9 in mm
%Cantilever - uniform - thick 
ymax2 = (w*.800^4)/(8*215*10^9*I*2);% Table A-9 in mm
%Cantilever - uniform - thin
ymax2end = (w*.150^4)/(8*215*10^9*endI*2);% Table A-9 in mm
totaly2 = (ymax1end+ymax2end)*100/2;
totaly = (ymax1+ymax2)*10/6;%convertion 
fprintf('Total deleflection cantilever thick end = %f m \n', totaly)
fprintf('Total deleflection cantilever thin end = %f m \n', totaly2)
fprintf('These are under .003 m \n')
fprintf('Peak stress is at peak moment, x = 2 m \n')
%dual supported 