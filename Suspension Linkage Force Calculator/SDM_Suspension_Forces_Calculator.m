%% SDM Suspension Linkage Force Calculator (Procedural version)
% Version: 2.6.1

%% About
% Author: Hayden Redmond
% Contact: hayden@redmondengineering.com

%% NOTICE
% For educational use only. Do not distribute. Usage permitted only with written consent.

%% User Inputs
clc; clear; close all

% Number of loadcases (must update switch statement in AppliedForces too)
Num_Loadcases = 10;

% Vehicle/Environment Parameters:
Load_Distribution = 0.485;
Downforce_Distribution = 0.49;
Coeff_Drag = 1.223;
Coeff_Downforce = 0.069;
Frontal_Area = 1643; % inches^2
Air_Density = 4.427*10^-5; 
Total_Mass = 620; % lbs; Mass of vehicle + driver
CG_Height = 11.2; % inches
Front_Trackwidth = 47.5; % inches
Rear_Trackwidth = 47; % inches
Wheelbase = 60.125; % inches
Mz_Max = 60; % ft lb
Mx_Max = 20; % ft lb            
Factor_of_Safety = 1.5;      

%% Main executable

% Import linkage dimensions and material properties
Linkage_Dimensions = table2array(readtable('Linkage_Properties.xlsx','Sheet',1,'Range','B3:G4'));
Material_Properties = table2array(readtable('Linkage_Properties.xlsx','Sheet',2,'Range','B2:C2'));

for big_loop = 1:4

% Import suspension points and organize
if big_loop == 1      % (FR)
    All_Points_Table = readtable('Suspension_Points.xlsx','Sheet',1);
    Inboard_Points = table2array(All_Points_Table(1:6,2:4));
    Outboard_Points = table2array(All_Points_Table(8:13,2:4));

    Moment_Arm = [0 0 0];
    Trackwidth = Front_Trackwidth;
    Half_Trackwidth = [0 Front_Trackwidth/2 0];

elseif big_loop == 2 % (FL)
    All_Points_Table = readtable('Suspension_Points.xlsx','Sheet',1);
    Inboard_Points = table2array(All_Points_Table(1:6,2:4));
    Outboard_Points = table2array(All_Points_Table(8:13,2:4));
    % FLIPPING SIGN OF Y COORDS
    Inboard_Points = [Inboard_Points(:,1) -Inboard_Points(:,2) Inboard_Points(:,3)];
    Outboard_Points = [Outboard_Points(:,1) -Outboard_Points(:,2) Outboard_Points(:,3)];

    Moment_Arm = [0 0 0];
    Trackwidth = Front_Trackwidth;
    Half_Trackwidth = [0 -Front_Trackwidth/2 0];

   
elseif big_loop == 3 % (RR)
    All_Points_Table = readtable('Suspension_Points.xlsx','Sheet',2);
    Inboard_Points = table2array(All_Points_Table(1:6,2:4));
    Outboard_Points = table2array(All_Points_Table(8:13,2:4));

    Moment_Arm = [-60.125 0 0];
    Trackwidth = Rear_Trackwidth;
    Half_Trackwidth = [0 Rear_Trackwidth/2 0];


elseif big_loop == 4 % (RL)
    All_Points_Table = readtable('Suspension_Points.xlsx','Sheet',2); 
    Inboard_Points = table2array(All_Points_Table(1:6,2:4));
    Outboard_Points = table2array(All_Points_Table(8:13,2:4));
    % FLIPPING SIGN OF Y COORDS
    Inboard_Points = [Inboard_Points(:,1) -Inboard_Points(:,2) Inboard_Points(:,3)];
    Outboard_Points = [Outboard_Points(:,1) -Outboard_Points(:,2) Outboard_Points(:,3)];

    Moment_Arm = [-60.125 0 0];
    Trackwidth = Rear_Trackwidth;
    Half_Trackwidth = [0 -Rear_Trackwidth/2 0];

end

% Call basic functions
UnitVectors = UnitVectorize(Inboard_Points, Outboard_Points);
Moment_Arm_Positions = MomentArmLength(Moment_Arm, Inboard_Points);
A_Cross_Products = A_CrossProduct(Moment_Arm_Positions, UnitVectors);

% Assemble unit vectors and moments into vector, A
A_UnitVectors = UnitVectors';
A_Moments = A_Cross_Products';
A = cat(1,A_UnitVectors,A_Moments);
  
% Call LoadcasesForces for each loadcase, store all outputs in one array
for i = 1:Num_Loadcases
    [Forces, G_Parameters] = AppliedForces(i, Load_Distribution, Downforce_Distribution, ...
    Coeff_Drag, Coeff_Downforce, Frontal_Area, Air_Density, Total_Mass, CG_Height, Trackwidth, ...
    Wheelbase);
    
    Temp_Forces_and_G_Parameters = cat(2,Forces, G_Parameters);
    Forces_and_G_Parameters(i,:) = Temp_Forces_and_G_Parameters;
    All_Forces(i,:) = Forces;
end

% Transpose All_Forces for later use in vector, b (applied forces and
% moments)
b_All_Forces = All_Forces';

% Call CrossProduct to find M from fx/fy/fz; seen in C54:E56 in excel calc
for i = 1:Num_Loadcases
    % M from Fx
    Fx_Forces = [All_Forces(i,1) 0 0];
    M_From_Fx(i,:) = b_CrossProduct(Fx_Forces, Half_Trackwidth);
    % M from Fy
    Fy_Forces = [0 All_Forces(i,2) 0];
    M_From_Fy(i,:) = b_CrossProduct(Fy_Forces, Half_Trackwidth);
    % M from Fz    
    Fz_Forces = [0 0 All_Forces(i,3)];
    M_From_Fz(i,:) = b_CrossProduct(Fz_Forces, Half_Trackwidth);     
end

% Organize all moments for vector, b (applied forces and moments)
for j = 1:Num_Loadcases
    b_All_Moments(j,:) = [M_From_Fz(j,1)+Mx_Max M_From_Fy(j,2) (M_From_Fx(j,3))+Mz_Max];
end

% Setup vector, b
b = cat(1, -b_All_Forces, b_All_Moments');

% Calculate vector, X
for k = 1:Num_Loadcases
    X(:,k) = A\b(:,k);
end

% Output vector, X to Forces_In_Linkages
if big_loop == 1      % (FR)
    writematrix(X,'Forces_In_Linkages.xlsx','Sheet',1,'Range','C3:L8')

elseif big_loop == 2 % (FL)
    writematrix(X,'Forces_In_Linkages.xlsx','Sheet',2,'Range','C3:L8')
   
elseif big_loop == 3 % (RR)
    writematrix(X,'Forces_In_Linkages.xlsx','Sheet',3,'Range','C3:L8')

elseif big_loop == 4 % (RL)
    writematrix(X,'Forces_In_Linkages.xlsx','Sheet',4,'Range','C3:L8')
end

end

%% Checking for yielding in linkages

% Import recently amended force calculations
FR_Linkage_Forces = table2array(readtable('Forces_In_Linkages.xlsx','Sheet',1,'Range','C3:L8'));
FL_Linkage_Forces = table2array(readtable('Forces_In_Linkages.xlsx','Sheet',2,'Range','C3:L8'));
RR_Linkage_Forces = table2array(readtable('Forces_In_Linkages.xlsx','Sheet',3,'Range','C3:L8'));
RL_Linkage_Forces = table2array(readtable('Forces_In_Linkages.xlsx','Sheet',4,'Range','C3:L8'));

% Calculate axial stress and factor of safety in all linkages for all loadcases
Row = 1; % Used for indexing within Stresses_and_Safety
for CORNER = 1:4 % loop through each corner
    for LINKAGE = 1:6 % loop through each linkage
        for LOADCASE = 1:10 % loop through each loadcase
            % Calculate area of specified linkage 
            Linkage_Area = pi*(Linkage_Dimensions(1,LINKAGE)/2)^2 - pi*((Linkage_Dimensions(1,LINKAGE)-Linkage_Dimensions(2,LINKAGE))/2)^2;   % changes for each LINKAGE
        
            switch CORNER
                case 1 % FR
                    Axial_Stress = abs(FR_Linkage_Forces(LINKAGE,LOADCASE)) / Linkage_Area; % Calculate axial stress
                    Stresses_and_Safety(Row,1) = CORNER; % Denote corner tested
                    Stresses_and_Safety(Row,2) = LINKAGE; % Denote linkage tested
                    Stresses_and_Safety(Row,3) = LOADCASE; % Denote loadcase tested
                    Stresses_and_Safety(Row,4) = Axial_Stress; % Save to Stresses_and_Safety_FR
                    FOS = Material_Properties(1) / Axial_Stress; % Factor of safety based on yielding
                    Stresses_and_Safety(Row,5) = FOS;

                    if FOS < Factor_of_Safety % Checking if calculated factor of safety is less than desired factor of safety
                        FOS_Check = 0;
                    else
                        FOS_Check = 1;
                    end
                    Stresses_and_Safety(Row,6) = FOS_Check;

                case 2 % FL
                    Axial_Stress = abs(FL_Linkage_Forces(LINKAGE,LOADCASE)) / Linkage_Area; % Calculate axial stress
                    Stresses_and_Safety(Row,1) = CORNER; % Denote corner tested
                    Stresses_and_Safety(Row,2) = LINKAGE; % Denote linkage tested
                    Stresses_and_Safety(Row,3) = LOADCASE; % Denote loadcase tested
                    Stresses_and_Safety(Row,4) = Axial_Stress; % Save to Stresses_and_Safety_FR
                    FOS = Material_Properties(1) / Axial_Stress; % Factor of safety based on yielding
                    Stresses_and_Safety(Row,5) = FOS;

                    if FOS < Factor_of_Safety % Checking if calculated factor of safety is less than desired factor of safety
                        FOS_Check = 0;
                    else
                        FOS_Check = 1;
                    end
                    Stresses_and_Safety(Row,6) = FOS_Check;

                case 3 % RR
                    Axial_Stress = abs(RR_Linkage_Forces(LINKAGE,LOADCASE)) / Linkage_Area; % Calculate axial stress
                    Stresses_and_Safety(Row,1) = CORNER; % Denote corner tested
                    Stresses_and_Safety(Row,2) = LINKAGE; % Denote linkage tested
                    Stresses_and_Safety(Row,3) = LOADCASE; % Denote loadcase tested
                    Stresses_and_Safety(Row,4) = Axial_Stress; % Save to Stresses_and_Safety_FR
                    FOS = Material_Properties(1) / Axial_Stress; % Factor of safety based on yielding
                    Stresses_and_Safety(Row,5) = FOS;

                    if FOS < Factor_of_Safety % Checking if calculated factor of safety is less than desired factor of safety
                        FOS_Check = 0;
                    else
                        FOS_Check = 1;
                    end
                    Stresses_and_Safety(Row,6) = FOS_Check;

                case 4 % RL
                    Axial_Stress = abs(RL_Linkage_Forces(LINKAGE,LOADCASE)) / Linkage_Area; % Calculate axial stress
                    Stresses_and_Safety(Row,1) = CORNER; % Denote corner tested
                    Stresses_and_Safety(Row,2) = LINKAGE; % Denote linkage tested
                    Stresses_and_Safety(Row,3) = LOADCASE; % Denote loadcase tested
                    Stresses_and_Safety(Row,4) = Axial_Stress; % Save to Stresses_and_Safety_FR
                    FOS = Material_Properties(1) / Axial_Stress; % Factor of safety based on yielding
                    Stresses_and_Safety(Row,5) = FOS;

                    if FOS < Factor_of_Safety % Checking if calculated factor of safety is less than desired factor of safety
                        FOS_Check = 0;
                    else
                        FOS_Check = 1;
                    end
                    Stresses_and_Safety(Row,6) = FOS_Check;
                     
            end
        Row = Row + 1; 
        end
    end
end

% Check Stresses_and_Safety for FOS_Check == 0, save to a table in which
% linkages in which corner for which loadcase that safety factor is not met
All_Yielding_Linkages = {};
for SAS_Index = 1:height(Stresses_and_Safety)
    if Stresses_and_Safety(SAS_Index,6) == 0
        switch Stresses_and_Safety(SAS_Index,1)
            case 1
                Corner_Identifier = 'Front Right';
            case 2
                Corner_Identifier = 'Front Left';
            case 3
                Corner_Identifier = 'Rear Right';
            case 4
                Corner_Identifier = 'Rear Left';
        end

        switch Stresses_and_Safety(SAS_Index,2)
            case 1
                Linkage_Identifier = 'Upper-Fore';
            case 2
                Linkage_Identifier = 'Upper-Aft';
            case 3
                Linkage_Identifier = 'Lower-Fore';
            case 4
                Linkage_Identifier = 'Lower-Aft';
            case 5
                Linkage_Identifier = 'Push/Pull';
            case 6
                Linkage_Identifier = 'Tie/Toe';
        end
             
    % Add linkages and their info that yield to a table
        % Corner, linkage type, loadcase, FOS
    Yielding_Linkage = {Corner_Identifier Linkage_Identifier Stresses_and_Safety(SAS_Index,3) Stresses_and_Safety(SAS_Index,5)};
    All_Yielding_Linkages = [All_Yielding_Linkages; Yielding_Linkage];
    end
end

if isempty(All_Yielding_Linkages) == 1
    All_Yielding_Linkages_Table = [];
else
    All_Yielding_Linkages_Table = cell2table(All_Yielding_Linkages, 'VariableNames', {'Corner', 'Linkage', 'Loadcase', 'Safety Factor'});
end

% Inform user if any linkages do not meet safety factor
if isempty(All_Yielding_Linkages_Table)
    fprintf('------------> Safety factor (yielding) of %.2f predicted in all linkages! <-------------\n', Factor_of_Safety)
else
    fprintf('------------> Safety factor of %.2f (yilelding) **not** predicted in all linkages! <-------------\n', Factor_of_Safety)
    disp(All_Yielding_Linkages_Table)
    % Create flashing red warning figure
    h = figure('Name','WARNING','NumberTitle','off',...
               'Color','r','MenuBar','none','ToolBar','none');
    ax = axes('Parent',h,'Position',[0 0 1 1]);
    axis off
    text(0.5,0.5,'⚠WARNING: Yielding Linkages Found!⚠', ...
         'Color','w','FontSize',18,'FontWeight','bold','HorizontalAlignment','center');

    % Flash red and black
    for t = 1:12
        if mod(t,2)==0
            set(h,'Color','r');
        else
            set(h,'Color','k');
        end
        pause(0.3);
    end
end


%% Functions

function Vectorized = UnitVectorize(Inboard_Points, Outboard_Points)
    % Calculate vectors through related points
    for i = 1:3         % Sort through columns
        for j = 1:6     % Sort through rows
            Vectors(j,i) = Outboard_Points(j,i) - Inboard_Points(j,i);
        end
    end
    
    % Unitize vectors
    for k = 1:3
        for m = 1:6
            Magnitude = sqrt((Vectors(m,1).^2) + (Vectors(m,2).^2) + (Vectors(m,3).^2));
            Vectorized(m,k) = Vectors(m,k) ./ Magnitude;
        end
    end

end

function MomentLengthed = MomentArmLength(Moment_Arm, Inboard_Points)
    for i = 1:3
        for j = 1:6
            MomentLengthed(j,i) = Inboard_Points(j,i) - Moment_Arm(i);
        end
    end
end

function A_CrossProducted = A_CrossProduct(Arm_Length, Unit_Vectors)
    for i = 1:6
        A_CrossProducted(i,1) = Arm_Length(i,2) .* Unit_Vectors(i,3) - Arm_Length(i,3) .* Unit_Vectors(i,2);
        A_CrossProducted(i,2) = Arm_Length(i,3) .* Unit_Vectors(i,1) - Arm_Length(i,1) .* Unit_Vectors(i,3);
        A_CrossProducted(i,3) = Arm_Length(i,1) .* Unit_Vectors(i,2) - Arm_Length(i,2) .* Unit_Vectors(i,1);
    end
end

function b_CrossProducted = b_CrossProduct(Arm_Vector, Vector)
    b_CrossProducted(1) = Arm_Vector(2) .* Vector(3) - Arm_Vector(3) .* Vector(2);
    b_CrossProducted(2) = Arm_Vector(3) .* Vector(1) - Arm_Vector(1) .* Vector(3);
    b_CrossProducted(3) = Arm_Vector(1) .* Vector(2) - Arm_Vector(2) .* Vector(1);
end

% Tire_Radius and Mz_Max not included as arguments: not used in these calcs?
function [Forces, G_Parameters] = AppliedForces(Loadcase_Number, Load_Distribution, Downforce_Distribution, ...
    Coeff_Drag, Coeff_Downforce, Frontal_Area, Air_Density, Total_Mass, CG_Height, Trackwidth, ...
    Wheelbase)
    
% Switch statement to determine G values and velocity for desired loadcase
    switch Loadcase_Number
        case 1
            Lateral_G = 0;
            Longitudinal_G = -1.4;
            Velocity = 30;
        case 2
            Lateral_G = 1.8;
            Longitudinal_G = 0;
            Velocity = 35;     
        case 3
            Lateral_G = 0;
            Longitudinal_G = 0;
            Velocity = 50;
        case 4
            Lateral_G = 0;
            Longitudinal_G = 0.98;
            Velocity = 35;
        case 5
            Lateral_G = 1.20;
            Longitudinal_G = -1.06;
            Velocity = 30;
        case 6
            Lateral_G = 0.94;
            Longitudinal_G = 0.52;
            Velocity = 35;
        case 7
            Lateral_G = 0;
            Longitudinal_G = 0; 
            Velocity = 35;
        case 8
            Lateral_G = 0;
            Longitudinal_G = -1.4;
            Velocity = 30;
        case 9
            Lateral_G = 1.6;
            Longitudinal_G = 0;
            Velocity = 30;
        case 10
            Lateral_G = 0;
            Longitudinal_G = 0;
            Velocity = 50;
    end

% Intermediate variables
Aero_Down = Coeff_Downforce * Velocity^2 * Load_Distribution / 2;
Aero_Drag = Coeff_Drag * Air_Density / 2 * Velocity^2 * Frontal_Area * Downforce_Distribution / 2;
WT_Lateral = Lateral_G * CG_Height / Trackwidth * Total_Mass * Load_Distribution;
WT_Long = Total_Mass * CG_Height  / Wheelbase * (-Longitudinal_G) * 0.5;
Fx_Car = Total_Mass * Longitudinal_G;
Fy_Car = Total_Mass * Lateral_G;

% Quasi-Intermediate variables
Fz = Aero_Down + WT_Lateral + WT_Long + Total_Mass * Load_Distribution * 0.5;
Fz_Percent = Fz / Total_Mass;

% Output variables
Fx = Fx_Car * Fz_Percent + Aero_Drag;
Fy = Fy_Car * Fz_Percent;

% Output arrays
Forces = [Fx Fy Fz];
G_Parameters = [Longitudinal_G Lateral_G Velocity];
  
end