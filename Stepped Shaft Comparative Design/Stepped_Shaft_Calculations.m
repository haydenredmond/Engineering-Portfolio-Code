%% Stepped Shaft Comparative Calculations
% Author: Hayden Redmond
clc; clear; close all

%% Shaft Physical Properties
Minor_Diameter = [1.25 1.75 2.25]; % Diameter of smaller part of shaft
Major_Diameter = [2 2.5 3]; % Diameter of larger part of shaft
Radius = [(Major_Diameter(1)-Minor_Diameter(1))/2 (Major_Diameter(2)-Minor_Diameter(2))/2 (Major_Diameter(3)-Minor_Diameter(3))/2]; % Radius of fillet between sections
Minor_Length = [2 2.5 3]; % Length of part of shaft with smaller diameter
Major_Length = [2.5 3 3.5]; % Length of part of shaft with larger diameter
Surface_Finish = 'Machined';

%% Material Properties (psi)
% Titanium alloy (Ti-6Al-4V)
Ti_Sy = 1.28e5; 
Ti_Sut = 1.38e5;
Ti_Suc = 1.41e5;
% Steel (4130)
St_Sy = 6.31e4;
St_Su = 9.72e4;
% Aluminum (7075-T6)
Al_Sy = 7.3e4;
Al_Su = 8.3e4;

%% Forces
F_Max = 14000; % lbf - Completely reversed
T_Max = 500; % lb ft - Completely reversed

%% Executable
% Prompt user to select stress concentration factors (K_t,K_ts) for each geometry set
for i = 1:3
    K_ts(i) = Stress_Conc_Factor('Torsion', Major_Diameter(i)/Minor_Diameter(i), Radius(i)/Minor_Diameter(i));
end

for j = 1:3
    K_t(j) = Stress_Conc_Factor('Bending', Major_Diameter(j)/Minor_Diameter(j), Radius(j)/Minor_Diameter(j));
end

% Calculate principal stresses for each geometry set
for IndexPrincipal = 1:3
    [Sigma_A_in Sigma_B_in] = PrincipalStresses(T_Max,F_Max,Minor_Diameter(IndexPrincipal),Minor_Length(IndexPrincipal),K_t(IndexPrincipal),K_ts(IndexPrincipal));
    All_Principal_Stresses(IndexPrincipal,1) = Sigma_A_in;
    All_Principal_Stresses(IndexPrincipal,2) = Sigma_B_in;
end

% Calculate (yielding) factors of safety based on failure theories
for YieldTheoryIndex = 1:3
    % Compute SF using Ductile Coulomb-Mohr for all relevant sets
    n_DCM_Ti(YieldTheoryIndex) = DuctileCoulombMohr(All_Principal_Stresses(YieldTheoryIndex,1),All_Principal_Stresses(YieldTheoryIndex,2),Ti_Sut,Ti_Suc); 
    % Compute SF using Distortion-Energy for all sets
    n_DE_St(YieldTheoryIndex,1) = DistortionEnergy(All_Principal_Stresses(YieldTheoryIndex,1),All_Principal_Stresses(YieldTheoryIndex,2),St_Sy); % DE for Steel
    n_DE_Al(YieldTheoryIndex,1) = DistortionEnergy(All_Principal_Stresses(YieldTheoryIndex,1),All_Principal_Stresses(YieldTheoryIndex,2),Al_Sy); % DE for Aluminum
    % Compute SF using Maximum Shear Stress for all sets
    n_MSS_St(YieldTheoryIndex,1) = MaximumShearStress(All_Principal_Stresses(YieldTheoryIndex,1),All_Principal_Stresses(YieldTheoryIndex,2),St_Sy); % MSS for Steel
    n_MSS_Al(YieldTheoryIndex,1) = MaximumShearStress(All_Principal_Stresses(YieldTheoryIndex,1),All_Principal_Stresses(YieldTheoryIndex,2),Al_Sy); % MSS for Aluminum
end

% Call BendingStress to get alternating and mean stress
for BendingStressIndex = 1:3
    [Sigma_x(BendingStressIndex),Sigma_Alternating(BendingStressIndex),Sigma_Mean(BendingStressIndex)] = BendingStress(F_Max,Minor_Diameter(BendingStressIndex),Minor_Length(BendingStressIndex),K_t(BendingStressIndex));
end

% Call FatigueSF to get (fatigue) factors of safety based on infinite life
for FatigueSFIndex = 1:3
    % Titanium
    [FatigueFOS_Ti(FatigueSFIndex,1),FatigueFOS_Ti(FatigueSFIndex,2),FatigueFOS_Ti(FatigueSFIndex,3),FatigueFOS_Ti(FatigueSFIndex,4)] = FatigueSF(Ti_Sut,Ti_Sy,Sigma_Alternating(FatigueSFIndex),Sigma_Mean(FatigueSFIndex),Surface_Finish,Minor_Diameter(FatigueSFIndex),'Titanium',FatigueSFIndex);
    % Steel
    [FatigueFOS_St(FatigueSFIndex,1),FatigueFOS_St(FatigueSFIndex,2),FatigueFOS_St(FatigueSFIndex,3),FatigueFOS_St(FatigueSFIndex,4)] = FatigueSF(St_Su,St_Sy,Sigma_Alternating(FatigueSFIndex),Sigma_Mean(FatigueSFIndex),Surface_Finish,Minor_Diameter(FatigueSFIndex),'Steel',FatigueSFIndex);
    % Aluminum
    [FatigueFOS_Al(FatigueSFIndex,1),FatigueFOS_Al(FatigueSFIndex,2),FatigueFOS_Al(FatigueSFIndex,3),FatigueFOS_Al(FatigueSFIndex,4)] = FatigueSF(Al_Su,Al_Sy,Sigma_Alternating(FatigueSFIndex),Sigma_Mean(FatigueSFIndex),Surface_Finish,Minor_Diameter(FatigueSFIndex),'Aluminum',FatigueSFIndex);
end

% Preallocating estimated cycle tables for clarity
Cycles_Ti(1:3,1:4) = 67;    % Default value to show cycles weren't calculated
Cycles_St(1:3,1:4) = 67;    % Default value to show cycles weren't calculated
Cycles_Al(1:3,1:4) = 67;    % Default value to show cycles weren't calculated

% Prompt user for fatigue strength fraction 
f_Ti = Fatigue_Factor(Ti_Sut);
f_St = Fatigue_Factor(St_Su);
f_Al = Fatigue_Factor(Al_Su);

% Check if any fatigue factors of safety < 1, call CycleLife for these
for LifeGeomIndex = 1:3
    for LifeTheoryIndex = 1:4
        if FatigueFOS_Ti(LifeGeomIndex,LifeTheoryIndex) < 1
            % Call CycleLife
            Cycles_Ti(LifeGeomIndex,LifeTheoryIndex) = CycleLife(Sigma_Alternating(LifeGeomIndex),Ti_Sut,Surface_Finish,Minor_Diameter(LifeGeomIndex),f_Ti);
        end
        
        if FatigueFOS_St(LifeGeomIndex,LifeTheoryIndex) < 1
            % Call CycleLife
            Cycles_St(LifeGeomIndex,LifeTheoryIndex) = CycleLife(Sigma_Alternating(LifeGeomIndex),St_Su,Surface_Finish,Minor_Diameter(LifeGeomIndex),f_St);
        end

        if FatigueFOS_Al(LifeGeomIndex,LifeTheoryIndex) < 1
            % Call CycleLife
            Cycles_Al(LifeGeomIndex,LifeTheoryIndex) = CycleLife(Sigma_Alternating(LifeGeomIndex),Al_Su,Surface_Finish,Minor_Diameter(LifeGeomIndex),f_Al);
        end        
    end
end

%% Functions

% Stress_Conc_Factor user prompt (CALLED)
function K_Value = Stress_Conc_Factor(type, D_over_d, r_over_d)

    labelStr = sprintf('D/d = %.4f, r/d = %.4f:', D_over_d, r_over_d);

    figWidth  = 900;
    figHeight = 450 + 120;
    fig = uifigure('Name', 'Enter Stress Concentration Factor', ...
        'Position', [300 200 figWidth figHeight]);

    if strcmpi(type,'Torsion')
        imgSrc = 'Figure A-15-8.png';
    elseif strcmpi(type,'Bending')
        imgSrc = 'Figure A-15-9.png';
    else
        error('Type must be Torsion or Bending');
    end

    uiimage(fig, 'ImageSource', imgSrc, 'Position', [0 120 900 450]);

    uilabel(fig, 'Text', labelStr, ...
        'Position', [20 70 350 30], 'FontSize', 14);

    % Numeric text box (empty by default)
    txtBox = uieditfield(fig, 'numeric', ...
        'Position', [380 70 450 30], ...
        'FontSize', 14, ...
        'Value', [], ...
        'AllowEmpty', true);

    % OK button
    uibutton(fig, 'Text', 'OK', ...
        'Position', [400 20 100 30], ...
        'ButtonPushedFcn', @(btn,event) closeAndReturn(fig));

    % Figure KeyPressFcn to detect Enter
    fig.KeyPressFcn = @(src,event) keyPressCallback(src, event, txtBox, fig);

    % Wait for user input
    uiwait(fig);

    % Get value after OK pressed
    if isvalid(fig)
        K_Value = txtBox.Value;  % directly read from the textbox
        delete(fig);
    else
        K_Value = [];
    end

end

% Fatigue_Factor user prompt (CALLED)
function f = Fatigue_Factor(S_ut) 
    labelStr = sprintf('S_ut = %.2f psi:', S_ut);

    figWidth  = 900;
    figHeight = 450 + 120;
    fig = uifigure('Name', 'Enter Fatigue Strength Fraction:', ...
        'Position', [300 200 figWidth figHeight]);

    imgSrc = 'Figure 6-23.png';

    uiimage(fig, 'ImageSource', imgSrc, 'Position', [0 120 900 450]);

    uilabel(fig, 'Text', labelStr, ...
        'Position', [20 70 350 30], 'FontSize', 14);

    % Numeric text box (empty by default)
    txtBox = uieditfield(fig, 'numeric', ...
        'Position', [380 70 450 30], ...
        'FontSize', 14, ...
        'Value', [], ...
        'AllowEmpty', true);

    % OK button
    uibutton(fig, 'Text', 'OK', ...
        'Position', [400 20 100 30], ...
        'ButtonPushedFcn', @(btn,event) closeAndReturn(fig));

    % Figure KeyPressFcn to detect Enter
    fig.KeyPressFcn = @(src,event) keyPressCallback(src, event, txtBox, fig);

    % Wait for user input
    uiwait(fig);

    % Get value after OK pressed
    if isvalid(fig)
        f = txtBox.Value;  % directly read from the textbox
        delete(fig);
    else
        f = [];
    end

end

function closeAndReturn(fig) % (NOT CALLED IN EXEC)
    uiresume(fig);  % resume uiwait
end

function keyPressCallback(~, event, txtBox, fig) % (NOT CALLED IN EXEC)
    % Only act if Enter is pressed AND numeric box is focused
    if strcmp(event.Key, 'return') && (gco == txtBox)
        closeAndReturn(fig);
    end
end

% BendingStress (CALLED)
function [Sigma_x,Sigma_Alternating,Sigma_Mean] = BendingStress(F_Max,Critical_Diameter,Minor_Length,K_t)
    % MOI
    I = (pi/64)*Critical_Diameter^4;

    % Max Moment, M_Max
    M_Max = F_Max * Minor_Length;
    
    % Max bending stress, Sigma_x
    Sigma_Max = K_t*(M_Max * (Critical_Diameter/2))/I;
    Sigma_x = Sigma_Max;
    
    % Min bending stress, Sigma
    Sigma_Min = 0;

    % Alternating
    Sigma_Alternating = abs((Sigma_Max-Sigma_Min)/2);           % DOES NOT INCLUDE FATIGUE STRESS CONCENTRATION
    
    % Mean
    Sigma_Mean = (Sigma_Max+Sigma_Min)/2;                       % DOES NOT INCLUDE FATIGUE STRESS CONCENTRATION

end

% TorisionalStress (NOT CALLED IN EXEC)
function [Tau,Tau_Alternating] = TorsionalStress(Max_Torque,Critical_Diameter,K_ts)
    % MOI
    J = (pi/32)*Critical_Diameter^4;
    % Tau (Not the same as Tau_Max)
    Tau = K_ts*(Max_Torque*(Critical_Diameter/2))/J;
    % Tau_Alternating 
    Tau_Alternating = 1;        % PLACEHOLDER
end

% PrincipalStresses (CALLED)
function [Sigma_A,Sigma_B] = PrincipalStresses(Max_Torque,F_Max,Critical_Diameter,Minor_Length,K_t,K_ts)
    % Call TorsionalStress
    [Tau,Tau_Alternating] = TorsionalStress(Max_Torque,Critical_Diameter,K_ts);
    % Call BendingStress
    [Sigma_x,Sigma_Alternating,Sigma_Mean] = BendingStress(F_Max,Critical_Diameter,Minor_Length,K_t); 
    % Calculate principal stresses
    Sigma_A = ((Sigma_x+0)/2) + sqrt(((Sigma_x-0)/2)^2+Tau^2);
    Sigma_B = ((Sigma_x+0)/2) - sqrt(((Sigma_x-0)/2)^2+Tau^2);
end

% --------------------------- Yielding ---------------------------

% DuctileCoulombMohr (CALLED)
function n_DCM = DuctileCoulombMohr(Sigma_A,Sigma_B,S_yt,S_yc)
    % (only for Titanium bc it's Syt =/= Syc)
    if Sigma_B <= Sigma_A && Sigma_B >= 0
        n_DCM = S_yt / Sigma_A;
    elseif 0 <= Sigma_A && 0 >= Sigma_B
        n_DCM = ((Sigma_A/S_yt)-(Sigma_B/S_yc))^(-1);
    elseif Sigma_A <= 0 && Sigma_A >= Sigma_B
        n_DCM = S_yc / Sigma_B;
    else
        error('Principal Stress Mismatch Error')
    end
end
    
% DistortionEnergy (CALLED)
function n_DE = DistortionEnergy(Sigma_A,Sigma_B,S_y)
    % (NOT for titanium bc it's Syt =/= Syc)
    Sigma_Prime = sqrt(Sigma_A^2-Sigma_A*Sigma_B+Sigma_B^2);
    n_DE = S_y / Sigma_Prime;
end

% MaximumShearStress (CALLED)
function n_MSS = MaximumShearStress(Sigma_A,Sigma_B,S_y)
    % (NOT for titanium bc it's Syt =/= Syc)
    if Sigma_B <= Sigma_A && Sigma_B >= 0
        Sigma_1 = Sigma_A;
        Sigma_3 = 0;
    elseif 0 <= Sigma_A && 0 >= Sigma_B
        Sigma_1 = Sigma_A;
        Sigma_3 = Sigma_B;
    elseif Sigma_A <= 0 && Sigma_A >= Sigma_B
        Sigma_1 = 0;
        Sigma_3 = Sigma_B;
    else
        error('Principal Stress Mismatch Error')
    end
    n_MSS = S_y / (Sigma_1-Sigma_3);
end

% --------------------------- Fatigue ---------------------------

% EnduranceLimit (NOT CALLED IN EXEC)
function [S_e, S_e_Prime] = EnduranceLimit(S_ut,Surface_Finish,Critical_Diameter)
    % Endurance limit (Se')
    if S_ut <= 200e3
        S_e_Prime = 0.5*S_ut;
    elseif S_ut > 200e3
        S_e_Prime = 100e3;
    else
        error('Endurance Limit Mismatch Error')
    end
    % Surface finish, k_a
    switch Surface_Finish
        case 'Ground'
            a = 1.21;
            b = -0.067;
        case 'Machined'
            a = 2;
            b = -0.217;
        case 'Hot-rolled'
            a = 11;
            b = -0.65;
        case 'As-forged'
            a = 12.7;
            b = -0.758;
        otherwise
            error('Surface Finish Naming Error')
    end
    k_a = a * (S_ut*10^3)^b;

    % Size factor, k_b
    if Critical_Diameter >= 0.3 && Critical_Diameter <= 2
        k_b = 0.879*Critical_Diameter^(-0.107);
    elseif Critical_Diameter > 2 && Critical_Diameter <= 10
        k_b = 0.91*Critical_Diameter^(-0.157);
    else
        error('Critical diameter out of bounds for size factor k_b')
    end

    S_e = k_a*k_b*S_e_Prime;
end   
   
% FatigueSF (CALLED)
function [Goodman,Soderberg,ASME,Gerber] = FatigueSF(S_ut,S_y,Sigma_Alternating,Sigma_Mean,Surface_Finish,Critical_Diameter,Material,Geom) 
    [S_e S_e_Prime] = EnduranceLimit(S_ut,Surface_Finish,Critical_Diameter);
    
    Not_Predicted_Msg = "Infinite life not predicted with";
    Space = " ";
    GeomSet = "Geometry Set";

    % Check if Sigma_Mean is less than zero (should be >= 0)
    if Sigma_Mean < 0
        error('MEAN STRESS IS LESS THAN ZERO!')
    end

    % Goodman
    Goodman = ((Sigma_Alternating/S_e)+(Sigma_Mean/S_ut)).^(-1);
    if Goodman < 1
        Criteria = "Goodman:";
        disp(Not_Predicted_Msg + Space + Criteria + Space + Material + Space + GeomSet + Space + Geom)
    end

    % Soderberg
    Soderberg = ((Sigma_Alternating/S_e)+(Sigma_Mean/S_y)).^(-1);
    if Soderberg < 1
        Criteria = "Soderberg:";
        disp(Not_Predicted_Msg + Space + Criteria + Space + Material + Space + GeomSet + Space + Geom)
    end

    % ASME Elliptical
    ASME = ((Sigma_Alternating/S_e).^2 + (Sigma_Mean/S_y).^2).^(-1/2);
    if ASME < 1
        Criteria = "ASME Elliptical:";
        disp(Not_Predicted_Msg + Space + Criteria + Space + Material + Space + GeomSet + Space + Geom)
    end

    % Gerber
    Gerber = (1/2)*(S_ut/Sigma_Mean).^2 * (Sigma_Alternating/S_e) * (-1+sqrt(1+((2*Sigma_Mean*S_e)/(S_ut*Sigma_Alternating)).^2));
    if Gerber < 1
        Criteria = "Gerber:";
        disp(Not_Predicted_Msg + Space + Criteria + Space + Material + Space + GeomSet + Space + Geom)
    end

end

% CycleLife (CALLED)
function N = CycleLife(Sigma_Alternating,S_ut,Surface_Finish,Critical_Diameter,f)
    % FOR USE WHEN INFINITE LIFE IS NOT PREDICTED
    % Call EnduranceLimit
    [S_e,S_e_Prime] = EnduranceLimit(S_ut,Surface_Finish,Critical_Diameter);
    
    % Convert S_ut units from psi to kpsi
    S_ut = S_ut * 10^(-3);

    % Calculate constants a and b
    a = ((f*S_ut)^2)/S_e_Prime;
    b = (-1/3)*log10((f*S_ut)/S_e_Prime);

    % Calculate Number of cycles
    N = (Sigma_Alternating/a)^(1/b);
end
