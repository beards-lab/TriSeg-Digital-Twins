function [Error,Amref, Vw, dimensions] = geom_0(LVEDV, RVEDV, left, use_Vw_LV, right, use_Vw_RV, LvSepR,inputs)
% getting a better guess for Amref_RV instead of assuming RV is a sphere,
% which it isn't. Will replace calculation of Amref_LV and Sep and well. 

% This function will calculate the initial geometry of the heart. i'm not
% sure how this will work with calc_xm_ym, I still don't understand that
% function

% left and right hold the information that we know about the walls of the
% left side (lw and sw) and the right side. can be wall volume or thickness

%% initial spherical calculations (need LVEDV and wall thicknesses)
% LV lumen radius (cm)
r_LV_and_SEP = (LVEDV * 3 / (4* pi))^(1/3); 
if(~use_Vw_LV)
    H_LW_and_SW = left; % if used echo info
    r_o_LV_and_SEP = r_LV_and_SEP + H_LW_and_SW; 
    V_o_LV_and_SEP = 4/3 * pi * r_o_LV_and_SEP^3; 
    Vw_LV_and_SEP = V_o_LV_and_SEP - LVEDV;
    Vw_LV  = Vw_LV_and_SEP * LvSepR; % 
    Vw_SEP = Vw_LV_and_SEP * (1 - LvSepR); 
else
    Vw_LV_and_SEP = left; % if used mri info
    Vw_LV = LvSepR * Vw_LV_and_SEP;
    Vw_SEP = Vw_LV_and_SEP - Vw_LV;
    H_LW_and_SW = (r_LV_and_SEP^3 + (3/(4*pi))*Vw_LV_and_SEP)^(1/3) - r_LV_and_SEP;
end
r_m_LV_and_SEP = ((1/2)*((r_LV_and_SEP + H_LW_and_SW)^3 + r_LV_and_SEP^3))^(1/3); % radius dividing Vw into shells of equal volume
%% Set axial and radial axes based on LvSepR
% Find a beta such that LvSepR is satisfied. This assumes midwall surface
% is shared between LV and SEP
beta = acos(2 * LvSepR - 1);
ym = r_m_LV_and_SEP * sin(beta);
xm_SEP = r_m_LV_and_SEP - (r_m_LV_and_SEP^2 - ym^2)^(1/2);
xm_LV = xm_SEP - 2*r_m_LV_and_SEP; % right of radial axis is defined to be negative. 

Vm_SEP = (pi / 6) * xm_SEP * (xm_SEP^2 + 3 * ym^2); 
SEP_CAP = Vm_SEP + 0.5 * Vw_SEP; % Volume of spherical cap of left heart to the left of the axial axis



%% Right side
syms h
if(~use_Vw_RV) % if you are calculating RV geometry based on H_RW
    if(length(right) == 1)
        H_RW = right;
        Vw_RV = (2 * pi / 3) * (2*ym^2 / (h^2 + ym^2)) * ...
            ((-(h^2 + ym^2) / (2*h) + (1/2)*(((H_RW^6 + 16*(-(h^2 + ym^2) / (2*h))^6)^(1/2) + 4*(-(h^2 + ym^2) / (2*h))^3)^(1/3) ...
            - (H_RW^2) / (((H_RW^6 + 16*(-(h^2 + ym^2) / (2*h))^6)^(1/2) + 4*(-(h^2 + ym^2) / (2*h))^3)^(1/3)) ...
            + H_RW - 2*(-(h^2 + ym^2) / (2*h))))^3 ...
            - (-(h^2 + ym^2) / (2*h) + (1/2)*(((H_RW^6 + 16*(-(h^2 + ym^2) / (2*h))^6)^(1/2) + 4*(-(h^2 + ym^2) / (2*h))^3)^(1/3) ...
            - (H_RW^2) / (((H_RW^6 + 16*(-(h^2 + ym^2) / (2*h))^6)^(1/2) + 4*(-(h^2 + ym^2) / (2*h))^3)^(1/3)) ...
            + H_RW - 2*(-(h^2 + ym^2) / (2*h))) - H_RW)^3);
        % equation definition
        eqn = RVEDV == (pi / 6) * (-ym^2 / h)  * ((-ym^2 / h)^2  + 3 * ym^2) - 0.5 * Vw_RV - SEP_CAP; 
        h = vpasolve(eqn,h,[-Inf 0]); % h < 0. put restriction on h_RV. I don't know how to set bounds as -Inf < h < 0
        if(~use_Vw_RV)
            Vw_RV = double(subs(Vw_RV)); %puts sym variable h that was just calculated into previous Vw_RV expression
        end
        h_RV = double(h);% convert from sym type to double type
        assert(h_RV < 0);
        
        % Extract outputs
        xm_RV = -ym^2 / h_RV;
        Error = NaN;
    elseif(length(right) == 2)
        P_RV = right(1);
        k_pas_RV_est = right(2);
        sigma_pas = 2.1585e-4;
        syms xm_RV Vw_RV

        eq1 = abs(k_pas_RV_est - P_RV / ((-2/3 * ((3 * ((-2 * xm_RV) / (xm_RV^2 + ym^2)) * Vw_RV) / (2 * (pi * (xm_RV^2 + ym^2)))) * (1 + ((3 * ((-2 * xm_RV) / (xm_RV^2 + ym^2)) * Vw_RV) / (2 * (pi * (xm_RV^2 + ym^2))))^2 / 3 + ((3 * ((-2 * xm_RV) / (xm_RV^2 + ym^2)) * Vw_RV) / (2 * (pi * (xm_RV^2 + ym^2))))^4 / 5)) * sigma_pas))/100000;
        eq2 = abs(RVEDV - (pi / 6) * (xm_RV)  * (xm_RV^2  + 3 * ym^2) - 0.5 * Vw_RV - SEP_CAP);
        eq1_f = matlabFunction(eq1, 'Vars', [Vw_RV, xm_RV]);
        eq2_f = matlabFunction(eq2, 'Vars', [Vw_RV, xm_RV]);
        objective = @(x) [eq1_f(x(1),x(2)); eq2_f(x(1), x(2))];
        BSA = sqrt((inputs.Height * inputs.Weight) / 3600);
        if inputs.Sex == 1
            x0 = [62.5592 6.5731]; %Vw, xm
            lb = x0(1)/1.9 *0.8363* BSA; 
            ub = x0(1)/1.9 *2.4330* BSA; %Vw, based on male ub: 240 from Andrew's data by linear
            % regresstion VwMRI = 0.6342VwTTE-0.003236. the result should be 152.2048
            % and scalar is 152.2048*1.055/66 = 2.4330
        else
            x0 = [45.4976 5.8212];
            lb = x0(1)/1.6 *0.8363* BSA;  %Vw, based on female lb 60 from Andrew's data by linear
            % regresstion VwMRI = 0.6342VwTTE-0.003236. the result should be 152.2048
            % and scalar is 38.0488*1.055/48 = 0.8363
            ub = x0(1)/1.6 *2.4330* BSA; %Vw, xm
        end
        % lb = [30,4];  %Vw, xm
        % ub = [200, 12]; %Vw, xm
        % rng default % Reproducible initial point
        % opts = optimoptions(@fmincon,'Algorithm','interior-point','Display','iter');
        % fminconstr = @(x) deal([], objective(x));
        % x = fmincon(@(x)0,x0,[],[],[],[],lb,ub,fminconstr,opts);       
        % Vw_RV = x(1);
        % xm_RV = x(2);
        opts = optimoptions('fsolve','Display','none',...
            'MaxIter',4e4,'MaxFunctionEvaluations',4e4,'Algorithm','levenberg-marquardt');
       [solution,E] = fsolve(objective,x0,opts);
        Error = E(1)^2+E(2)^2;
        Vw_RV = double(solution(1));
       
        if Vw_RV < lb
            Vw_RV = lb;
        end
        if Vw_RV > ub
            Vw_RV = ub;
        end
        xm_RV = double(solution(2));
    end

else % if you are calculating RV geometry based on Vw_RV
    assert(use_Vw_RV);
    Vw_RV = right;
    % equation definition
    eqn = RVEDV == (pi / 6) * (-ym^2 / h)  * ((-ym^2 / h)^2  + 3 * ym^2) - 0.5 * Vw_RV - SEP_CAP; 
    h = vpasolve(eqn,h,[-Inf 0]); % h < 0. put restriction on h_RV. I don't know how to set bounds as -Inf < h < 0
    if(~use_Vw_RV)
        Vw_RV = double(subs(Vw_RV)); %puts sym variable h that was just calculated into previous Vw_RV expression
    end
    h_RV = double(h);% convert from sym type to double type
    assert(h_RV < 0);
    
    % Extract outputs
    xm_RV = -ym^2 / h_RV;
    Error = NaN;
end

% if(~use_Vw_RV)
%     r_m_RV = (xm_RV - h_RV) / 2; 
%     %^ radius of RV spherical cap is c (x-coord of center) minus h_RV (distance of right point of whole RV sphere on axial axis) 
%     delta_rw = 0.5 * (((H_RW^6 + 16*r_m_RV^6)^(1/2) + 4*r_m_RV^3)^(1/3) - (H_RW^2 / ((H_RW^6 + 16*r_m_RV^6)^(1/2) + 4*r_m_RV^3)^(1/3)) + H_RW - 2*r_m_RV);
%     Vw_RV = (2 * pi / 3)*((2 * h_RV^2) / (h_RV^2 + ym^2))*((r_m_RV + delta_rw)^3 - (r_m_RV + delta_rw - H_RW)^3);
% end


%% Outputs
Amref_LV  = pi * (xm_LV^2  + ym^2);
Amref_SEP = pi * (xm_SEP^2  + ym^2);
Amref_RV  = pi * (xm_RV^2  + ym^2);

Amref = [Amref_LV, Amref_SEP, Amref_RV];
Vw = [Vw_LV, Vw_SEP, Vw_RV];
dimensions = [xm_LV, xm_SEP, xm_RV, ym]; % initial conditions for calc_xm_ym
end