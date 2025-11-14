function [PVALV, PVARV, V0LV, V0RV,PVloop_area_LV,PVloop_area_RV,A_sol,B_sol,A_sol_rv,B_sol_rv]  = Calculate_PVA(init, params,fix_V0LV,fix_V0RV,V0LV_fixed,V0RV_fixed,coeffLV,coeffRV)
if isnan(coeffLV)
    volScales = [0.75 0.875 1.00 1.125 1.25];

    LVESP_points = zeros(length(volScales), 2);
    LVEDP_points = zeros(length(volScales), 2);

    RVESP_points = zeros(length(volScales), 2);
    RVEDP_points = zeros(length(volScales), 2);

    init0  = init;
    for change = 1:length(volScales)
        scale = volScales(change);

        % change the loading condition
        init.V_SV_s = init0.V_SV_s * scale;
        init.V_PV_s = init0.V_PV_s * scale;

        RunSimWithinPVA;
        if change == 3
            % ----------- LV ----------
            LVEDV = V_LV(1);
            LVEDP = P_LV(1);
            VLV = V_LV;
            PLV = P_LV;
            V0_Klotz_LV = LVEDV * (0.6 - 0.006 * LVEDP);
            V30_LV = V0_Klotz_LV  + ((LVEDV - V0_Klotz_LV) / ((LVEDP / 27.78)^(1/2.76)));
            V15_LV = 0.8 * (V30_LV - V0_Klotz_LV) + V0_Klotz_LV;
            if LVEDP <= 22
                beta_LV = log(LVEDP / 30) / log(LVEDV / V30_LV);
                alpha_LV = 30 / (V30_LV^beta_LV);
            else
                beta_LV = log(LVEDP / 15) / log(LVEDV / V15_LV);
                alpha_LV = LVEDP / (LVEDV^beta_LV);
            end
            % ----------- RV ----------
            RVEDV = V_RV(1);
            RVEDP = P_RV(1);
            VRV = V_RV;
            PRV = P_RV;
            V0_Klotz_RV = RVEDV * (0.6 - 0.006 * RVEDP);
            V30_RV = V0_Klotz_RV  + ((RVEDV - V0_Klotz_RV) / ((RVEDP / 27.78)^(1/2.76)));
            V15_RV = 0.8 * (V30_RV - V0_Klotz_RV) + V0_Klotz_RV;
            if RVEDP <= 22
                beta_RV = log(RVEDP / 30) / log(RVEDV / V30_RV);
                alpha_RV = 30 / (V30_RV^beta_RV);
            else
                beta_RV = log(RVEDP / 15) / log(RVEDV / V15_RV);
                alpha_RV = RVEDP / (RVEDV^beta_RV);
            end
        end
        LVESP_points(change,:) = [V_LV(Qa_pos_end), P_LV(Qa_pos_end)];
        LVEDP_points(change,:) = [V_LV(1), P_LV(1)];

        RVESP_points(change,:) = [V_RV(Qp_pos_end), P_RV(Qp_pos_end)];
        RVEDP_points(change,:) = [V_RV(1), P_RV(1)];
    end

    % ------- LV EDPVR 拟合 --------
    Ved_lv = LVEDP_points(:,1);
    Ped_lv = LVEDP_points(:,2);
    edpvr_fun = @(b, V) b(1) * (exp(b(2)*(V-b(3)))-1);
    init0_lv = [beta_LV, alpha_LV, V0_Klotz_LV];
    opts = optimset('MaxFunEvals', 2000, 'Display','off');
    fit_lv = lsqcurvefit(edpvr_fun, init0_lv, Ved_lv, Ped_lv, [], [], opts);

    A_fit_lv = fit_lv(1);
    B_fit_lv = fit_lv(2);
    V0_fit_lv = fit_lv(3);
    P0 = 1; % mmHg

    if fix_V0LV == 0
        V0LV = V0_fit_lv + (1/B_fit_lv)*log(P0/A_fit_lv + 1);
        if V0LV < 5
            V0LV = 5; % Set a minimum value for V0LV
        elseif V0LV > 0.6*min(VLV)
            V0LV = 0.6*min(VLV);
        end
    else
        V0LV = V0LV_fixed;
    end
    V1 = LVEDP_points(1,1); P1 = LVEDP_points(1,2);
    V3 = LVEDP_points(3,1); P3 = LVEDP_points(3,2);
    funA = @(A) (V3 - V0LV)*log(P1/A + 1) - (V1 - V0LV)*log(P3/A + 1);

    A_min = eps;
    A_max = min([P1, P3])-eps;
    try
        A_sol = fzero(funA, [A_min, A_max]);
    catch
        A_sol = A_fit_lv;
    end
    B_sol = log(P3/A_sol + 1)/(V3 - V0LV);

    PEDV = linspace(V0LV, LVEDV, 100);
    PEDP = A_sol * (exp(B_sol .* (PEDV - V0LV)) - 1);
    % ---- LV ESPVR ---
    % Follow Dan's suggestion, go to tangent
    k_LV = PLV(1:round(end/2))./(VLV(1:round(end/2))-V0LV);
    [~,LVESP_Loci] = max(k_LV);
    PESV = [VLV(LVESP_Loci) V0LV];
    PESP = [PLV(LVESP_Loci) 0];
    x_lv = [VLV(1:LVESP_Loci); PESV'; PEDV'];
    y_lv = [PLV(1:LVESP_Loci); PESP'; PEDP'];
    PVALV = polyarea(x_lv, y_lv);

    % ------- RV EDPVR  --------

    % 07/11/2025 comments: I don't know why RV fit is so bad.
    Ved_rv = RVEDP_points(:,1);
    Ped_rv = RVEDP_points(:,2);
    edpvr_fun_rv = @(b, V) b(1) * (exp(b(2)*(V-b(3)))-1);
    init0_rv = [beta_RV*10, alpha_RV, V0_Klotz_RV];
    fit_rv = lsqcurvefit(edpvr_fun_rv, init0_rv, Ved_rv, Ped_rv, [], [], opts);

    A_fit_rv = fit_rv(1);
    B_fit_rv = fit_rv(2);
    V0_fit_rv = fit_rv(3);

    if fix_V0RV == 0
        V0RV = V0_fit_rv + (1/B_fit_rv)*log(P0/A_fit_rv + 1);
        if V0RV < 1
            V0RV = 1;
        elseif V0RV > 0.6*min(VRV)
            V0RV = 0.6*min(VRV);
        end
    else
        V0RV = V0RV_fixed;
    end
    V1_rv = RVEDP_points(1,1); P1_rv = RVEDP_points(1,2);
    V3_rv = RVEDP_points(3,1); P3_rv = RVEDP_points(3,2);
    funA_rv = @(A) (V3_rv - V0RV)*log(P1_rv/A + 1) - (V1_rv - V0RV)*log(P3_rv/A + 1);

    A_min_rv = eps;
    A_max_rv = min([P1_rv, P3_rv])-eps;
    try
        A_sol_rv = fzero(funA_rv, [A_min_rv, A_max_rv]);
    catch
        A_sol_rv = A_fit_rv;
    end
    B_sol_rv = log(P3_rv/A_sol_rv + 1)/(V3_rv - V0RV);

    PEDV_rv = linspace(V0RV, RVEDV, 100);
    PEDP_rv = A_sol_rv * (exp(B_sol_rv .* (PEDV_rv - V0RV)) - 1);
    % ---- RV ESPVR ----

    k_RV = PRV(1:round(end/2))./(VRV(1:round(end/2))-V0RV);
    [~,RVESP_Loci] = max(k_RV);
    PESV_rv = [VRV(RVESP_Loci) V0RV];
    PESP_rv = [PRV(RVESP_Loci) 0];
    x_rv = [VRV(1:LVESP_Loci); PESV_rv'; PEDV_rv'];
    y_rv = [PRV(1:LVESP_Loci); PESP_rv'; PEDP_rv'];
    PVARV = polyarea(x_rv, y_rv);
else
    RunSimWithinPVA;
    LVEDV = V_LV(1);
    LVEDP = P_LV(1);
    VLV = V_LV;
    PLV = P_LV;
    A_sol = coeffLV(1);
    B_sol = coeffLV(2);
    V0LV = coeffLV(3);
    PEDV = linspace(V0LV, LVEDV, 100);
    PEDP = A_sol * (exp(B_sol .* (PEDV - V0LV)) - 1);
    k_LV = PLV(1:round(end/2))./(VLV(1:round(end/2))-V0LV);
    [~,LVESP_Loci] = max(k_LV);
    PESV = [VLV(LVESP_Loci) V0LV];
    PESP = [PLV(LVESP_Loci) 0];
    x_lv = [VLV(1:LVESP_Loci); PESV'; PEDV'];
    y_lv = [PLV(1:LVESP_Loci); PESP'; PEDP'];
    if x_lv(1) ~= x_lv(end) || y_lv(1) ~= y_lv(end)
        x_lv = [x_lv; x_lv(1)];
        y_lv = [y_lv; y_lv(1)];
    end
    PVALV = polyarea(x_lv, y_lv);

    RVEDV = V_RV(1);
    RVEDP = P_RV(1);
    VRV = V_RV;
    PRV = P_RV;
    A_sol_rv = coeffRV(1);
    B_sol_rv = coeffRV(2);
    V0RV = coeffRV(3);
    PEDV_rv = linspace(V0RV, RVEDV, 100);
    PEDP_rv = A_sol_rv * (exp(B_sol_rv .* (PEDV_rv - V0RV)) - 1);

    k_RV = PRV(1:round(end/2))./(VRV(1:round(end/2))-V0RV);
    [~,RVESP_Loci] = max(k_RV);
    PESV_rv = [VRV(RVESP_Loci) V0RV];
    PESP_rv = [PRV(RVESP_Loci) 0];
    x_rv = [VRV(1:RVESP_Loci); PESV_rv'; PEDV_rv'];
    y_rv = [PRV(1:RVESP_Loci); PESP_rv'; PEDP_rv'];
    if x_rv(1) ~= x_rv(end) || y_rv(1) ~= y_rv(end)
        x_rv = [x_rv; x_rv(1)];
        y_rv = [y_rv; y_rv(1)];
    end
    PVARV = polyarea(x_rv, y_rv);
end
% ----------- LV -----------
Nbeat_LV = round(length(VLV) / 2);  
PV_x_LV = VLV(1:Nbeat_LV);
PV_y_LV = PLV(1:Nbeat_LV);

if PV_x_LV(1) ~= PV_x_LV(end) || PV_y_LV(1) ~= PV_y_LV(end)
    PV_x_LV = [PV_x_LV; PV_x_LV(1)];
    PV_y_LV = [PV_y_LV; PV_y_LV(1)];
end

PVloop_area_LV = polyarea(PV_x_LV, PV_y_LV);

% ----------- RV -----------
Nbeat_RV = round(length(VRV) / 2);
PV_x_RV = VRV(1:Nbeat_RV);
PV_y_RV = PRV(1:Nbeat_RV);

if PV_x_RV(1) ~= PV_x_RV(end) || PV_y_RV(1) ~= PV_y_RV(end)
    PV_x_RV = [PV_x_RV; PV_x_RV(1)];
    PV_y_RV = [PV_y_RV; PV_y_RV(1)];
end
PVloop_area_RV = polyarea(PV_x_RV, PV_y_RV);
figure(888);plot(PV_x_RV, PV_y_RV);hold on;
plot(PV_x_LV, PV_y_LV);
plot(x_lv,y_lv);
plot(x_rv,y_rv);
end
