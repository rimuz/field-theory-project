% Compendium page 125
function [S11_LR, S13_LR, S31_LR, S33_LR] = ...
        combineLR(S11_L, S12_L, S21_L, S22_L, S22_R, S23_R, S32_R, S33_R)
    
    N = length(S22_R);
    S11_LR = S11_L + S12_L*inv(eye(N) - S22_R*S22_L)*S22_R*S21_L;
    S13_LR = S12_L*inv(eye(N) - S22_R*S22_L)*S23_R;
    S31_LR = S32_R*inv(eye(N) - S22_L*S22_R)*S21_L;
    S33_LR = S33_R + S32_R*inv(eye(N) - S22_L*S22_R)*S22_L*S23_R;
end