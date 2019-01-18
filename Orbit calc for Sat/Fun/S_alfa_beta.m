function S_alfa_beta_ = S_alfa_beta(alfa, beta)
    % alfa [rad]
    % beta [rad]
    
    if abs(alfa) > (pi/180) & abs(alfa) < (179*pi/180) & abs(beta) > (pi/180) & abs(beta) < (179*pi/180)
        Sx_ = (1/tan(alfa))/(1+(1/(tan(alfa)*tan(alfa)))+(1/(tan(beta)*tan(beta))))^0.5;
        Sy_ = (1/tan(beta))/(1+(1/(tan(alfa)*tan(alfa)))+(1/(tan(beta)*tan(beta))))^0.5;
        Sz_ = 1/(1+(1/(tan(alfa)*tan(alfa)))+(1/(tan(beta)*tan(beta))))^0.5;
        S_alfa_beta_ = [Sx_ Sy_ Sz_];
    else
        S_alfa_beta_ = [1 0 0];        
    end
end