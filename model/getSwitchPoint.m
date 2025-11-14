function [sw_b, sw_d] = getSwitchPoint(D,subjectNum)

    k_min = readParam(subjectNum,'k_min');
    k_max = readParam(subjectNum,'k_max');
    b = readParam(subjectNum,'b');
    d = readParam(subjectNum,'d');

    load('parameta_b.mat');
    beta = betaSol;
    
    b1 = beta(1);
    b2 = beta(2);
    b3 = beta(3);
    b4 = beta(4);
    b5 = beta(5);
    b6 = beta(6);
    
    r_b = b1*exp(b2*k_max) + b3*exp(b4*k_min) + b5*b + b6;
    
    
    load('parameta_d.mat');
    beta = betaSol;

    b1 = beta(1);
    b2 = beta(2);
    b3 = beta(3);
    b4 = beta(4);
    b5 = beta(5);
    b6 = beta(6);
    
    r_d = b1*exp(b2*k_max) + b3*exp(b4*k_min) + b5*d + b6;
    
    sw_b = b*D - r_b*D*(1-b);
    sw_d = d*D + r_d*D*d;
end