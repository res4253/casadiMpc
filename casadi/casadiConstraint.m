function [formula, const] = constraint(t, x, u)  
    % example
    %
    % x(1)+2x(2) = 0, 1 <= x(3)-x(4) <= 4;　なら
    % 
    % formula = [x(1)+2x(2)
    %            x(3)-x(4)]
    % const   = [0, 0
    %            1, 4]
    % 時変の制約も可能です
    
    formula = [];
    const = [];
end
