function [sub] = MEG_num2Sub(num)
% convert subject label to a subject number

switch num
    case 1
        sub = 'meg_081913';
    case 2
        sub = 'meg_082013';
    case 3
        sub = 'meg_082113';
end
