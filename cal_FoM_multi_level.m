function [F,G] = cal_FoM_multi_level(T, grad_adjs, OptParm) 

if length(OptParm.EFs) == 2
   
   [F,G] = cal_FoM(T, grad_adjs, OptParm.b_coefficient);
    
elseif length(OptParm.EFs) == 3
    
     [F,G] = cal_FoM_three_level(T, grad_adjs, OptParm.b_coefficient);
    
elseif length(OptParm.EFs) == 4
    
    [F,G] = cal_FoM_four_level(T, grad_adjs);
    
end






end