function [T, grad_adjs, varargout] = multi_level_beam_switch_objf(img, OptParm)
    nout = max(nargout,1)-2;
    
    OptParm.img = img;
OptParm = construct_OptParm(OptParm);
 
    num_level = length(OptParm.EFs);

 OptParm.ordernumberset = repmat(OptParm.diffraction_channels, [1 num_level]);
 
OptParm.n_grset= reshape(repmat(OptParm.n_grs, [length(OptParm.diffraction_channels) 1]),[1 length(OptParm.diffraction_channels)*num_level]);

%  for i =1 : 1000
     
OptParm_temp = OptParm;   






    parfor j=1:length(OptParm.diffraction_channels)*num_level
%         j = 1;
         [prv,vmax]=retio([],inf*1i);
         
   OptParm = struct();
   OptParm = OptParm_temp;
   OptParm.ordernumber = OptParm.ordernumberset(j);
   OptParm.n_gr = OptParm.n_grset(j);
         
         
if nout == 0
 [T0, grad_adj] = Eval_Eff_1D_adjoint_augmented(OptParm);
elseif nout == 1
[T0, grad_adj, n_eff] = Eval_Eff_1D_adjoint_augmented(OptParm);
n_effs(j,:) = n_eff;
elseif nout == 2
[T0, grad_adj, n_eff, t] = Eval_Eff_1D_adjoint_augmented(OptParm);
n_effs(j,:)  = n_eff;
ts(j) = t;
end

% end

% [T0, grad_adj, n_eff, t] = 


T(j) = T0;
grad_adjs(j,:) = grad_adj;

%      end
    end
    
 T = reshape(T, length(OptParm.diffraction_channels), num_level,[]).';
grad_adjs = reshape(grad_adjs, length(OptParm.diffraction_channels), num_level, []);
 grad_adjs = permute(grad_adjs, [2 1 3]);
%     grad_adjs = reshape(grad_adjs, num_level, length(OptParm.diffraction_channels),[]);
    
    
   if nout == 1
       n_effs = reshape(n_effs, length(OptParm.diffraction_channels), num_level,[]);
       n_effs = permute(n_effs, [2 1 3]);
varargout{1} = n_effs;
elseif nout == 2
       n_effs = reshape(n_effs, length(OptParm.diffraction_channels), num_level,[]);
       n_effs = permute(n_effs, [2 1 3]);
 ts = reshape(ts.',length(OptParm.diffraction_channels), num_level,  []).';

           
varargout{1} = n_effs;
varargout{2} = ts;
end

%    figure(1);
%    subplot(321);
%    plot(grad_adjs(1,:));
%    subplot(322);
%    plot(grad_adjs(2,:));
%    subplot(323);
%    plot(grad_adjs(3,:));
%    subplot(324);
%    plot(grad_adjs(4,:));
%    subplot(325);
%    plot(grad_adjs(5,:));
%    subplot(326);
%    plot(grad_adjs(6,:));
% F = T(1) + T(6);
% G = grad_adjs(1,:) + grad_adjs(6,:);

% disp(img);
 end