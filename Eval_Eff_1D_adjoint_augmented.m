function [ abseff , gradient, varargout] = Eval_Eff_1D_adjoint_augmented(OptParm)
 
nout = max(nargout,1)-2;



OptParm = create_texture_profile(OptParm);
 
if nout ==0 
[t, E_f_x, E_f_y, E_f_z, top_angle, Eff] = Forward_sim(OptParm);
elseif nout == 1
[t, E_f_x, E_f_y, E_f_z, top_angle, Eff, n_eff] = Forward_sim(OptParm);
varargout{1} = n_eff;
elseif nout == 2
[t, E_f_x, E_f_y, E_f_z, top_angle, Eff, n_eff] = Forward_sim(OptParm);
varargout{1} = n_eff;
varargout{2} = t;
end
 
OptParm.theta0 = -top_angle;
[E_a_x, E_a_y, E_a_z] = Adjoint_sim(OptParm);
 
 
 
if strcmp(OptParm.gradient_type,'grayscale')
    dV = OptParm.height*OptParm.period/OptParm.N;
    if OptParm.polarization == -1
        grad_adj = dV*(-1j*conj(t)*(E_f_x.*E_a_x + E_f_z.*E_a_z));
    else
        grad_adj = dV*(1j*conj(t)*(E_f_y.*E_a_y));
    end
    grad_adj = grad_adj *(2*pi)/(OptParm.period*OptParm.wavelength);
    %  grad_adj = dV*real(-1j*exp(1j*angle(conj((t))))*(E_f_x.*E_a_x + E_f_z.*E_a_z)) * (n_Si^2 - n_Air^2);
    grad_adj = mean(grad_adj(2:end-1,:),1);
    grad_adj = reshape(grad_adj, 21, []);
    if size(OptParm.nvec,1)==1
        OptParm.nvec = OptParm.nvec';
    end
    grad_adj = real(mean(grad_adj, 1) .* transpose(OptParm.nvec)* 2* (OptParm.n_grating - OptParm.n_air));
    % grad_adj = (mean(grad_adj, 1)* (OptParm.n_grating^2 - OptParm.n_air^2));
    % grad_adj = 2*dV*real(conj(t)*(E_f_x.*(E_a_x) + E_f_z.*(E_a_z))) * (n_Si^2 - n_Air^2);
    
    T0 = abs(t)^2;
    abseff = T0;
    gradient = grad_adj;
else
    if OptParm.polarization == -1
        ME = MException('NOT IMPLEMENTED', 'TM polarization not yet implemented');
        throw(ME);
        
    else
        grad_adj = OptParm.height*(1j*conj(t)*(E_f_y.*E_a_y));
    end
    grad_adj = grad_adj *(2*pi)/(OptParm.period*OptParm.wavelength);
    grad_adj = mean(grad_adj(2:end-1,:),1);
    grad_adj = real(grad_adj*(OptParm.n_grating^2-OptParm.n_air^2))*(OptParm.period/(OptParm.N));
    T0 = abs(t)^2;
    abseff = T0;
    gradient = grad_adj;
    
end
 
end
 

