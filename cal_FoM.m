function [F,G] = cal_FoM(T, grad_adjs, b_coefficient)

% 
% inds = [1 3];
% 
% for i=1:size(T,1) 
%     FoM_deno = sum(T(i,:),2);
%     FoM_nume = T(i,inds(i));
% 
%     grad_deno = sum(grad_adjs(i,:,:),2);
%     grad_nume = grad_adjs(i,inds(i),:);
% 
% %    F(i) = (FoM_nume/FoM_deno);
% % G(:,i) =  ( (grad_nume*FoM_deno - FoM_nume*grad_deno)/ FoM_deno^2 );
% 
%   D(i) = (FoM_nume/FoM_deno);
% D_adjs(:,i) =  ( (grad_nume*FoM_deno - FoM_nume*grad_deno)/ FoM_deno^2 );
% DE(i) = FoM_nume;
% DE_adjs(:,i) = grad_nume;
% 
% % F(i) = D(i)*DE(i);
% % G(:,i) = D_adjs(:,i)*DE(i) + D(i)*DE_adjs(:,i);
% 
% 
% %   F(i) = FoM_nume;
% % G(:,i) = grad_nume;
% 
% F(i) = D(i);
% G(:,i) = D_adjs(:,i);
% 
% end

F = -sum(T,'all') + (1+b_coefficient)*(T(1,1)+T(2,3));
G = squeeze(-sum(sum(grad_adjs,1),2) +(1+b_coefficient)*(grad_adjs(1,1,:)+grad_adjs(2,3,:)));


% 
% F = sum(F);
% G = sum(G,2);



% FoM_on_deno = T(1)+T(2)+T(3);
% FoM_on_nume = T(1);
% FoM_off_deno = T(6)+T(4)+T(5);
% FoM_off_nume = T(6);
% 
% grad_on_deno = grad_adjs(1,:)+grad_adjs(2,:)+grad_adjs(3,:);
% grad_on_nume = grad_adjs(1,:);
% grad_off_deno = grad_adjs(6,:)+grad_adjs(4,:)+grad_adjs(5,:);
% grad_off_nume = grad_adjs(6,:);
% 
% 
% G(:,1) =  ( (grad_on_nume*FoM_on_deno - FoM_on_nume*grad_on_deno)/ FoM_on_deno^2 );
% G(:,2) = ((grad_off_nume*FoM_off_deno - FoM_off_nume*grad_off_deno)/FoM_off_deno^2);
% 
% % G(:,1) =  (( (grad_on_nume*FoM_on_deno - FoM_on_nume*grad_on_deno)/ FoM_on_deno^2 ) + ...
% %              ((grad_off_nume*FoM_off_deno - FoM_off_nume*grad_off_deno)/FoM_off_deno^2))/2;
% 
% 
% % G(:,3) = grad_adjs(1,:);
% % G(:,4) = grad_adjs(6,:);
% % G  = 10*G;
% % grad = grad_adjs(1,:)+grad_adjs(5,:);
% 
% F(1) = (FoM_on_nume/FoM_on_deno);
% F(2) = (FoM_off_nume/FoM_off_deno);
% 
% % F(1) = ((FoM_on_nume/FoM_on_deno) + (FoM_off_nume/FoM_off_deno))/2;
% 
% % F(3) = T(1);
% % F(4) = T(6);
disp(T);
disp(F);

F = -F;
G = -G;

end