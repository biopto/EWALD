fonts = 15;
show_axis_in_microns = 0;

%-----------------------------------------------------------------------
% REC = RECON;
if ~exist('n_immersion','var'); n_immersion = n_imm; end
%============================================NEW
if exist('lambda','var') && (numel(n_immersion) == size(sino_params,2))
    n_immersion = unique(n_immersion(round(sino_params(3,:), 10) == lambda));
end
%============================================
if ~exist('dxo','var') && exist('dx','var')
    dxo_ = dx; 
elseif exist('dxo','var')  
    dxo_ = dxo;
else; error('No dx or dxo!');
end 
%
if show_axis_in_microns
    Xlab = 'X (\mum)'; Ylab = 'Y (\mum)'; Zlab = 'Z (\mum)';
else
    dxo_ = 1;
    Xlab = 'X'; Ylab = 'Y'; Zlab = 'Z';
end
%
x = (1:size(RECON, 1)).*dxo_;
y = (1:size(RECON, 2)).*dxo_;
z = (1:size(RECON, 3)).*dxo_; 
%-----------------------------------------------------------------------

scale = RECON(abs(RECON(:)-n_immersion) > 0.005); 
scale = sort(scale); 
try
    scale = [scale(round(0.001*end))-0.006 scale(round(0.999*end))+0.006];
catch
    scale = [-0.015 0.015] + n_immersion;
end

COM = round(size(RECON)/2);
figure;set(gcf, 'WindowState','maximized');
subplot(231);imagesc(x,z,permute(squeeze(RECON(COM(1),:,:)),[2 1]),scale); axis image; title('y = 50%');  xlabel(Xlab); ylabel(Zlab); set(gca,'fontsize',fonts);
line([COM(3),COM(3)]*dxo_,[1,size(RECON,3)]*dxo_,'LineStyle','--','Color','w','LineWidth',1);
line([1,size(RECON,2)]*dxo_,[COM(2),COM(2)]*dxo_,'LineStyle','--','Color','w','LineWidth',1)

subplot(232);imagesc(y,z,permute(squeeze(RECON(:,COM(2),:)),[2 1]),scale); axis image; title('x = 50%'); xlabel(Ylab); ylabel(Zlab); set(gca,'fontsize',fonts);
line([COM(3),COM(3)]*dxo_,[1,size(RECON,3)]*dxo_,'LineStyle','--','Color','w','LineWidth',1);
line([1,size(RECON,1)]*dxo_,[COM(1),COM(1)]*dxo_,'LineStyle','--','Color','w','LineWidth',1)

subplot(233);imagesc(x,y,squeeze(RECON(:,:,COM(3))),scale); axis image; title('z = 50%'); xlabel(Xlab); ylabel(Ylab); set(gca,'fontsize',fonts);
line([COM(2),COM(2)]*dxo_,[1,size(RECON,2)]*dxo_,'LineStyle','--','Color','w','LineWidth',1);
line([1,size(RECON,1)]*dxo_,[COM(1),COM(1)]*dxo_,'LineStyle','--','Color','w','LineWidth',1)
COM = round(centerOfMass((abs(RECON-n_immersion)).^2));

subplot(234);imagesc(x,z,permute(squeeze(RECON(COM(1),:,:)),[2 1]),scale); axis image; title(['y = ' num2str(COM(1))]); xlabel(Xlab); ylabel(Zlab); set(gca,'fontsize',fonts);
line([COM(2),COM(2)]*dxo_,[1,size(RECON,2)]*dxo_,'LineStyle','--','Color','w','LineWidth',1);
line([1,size(RECON,3)]*dxo_,[COM(3),COM(3)]*dxo_,'LineStyle','--','Color','w','LineWidth',1)

subplot(235);imagesc(y,z,permute(squeeze(RECON(:,COM(2),:)),[2 1]),scale); axis image; title(['x = ' num2str(COM(2))]); xlabel(Ylab); ylabel(Zlab); set(gca,'fontsize',fonts);
line([COM(1),COM(1)]*dxo_,[1,size(RECON,1)]*dxo_,'LineStyle','--','Color','w','LineWidth',1);
line([1,size(RECON,3)]*dxo_,[COM(3),COM(3)]*dxo_,'LineStyle','--','Color','w','LineWidth',1)

subplot(236);imagesc(x,y,squeeze(RECON(:,:,COM(3))),scale); axis image; title(['z = ' num2str(COM(3))]); xlabel(Xlab); ylabel(Ylab); set(gca,'fontsize',fonts);
line([COM(2),COM(2)]*dxo_,[1,size(RECON,2)]*dxo_,'LineStyle','--','Color','w','LineWidth',1);line([1,size(RECON,1)]*dxo_,[COM(1),COM(1)]*dxo_,'LineStyle','--','Color','w','LineWidth',1)

pos = get(subplot(2,3,6),'Position');pos2 = get(subplot(2,3,3),'Position');
colorbar('Position', [pos(1)+pos(3)+0.005 pos(2) 0.02 pos2(2)-pos(2)+pos(4)]);
pause(0.2);