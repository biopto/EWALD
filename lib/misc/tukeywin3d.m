function [tuk] = tukeywin3d(siz, rs)
% % -------------------------
% siz = [4 8 6];
% rs = [0.25 0.25 0.25];
% % -------------------------

cn = num2cell(siz);
cr = num2cell(rs);
[n1, n2, n3] = cn{:};
[r1, r2, r3] = cr{:};

t1 = single(tukeywin(n1,r1));
t2 = single(tukeywin(n2,r2));
t3 = single(tukeywin(n3,r3));

%tuk13 = t1*t3';
tuk13 = bsxfun(@times, t1, t3');% Octave: 8% faster
tuk = single(zeros(siz));
for k=1:size(tuk13,2)
%    tuk(:,:,k) = tuk13(:,k)*t2';
    tuk(:,:,k) = bsxfun(@times, tuk13(:,k), t2'); % 8% faster!!! (Octave)
end

% % -------------------------
% size(tuk)
% 
% figure
% subplot(311); imagesc(squeeze(tuk(:,:,end/2+1))); axis image; colorbar
% subplot(312); imagesc(squeeze(tuk(:,end/2+1,:))); axis image; colorbar
% subplot(313); imagesc(squeeze(tuk(end/2+1,:,:))); axis image; colorbar
% % -------------------------
