% Generate OCT test data that are to be combined with ODT Sinogram
% and optionally save it to MAT file 
%
% CompAmp:  select complex amplitude data type: 'ones_array' - ones-padded array, 'UMK' - OCT data from .complex file
% SINOdim:  dimension of the SINOamp or SINOph
% k_slice:  select 2D complex amplitude (xy) for 'k_slice' from the range of k1..kN (kN = Z; see the code below)      
function [octSINOamp, octSINOph, oct_lambda, octNA, octRayXY] = GenOctData(CompAmp, SINOdim, SINOdx, k_slice)

if nargin<4 && strcmp(CompAmp, 'UMK')
    error('Specify the fourth argument of function GenOctData: k_slice');
end

save_mat_file = 0;  % save OCT data to MAT file
show_plots = 0;  % 1: Complex amplitude plots, 2: add Fp plot
rescale = true;  % rescale octSINOamp and octSINOph to fit SINOamp and SINOph size 
% A path to OCT meaurement data file (.complex):
OCTfilePath = 'D:\PW\Programowanie\Dane OCT Torun\2018-09-11\ComplexData_Bs500xAs500xPx70.complex';

%===============================================================
% OCT system parameters:
octNA = 1.3;  % NA of the OCT system
% A wavelength that corresponds to k (indicated by k_slice) for which the OCT complex amplitude 
% was selected from 3D array of spectra:
oct_lambda = 0.6328;  % OCT complex amplitude wavelength (according to k); 
                                  % should be equal or very close to ODT light source wavelength (=lambda)
octM = 48;  % OCT system lateral magnification = ODT system magnification
oct_downsamp = 3.0808;  % not known yet whether downsampling is to be performed for OCT data
oct_cam_pix = 3.45;  % camera pixel size in microns - the value not known yet!!!
%
% Camera pixel size in object space (takes into account downsampling): 
%oct_dx = oct_cam_pix * oct_downsamp / octM;  % the value not known yet!!!
%=====================================
oct_dx = SINOdx;  % added temporarly for tests
%=====================================
octRayXY = [0.0; 0.0];  % projection illumination direction vector (on-axis)
%octPeak =  [0.0; 0.0];  % not impartant!
%===============================================================
% 
if strcmp(CompAmp, 'UMK')  % use OCT data provided by UMK
   %k_slice = 40;  % get x-y plane of 3D spectral OCT data for k=k_slice

    % Read processed measurement data from Optical Coherence Microscopy.
    % Imported array consists of 1D data used further to reconstruct a 3D array of the OCT spectral data (processed spectral fringes)
    %
    % Open .complex file:
    oct_file = fopen(OCTfilePath,'r');
    header = fread(oct_file, 3, 'int');
    X = header(1); Y = header(2); Z = header(3);  % X,Y,Z indicate OCT spectral data dimensions for 3D array
    dim = X*Y*Z; % dimension of 1D array
    % Raw data from .complex file stored in 1D array of double
    data = fread(oct_file, dim*2, 'float32=>single'); 
    fclose(oct_file);
    % Combine real and imag parts from 'data' array to form 1D spectral data:
    data_comp = complex(zeros(1,dim)); % 1D array of the OCT compex signal data (spectra after IFFT)
    for ii = 1:2:dim*2
        data_comp((ii+1)/2) = data(ii) + 1i*(data(ii+1)); 
    end
    clear data;
    % Retrieve original 3D array dimensions of the spectral OCT data
    data_comp = permute(reshape(data_comp, [Z X Y]), [2 3 1]);  % swap 1st dimension with the 3rd one
    %data_comp = flip(data_comp, 2);
   
    octSINOamp = single(abs(squeeze(data_comp(:,:,k_slice))));
    octSINOph = single(angle(squeeze(data_comp(:,:,k_slice))));
    clear data_comp;
    
    %=========================================================
    % Resize octSINOamp and octSINOph to fit SINOamp and SINOph
    if rescale && X ~= SINOdim && oct_dx ~= SINOdx
    %if rescale
        % assuming square projection
        octSINOamp = match_resolution(octSINOamp,oct_dx, SINOdx, SINOdim);
        octSINOph = match_resolution(octSINOph, oct_dx, SINOdx, SINOdim);
    end
    %=========================================================
    
    if show_plots>1
        Up = single(log(octSINOamp)+1i*octSINOph);  % Rytov
        Fp = fftshift(fft2(ifftshift(Up))) * oct_dx^2;
    end
    
    if show_plots
        % x-y plane for k number k_slice across the OCT 3D spectral data
        figure(1001); %set(figure(1001),'Units', 'Normalized', 'OuterPosition', [0 0 1 1]); 
        sgtitle(['XY plane through 3D spectral data for k #',num2str(k_slice), ' of ', num2str(Z)],  'FontSize', 14)
        %colormap gray
        %
        subplot(1,2,1)
        imagesc(octSINOamp)
        set(gca,'YDir','normal')
        colorbar
        axis image
        %ylabel('\fontname{Times New Roman}y', 'FontSize', 16)
        ylabel('y', 'FontSize', 14)
        xlabel('x', 'FontSize', 14)
        title('Amplitude', 'FontSize', 12)
        %
        subplot(1,2,2)
        imagesc(octSINOph)
        set(gca,'YDir','normal')
        colorbar
        axis image
        ylabel('y', 'FontSize', 14)
        xlabel('x', 'FontSize', 14)
        title('Phase', 'FontSize', 12)
        
        if show_plots>1
            figure(1002);
            imagesc(log10(0.01+abs(Fp.'))); axis image;
            xlabel('Kx'); ylabel('Ky'); title('Fp');
            colormap jet
        end
    end
%
elseif strcmp(CompAmp, 'ones_array') && nargin<4  % use zero-padded array
    octSINOamp = single(ones(SINOdim));
    octSINOph = single(ones(SINOdim));
else
    error('Use ''UMK'' or ''ones_array'' option for ''CompAmp'' parameter of GenOctData function.');
end
if save_mat_file
    save(['OctTestData_',CompAmp,'.mat'], 'octSINOamp', 'octSINOph', 'oct_lambda', 'octNA', 'octRayXY');
end
    

