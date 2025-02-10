% Calculate a diffraction PSF
%User inputs
pixel_pitch = 5e-6;    %[m] pixel pitch
f_number    = 5;       %[fl/D] Fnumber
wavelength  = 0.555e-6 %[m] wavelength
osf         = 17;       %[#] Over sample factor
buffer      = 5;        %[pixels] Add a buffer to the PRP for motion

% calculate a sum normalized psf
[psf, x, y] = diffraction_psf(pixel_pitch/osf, f_number,wavelength);
psf = psf/sum(psf(:));

%turn it into a prp
prp = conv2(psf,ones(osf),'same')

%pad array for motion
prp_padded = padarray(prp,[1 1]*round(buffer*osf/2),0,'both')
extend = [diff(x(1:2)):diff(x(1:2)):(round(buffer*osf/2)*diff(x(1:2)))]+max(x)
x2 = [-fliplr(extend),x,extend]
y2 =  x2


%Plot the PSF
figure;
subplot(1,2,1)
imagesc(x*1e6, y*1e6, psf); % Extent for um
axis equal; %important to make it look circular
title(['Diffraction-limited PSF (Pixel Pitch: ' num2str(pixel_pitch*1e6) ' µm, f/' num2str(f_number) ')']);
xlabel('X Position (µm)');
ylabel('Y Position (µm)');
axis tight
cb=colorbar;
ylabel(cb,'Power Ratio','FontSize',10,'Rotation',270)


%Plot the PSF
subplot(1,2,2)
imagesc(x2*1e6, y2*1e6, prp); % Extent for um
axis equal; %important to make it look circular
title(['Pixel Relative Power LUT']);
xlabel('X Position (µm)');
ylabel('Y Position (µm)');
cb=colorbar;
axis tight
ylabel(cb,'Pixel Relative Power','FontSize',10,'Rotation',270)




%Example of getting the FWHM (Full Width at Half Maximum)
fwhm_val = calculate_fwhm(psf, x);
disp(['FWHM: ' num2str(fwhm_val*1e6) ' µm']);


%% Generate a PSF
function [psf, x, y] = diffraction_psf(pixel_pitch, f_number, wavelength)
% Generates the diffraction-limited point spread function (PSF).

% Input arguments:
%   pixel_pitch: Pixel pitch in meters.
%   f_number: f-number of the optical system.
%   wavelength (optional): Wavelength of light in meters. Defaults to 550nm.

% Output arguments:
%   psf: The PSF array.
%   x: The x-coordinates of the PSF grid.
%   y: The y-coordinates of the PSF grid.

if nargin < 3
    wavelength = 550e-9; % Default wavelength: green light
end

% Calculate the Airy disk radius (first minimum)
airy_radius = 1.22 * wavelength * f_number;

% Create a grid for the PSF
psf_size = ceil(6 * airy_radius / pixel_pitch); % Size of the PSF array
if mod(psf_size, 2) == 0  % Ensure psf_size is odd for centering
    psf_size = psf_size + 1;
end

x = linspace(-psf_size/2, psf_size/2, psf_size) * pixel_pitch;
y = linspace(-psf_size/2, psf_size/2, psf_size) * pixel_pitch;
[X, Y] = meshgrid(x, y);
r = sqrt(X.^2 + Y.^2);

% Calculate the PSF using the Airy disk formula
psf = zeros(size(r));
nonzero_r = r > 0;
psf(nonzero_r) = (2 * besselj(1, pi * r(nonzero_r) / airy_radius) ./ (pi * r(nonzero_r) / airy_radius)).^2;
psf(r == 0) = 1.0; % Central maximum

% Normalize the PSF (integral should be 1)
psf = psf / sum(psf(:));
end


function fwhm_val = calculate_fwhm(psf, x)
    % Calculates the FWHM of the PSF.
    [~,center_index] = find(psf == max(psf(:)));
    half_max = max(psf(:)) / 2;
    
    % Find indices where PSF is closest to half max (left and right)
    [~, left_index] = min(abs(psf(center_index,1:center_index) - half_max));
    [~, right_index] = min(abs(psf(center_index,center_index:end) - half_max));
    right_index = right_index + center_index - 1; % Adjust right index
    
    fwhm_val = x(right_index) - x(left_index);
end