function LMS_Feat = SRGB2LMS(RGB_Feat)

XYZ_RGB = [0.4124, .3576, .1805; 0.2126, 0.7152, 0.0722; 0.0193, 0.1192, 0.9505];
XYZ_RGB = XYZ_RGB ./ repmat(sum(XYZ_RGB, 2), [1, 3]);
LMS_XYZ = [0.7328, 0.4296, -0.1624; -0.7036, 1.6975, 0.0061; 0.0030, 0.0136, 0.9834];
LMS_RGB = LMS_XYZ * XYZ_RGB; %CIECAM02


%see if we need to reshape the values
if size(RGB_Feat, 3) == 3 || size(RGB_Feat, 2) ~= 3
    is_image = true;
    [sy, sx, sz] = size(RGB_Feat);
    R = RGB_Feat(:, :, 1);
    G = RGB_Feat(:, :, 2);
    B = RGB_Feat(:, :, 3);
    RGB_Feat = [R(:), G(:), B(:)];
else
    is_image = false;
end
RGB_Feat = single(RGB_Feat);

% normalize:
if any(RGB_Feat(:) > 1)
    RGB_Feat = RGB_Feat ./ 255;
end

% remove sRGB gamma nonlinearity:
mask = RGB_Feat <= 0.04045;
RGB_Feat(mask) = RGB_Feat(mask) ./ 12.92;
RGB_Feat(~mask) = ((RGB_Feat(~mask) + 0.055) ./ 1.055) .^ 2.4;

%now do the conversion
LMS_Feat = (LMS_RGB * RGB_Feat')';

LMS_Feat = max(LMS_Feat, 0);
LMS_Feat = min(LMS_Feat, 1);

%turn back into an image if necessary
if is_image
    L = reshape(LMS_Feat(:, 1), [sy, sx]);
    M = reshape(LMS_Feat(:, 2), [sy, sx]);
    S = reshape(LMS_Feat(:, 3), [sy, sx]);
    LMS_Feat = zeros(sy, sx, 3, 'single');
    LMS_Feat(:, :, 1) = L;
    LMS_Feat(:, :, 2) = M;
    LMS_Feat(:, :, 3) = S;
end
