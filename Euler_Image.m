%
%           Euler_Image
%
%   Description: This function returns the orientation data in a format
%   that can be read in by a texture analysis software.
%
%   IN:  imageAddress   This a string 'exampleImage.tif' that points to the
%   image (Tiff files here) to be read in for orientation analysis.
%
%   OUT: fit_res = [phi1, PHI, phi2, ncount, fract]  A vector containing the three Euler angles (in
%   degrees) using Bunge notation, passiver movements (lab frame moves)
%   [phi1, PHI, phi2] are the Euler angles
%   ncount  is the number of identified Bragg peaks fit by single
%           orientation diffraction pattern
%   fract   is the fraction of identified Bragg peaks fit by single
%           orientation diffraction pattern

function [ph1, PHC, ph2, nfits, fract, za_h, za_k, za_l] = Euler_Image(imageAddress, X_cen, Y_cen, pxlsperq)

img = imread(imageAddress);  % Read in image. img is a Matrix of intensities
    img = double(img);
       img(585,513)=0;   % Reset of a bad pixel on the detector 
       img(887,598)=0;
display_img_fits = 0;    % Can see the fits in action if display_img_fits = 1
% X_cen = 431.0;      % The column where the beam center is located
% Y_cen = 723.8;      % The row where the beam center is located
% pxlsperq = 34.0;    % How many pixels correspond to q100   ???
                    % In the case of double gyroid, primary peak is q211,
                    % so then pxlsperq = q100 = q211/(8)^0.5
                    % (8=2^1+1^2+1^2)
% Bragg Peak Characteristics to help fit a 2D Gaussian
r_width = 4;    % A guess at how many pixels corresponds to the std dev about q
theta_width = 3/180 * pi; % Radians corresponding to standard dev along azimuth

    %  Here we blot out the high intensity pixels from x-rays that leak out
    %  around the beam stop. All points below q200 will be made zero. This is
    %  because the primary scattering peak is q211
    BetterBeamStop = ones(size(img));  % It's effectively putting a black at the beam center         
    for h = 1:size(img,2)
        for i = 1:size(img,1)
            rad = sqrt( (i-Y_cen)^2+(h-X_cen)^2 );
            if rad < 2*pxlsperq   %2*q100 = q200, which is the radius of the spot.
                BetterBeamStop(i,h) = 0; %Between q=root(8) and q=root(10) is q=3
            end
        end
    end
    img = img.*BetterBeamStop+1;   % Now apply the better beam stop to image

% Now it's time to initialize variables that will be found
 peaks.r = []; peaks.theta = []; peaks.intensity = []; % list of peak attributes
 peaks.x = []; peaks.y = [];
 
 img_log = real(log10(img));
[xm, ym] = find(img_log == max(img_log(:)));  % Find x,y coords of max intensity (Bragg peak)
xm = xm(1); ym = ym(1);  % If there are two identical intensity piixels, pick one
maxint = img(xm,ym);   % The entensity of that max intensity pixe

thresh_loq = maxint/5;  % or maint/10 ...threshold for diffraction peak fitting at low q (211 and 220)
thresh_hiq = 30;%25;  %for 96k and 71k polymers...The high q threshold should be the background intensity maximum outside low q
thresh = thresh_loq;    %  We start witht the larger, loq intensity threshold

%  LOOP FINDING LOCATION OF DIFFRACTION SPOTS STARTS HERE
loopn = 0 ; loopthresh = 100;  %Count number of diffraction spots fit, if >100: quit
     while maxint > thresh
        img_log = real(log10(img));
        [xm, ym] = find(img_log == max(img_log(:)));  % Find x,y coords of max intensity (Bragg peak)
        xm = xm(1); ym = ym(1);  % If there are two identical intensity piixels, pick one
        maxint = img(xm,ym);   % The intensity of that max intensity pixel
        % FIT DIFFRACTION SPOT finds the parameters for the peak at x,y
        result = Fit_Diffraction_Spot(img, xm,ym,r_width, theta_width, X_cen, Y_cen);
        peaks.x = [peaks.x, result.x]; peaks.y = [peaks.y, result.y];
        peaks.r = [peaks.r, result.radius]; peaks.theta = [peaks.theta, result.theta];
        peaks.intensity = [peaks.intensity, result.intensity];
        % Replace the high intensity peak with a mask of zeroes where high
        % intensity pixels previously existed
        fitmask = result.data; % This data replaces the Bragg peak region in the image
        datacut = 0.02;  % All data greater than 2.5% of the main peak is zeroed
        [dimx,  dimy]= size(result.fit);
        for i = 1:dimx
            for j = 1:dimy
                if result.data(i,j) > datacut*result.amplitude
                    fitmask(i,j) = 1;
                end
            end
        end
        img(round(result.x-dimx/2):round(result.x+dimx/2)-1, round(result.y-dimy/2):round(result.y+dimy/2)-1) = fitmask;
        % plot(fitmask) 
               
        % First test to see if the brightest pixel is not bright enough to
        % be low q scattering (the high intensity 211 220 peaks)
        if maxint < thresh   % If the intensity of the main peak is less than the threshold requirement
            BlackSpotMask = ones(size(img));   % After mask (binary) used, will set threshold to lower value
            for h = 1:size(img,2)
                for i = 1:size(img,1)
                    rad = sqrt( (i-Y_cen)^2+(h-X_cen)^2 );
                    if rad < 3.2*pxlsperq % equal to 3x a* (unit inverse vectors: (a*, b*, c*))
                        BlackSpotMask(i,h) = 0; %Between q=root(8) and q=root(10) is q=3
                    end
                end
            end
            img = img.*BlackSpotMask+1;  % take away everything near center
            thresh = thresh_hiq;   % New threshold since we are past the low q threshold
        end    
       
        % This is to prevent fitting more than 100 Bragg peaks
        loopn = loopn + 1;
        if loopn > loopthresh
            thresh = maxint;  % Will cause while loop to end
        end   
     
        if display_img_fits == 1    %  FIGURE 1: 1 x 2 subplot with img, data, and fit
            figure(1); clf; % Open figure
               img_log = real(log10(img));
               wid = 38^0.5*pxlsperq;
               imagesc(img_log(round((Y_cen-wid):(Y_cen+wid)),...
                               round((X_cen-wid):(X_cen+wid))), [0, 4.5]); 
                             % Display Image scaling intensity  from 1 to 10^4.5
               axis image; colormap('jet');  colorbar ; 
        end
     
     end  % END WHILE MAXINT > THRESH LOOP  % % % % % ALL PEAKS IDENTIFIED

% Now scale the x,y peak locations relative to q100
datafit = zeros(length(peaks.x),2);
for j = 1: length(peaks.x)   
    datafit(j,2)=(peaks.y(j)-X_cen)/pxlsperq; %  q_x list relative to q100
    datafit(j,1)=(peaks.x(j)-Y_cen)/pxlsperq; %  q_y
end
 
[ph1, PHC, ph2, nfits, za_h, za_k, za_l] = Fit_Orientation_Bunge(datafit); %Now adds total number of fits found (loopn) to compare
fract = nfits/length(peaks.x);



end



       