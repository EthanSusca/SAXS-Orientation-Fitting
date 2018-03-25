%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       WriteEBSDtxt
%           The output of this script is a text file formatted in such a
%           such that it can be easily read-in by electron backscatter
%           diffraction analysis software. Each row contains orientation
%           information for a single pixel. Pixel location
%
%       OUT: m-by-12 text file where each row contains:
%               [ index (row number), phase (arbitrary), 
%                   x location (zero indexed), y location (zero indexed),
%                       phi1,PHI,phi2, n_fits, fract, - za_h, za_k, za_l]
%       index = the data point, or the row number
%       phase = 1 for gyroid, this column is for post-processing, so that
%           the user can define constraints on what combination of n_fits
%           and fract constitute a succesfull fit to calculated data.
%       x = integer from 0 to number of columns-1
%       y = integer from 0 to number of rows-1
%       phi1,PHI,phi2 = the three Euler angled defined using Bunge notation 
%           with passive rotations
%       n_fits = number of identified Bragg peaks (from experimental data)
%           that are accounted for by a single orientation
%       fract = n_fits/(total number of identified Bragg peaks 
%

rows = 1:121;  % Number of rows scanned 
cols = 1:181;  % Number of columns scanned

X_cen = 431.0;      % The column where the beam center is located
Y_cen = 723.8;      % The row where the beam center is located
pxlsperq = 34.0;    % How many pixels correspond to q100   ???

problempixls = [0 0 0];  % To record the failed index, row, and column
ebsd443n1 = zeros((length(rows))*(length(cols)),12);  % initiate text file
for n = 1: length(rows)
    for m = 1:length(cols)
        count = (n-1)*length(cols) + m;  %
        
            if cols(m) < 10     % import image address, format depends on image column number
                imageAddress = sprintf('EMS443-n1/scan_%d_master0000%d.tif',rows(n),cols(m));
            elseif cols(m) < 100 && cols(m) > 9
                imageAddress = sprintf('EMS443-n1/scan_%d_master000%d.tif',rows(n),cols(m));
            elseif cols(m) > 99
                imageAddress = sprintf('EMS443-n1/scan_%d_master00%d.tif',rows(n),cols(m));
            end
            
            % Get the Euler angles, number of successful fits, and fraction
            % fit (out of all Bragg peaks identified) with Euler_Image.
            try
            [phi1, PHI, phi2, nfits, fract, za_h, za_k, za_l] = Euler_Image(imageAddress, X_cen, Y_cen, pxlsperq) ;
            catch
            phi1=0; PHI=0; phi2=0; nfits=0; fract=0; za_h=0; za_k=0; za_l=0;
            problempixls = [problempixls; count rows(n) cols(m)]; 
            dlmwrite('prob_pxls.txt', problempixls, '\t');
            end
        ebsd443n1(count,1:12) = [count, 1, m-1,n-1, phi1, PHI, phi2, nfits, fract, za_h, za_k, za_l];
        dlmwrite('ebsd443n1_piD.txt', ebsd443n1, '\t');   % record every loop for data safety
        disp([rows(n) cols(m) ebsd443n1(count,8:12)]);% This is to keep track of the progress processing 
    end
end

%convert the m*10 matrix into txt file
dlmwrite('ebsd443n1_piD.txt', ebsd443n1, '\t');