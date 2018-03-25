%
%               Fit the Peak Locations Scaled (to q100) to return  
%                   the 3 Euler Angles, Bunge Notation, Passive Rotation
%   
%   IN:   an m x 2  matrix consisting of two colums: x and y locations
%           of Bragg peaks
%
%   OUT: The best fit phi1, PHI, phi2 Euler Angles using Bunge Notation
%
%   SETTINGS:
%       N_order: an integer corresponding to the largest h,k, or l
%               ie N_order = 5 means that possible diffraction spots from 
%                (0 0 1) through (5 5 5) will be explored as possible
%                reflections in recipricol space
%       ctoff:  number of radians corresponding to the maximum angle from
%                the Ewald plane for a reflection to be located in order
%                for it to still count as a diffraction spot
%
%   Dependencies
%       rmatrixbunge.m
%       spots_new.m
%       error_hkl.m
%

function [ph1, PHC, ph2, nfits, za_h, za_k, za_l ] = Fit_Orientation_Bunge(data)

% res = [phi1, PHI, ph2, n_fits]


% SETTINGS
display_plot = 0; % 0 means no fitting plot, 1 means it will display plot
N_order=5;        % scan h,k,l each from 0:5
ctoff=8*pi/180; %8*pi/180; % generate spot (on plot) if within 8 degrees of the Ewald plane
                % This should be about the full width at half max of a
                % Bragg peak in azimuthal directions (constant |q|)
criteria2=5;    % sum of the squares of the error in angle and radius 
                % between a found and calculated diffraction spot
                % Necessary condition for a Bragg reflex to be consistent 
                % with a calculated diffraction spot is that error < crit2
                % crit2 is relative to q = 0.05 (see error_hkl.m) where
                % error = (error_angle/0.05)^2+(error_radius/0.05)^2 
                % so criteria2=5 means max angular error = 6 degrees and
                % max radial error is q = 0.11 (this means 1/10th of q100

%% COARSE SCAN start by initially scanning angles in COARSE  6 degree increments to 90 degrees
%     (high symmetry system means this samples redundant orientations. 
k=1;  % start indexing for preresult at 1...  a loop counter to sort through to find best fit at end of search
for phi1=1:5:86         % 1,6,11,16,26,...81,86
    for PHI=2:5:87         % 18^3 = 2197 diffferent COARSE orientations
        for phi2=3:5:88      
                % rmatrixbunge returns the 3x3 rotation matrix for
                % phi1,PHI,phi2
            resmtrix=rmatrixbunge(phi1,PHI,phi2); 
                % Then generate h,k,l values that are consistent with the
                % rotation matrix and symmetry (see spots_new.m)
            res = spots_new(resmtrix.rmtrixbunge, N_order,ctoff); %result
                    % is a list of spot locations, x,y to compare with
                    % the experimentally fit Bragg peak x, y
            ncount=0;
            minm=0;    % initialize metrics for how good a fit it is
            erro=0;
            for i=1:length(data)
                for j=1:length(res)
                erro=error_hkl(res(j,4:6),data(i,:)); % criteria 2 comparator
                if(erro<criteria2)                      
                    ncount=ncount+1;  % congrats it's a hit!
                    minm=minm+erro;   % record the distance between experiment and calcuation
                end
                end
            end

            preresult(k).phi1= phi1; 
            preresult(k).PHI= PHI;  
            preresult(k).phi2= phi2;
            preresult(k).ncount=ncount;
            preresult(k).min=minm;
            k=k+1;
        end
    end
end

%find the orientation with the most hits and then
max=preresult(1).ncount;
idx=1;
for i=2:length(preresult)
    if(preresult(i).ncount>max)
    max=preresult(i).ncount;
    idx=i;
    end
end
% end
% idx1=find([preresults.ncount]>preresults(idx).ncount-1); 
% [c,idx]=min([preresults(idx1).min]);

%%  FINE SCAN Refine peak locations now
% get down to nearest 0.5 degree in all rotations
k=1;
for phi_1 = (preresult(idx).phi1-3.0):0.5:(preresult(idx).phi1+3.0)    
    for PHI_C = (preresult(idx).PHI-3.0):0.5:(preresult(idx).PHI+3.0)   
        for phi_2 = (preresult(idx).phi2-3.0):0.5:(preresult(idx).phi2+3.0)
            resmtrix=rmatrixbunge(phi_1,PHI_C,phi_2);
            res = spots_new(resmtrix.rmtrixbunge, N_order,ctoff);
            ncount=0;
            minm=0;
            erro=0;
            for i=1:length(data)
                for j=1:length(res)
                        erro=error_hkl(res(j,4:6),data(i,:)); 
                    if(erro<criteria2)                   
                        ncount=ncount+1;
                        minm=minm+erro;   
                    end
                end
            end
            results(k).phi_1=phi_1;   
            results(k).PHI_C=PHI_C;
            results(k).phi_2=phi_2;
            results(k).ncount=ncount;
            results(k).min=minm;
            k=k+1;
        end
    end
end
%%    END second set of for loop for peak refinement to 0.5 degrees

% Below: find the largest match number, if two have the same number of "hits", find
% the fit with the minimum residue
max=results(1).ncount;
idx=1;
for i=2:length(results)
    if(results(i).ncount>max)
    max=results(i).ncount;
    idx=i;
    end
end
idx1=find([results.ncount]>results(idx).ncount-1);
[c,idxmin]=min([results(idx1).min]);

%%  Generate Rotation Matrix if you would like to find the zone axis
    aa=rmatrixbunge(results(idx1(idxmin)).phi_1,results(idx1(idxmin)).PHI_C,results(idx1(idxmin)).phi_2);
    if display_plot == 1
        figure(2);clf;
        show_dots_new(data,aa.rmtrixbunge,ctoff); axis image;
        pause
    end
% THE OUTPUTS BELOW
ph1 = results(idx1(idxmin)).phi_1;
PHC = results(idx1(idxmin)).PHI_C;
ph2 = results(idx1(idxmin)).phi_2;
nfits = results(idx1(idxmin)).ncount;
        zoneaxis = aa.rmtrixbunge\[0 0 1]'; % left multiply the inverse(rotation matrix)
        za = sort(abs(zoneaxis));   % in ascending order za(1), za(2), za(3)
        za_h = za(3);
        za_k = za(2);
        za_l = za(1);
% below will return nearest h, k,l integer relative to max integer/index
%{
        max_ind = 20;  % so 2 1 1 = 10 5 5  or 1 0 0  = 10 0 0

za_h = round( max_ind*za(3) );  % my convention is largest index first
za_k = round( max_ind*za(2) );
za_l = round( max_ind*za(1) );
%}

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       dependent sub-routines: rmatrixbunge, spots_new, error_hkl
%
function res=rmatrixbunge(phi1,PHI,phi2)
        %generate the rotation matrix (res.rmtrixbung) based on Bunge convention
        %Bunge are the angles about Z, X', Z' in that order
res.rmtrix0=[1, 0, 0; 0, cosd(PHI), sind(PHI);0,-sind(PHI), cosd(PHI)]*[cosd(phi1),sind(phi1),0;-sind(phi1),cosd(phi1),0;0,0,1];
z2axis=[cosd(phi2),sind(phi2),0;-sind(phi2),cosd(phi2),0;0,0,1];
res.rmtrixbunge=z2axis*res.rmtrix0;
end




function res = spots_new(a, N_order,ctoff)
N=N_order;
%algyroid=[2 6 8 10 12 14 16 18 20 22 24 26 30 32 34];
gyroid=[6 8 14 16 20 22 24 26 30 32]; %38 40 42 %46 48]; % note compression (330) = 18 peak appears
%plumber=[2 4 6 8 10 12 14 16 18 20 22 24 26 30 32 34];
%ddiamond=[2 3 4 6 8 9 10 11 12 14 16 17 18 19 20 21 22 24 26 27 29 30 32 33 34];
%index=[1 2 3 6 8 10 12 14 16 18 20 22 24 26 30 32 34];
%algyroid=[2 6 8 10 12 14 16 18 20 22 24 26];

% Initialize the data set.
Nmax = (2*N+1)^3; % Maximum possible number of reflections.
res = zeros(Nmax,7);
Nref = 0;

phase_array=gyroid;
for h=-1.0*N:N
    for kay=-1.0*N:N
        for l=-1.0*N:N
             flag=0;
             square=h*h+kay*kay+l*l;
             for i=1:length(phase_array)
                if (square==phase_array(i)) flag=1; end
             end
             if (flag==1)
             %if (mod(round(h+k+l),2)==0)
              test = a*[ h kay l]';
              angl = norm(test(3))/ (norm(test) + (norm(test)==0)); % norm() calculates abs() for a scalor and distance sqrt(x^2+y^2+z^2) for a vector 
              angl = acos(angl)-pi/2;
              angl=norm(angl);
              if (angl<ctoff)
              Nref = Nref+1;
              res(Nref,:) = [h, kay, l , test', angl ]; 
              end
            end  
        end
    end
end
% Resize array to just include valide points.
res = res([1:Nref],:);
% Sort by angle.
[tmp , idx] = sort(res(:,7));
res = res(idx,:);
end

function  result = error_hkl(spot, q)
    std_angle = 0.05; std_q = 0.05;
        rs = norm(spot);
        rq = norm(q);
        if ((rs*rq)>0)  % This should always be the case unless norm() = 000
            angle_err = acos( (spot(1)*q(1)+spot(2)*q(2))/(rs*rq));
        else
            angle_err=0;
        end
        radius_err = rs-rq;
        
        result = (angle_err/std_angle).^2 + (radius_err/std_q).^2;
end    

% function show_dots(data,mtrix,ctoff)
%
%  Plots to figure the diffraction dots in data.
%  Computes (using spots()) the lattice of allowed reflections from a crytal with unit axes mtrix(:,j).
%  The allowed mosaicity is specified by ctoff in radians.  Typical value is 6 degrees or ctoff=0.1.
%
%  Has three constants you may need to adjust
%               N_order - Largest lattice sites considered [N_order, N_order, N_order ]
%                         This is an integer.  eg. if 4 only peaks out to q^2=16 are covered properly.
%               Range -   Range of the axes in q_x/y 
%                         6 or 10 is a good value for our data set.  
%               cmap - Determines the colours of each diffraction order.
%
%  Dependencies - spots()

function show_dots_new(data, mtrix, ctoff)

    N_order=5; %since the largest h^2+k^2+l^2 is 22.
    Range=7;%  ????
    
    %plot(data(:,1), data(:,2),'b+'); % Plot the measured diffraction spots on as crosses.
    plot(data(:,2), data(:,1),'b+'); %plot in image orientation
%    axis image; 
    axis([ -1 1 -1 1]*Range); % Set the axes.
    hold on;
    
    % Choose colour scheme of diffraction order.
    % Colour coding is handy if you want to match a [2 2 2] peak to the green rings.
    %  Black is the simpler to see, though.
    %cmap={'bo','go','co','mo','yo','ko'};
    cmap={'ko','ko','ko'};
    
    set(gca,'Ydir','reverse');
    
    u = spots_new(mtrix, N_order,ctoff);  % Compute the lattice for a crystal with unit axes mtrix(:,j)
    xlabel('q_x'); ylabel('q_y');
    
    % Plot all dots in the lattice that are with ctoff radians of the x-y plane (Small Angle Approx to the Ewald Sphere)
    for j=1:length(u)
        
             idx=(u(j,1)*u(j,1)+u(j,2)*u(j,2)+u(j,3)*u(j,3))/2; % index for the color code
             idx=floor(mod(idx,length(cmap)))+1;
                  %plot(u(j,4), u(j,5),cmap{idx}); 
                  plot(u(j,5), u(j,4),cmap{idx});  %plot image orientation      
    end
end
