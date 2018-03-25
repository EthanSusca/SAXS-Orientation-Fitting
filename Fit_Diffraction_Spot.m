% Fits a gaussian peak to a diffraction spot in an image.
%
%   Inputs - imag - array containing the diffraction pattern
%            x,y - row and column ( down, right) position of peak in pixels.
%            global X_cen, Y_cen -> row, column position of beam centre
%            r_width - defines radial extent of peak (pixels)
%            theta_width - defines angular width of peak (radians)
%                 Points considered when fitting the peak must have a
%                 radius within r_width pixels of peak and angle within
%                 theta_width of the peak.
%
%   Outputs - 
%            result.x = row position of peak.
%            result.y = column position of peak.
%            result.radius = radius of peak from beam centre (X_cen, Y_cen)
%            in pixels.
%            result.theta =  angle of peak from down direction
%                       (counterclockwise)
%            result.radius_width = radial width of peak (pixels)
%            result.theta_width =  angular width of peak (radians)
%            result.amplitude = height of gaussian.
%            result.background = constant background level
%            result.intensity = Integrated scattering power in reciprocal space
%            result.error = root mean square deviation over area where peak
%                           should fit
%            result.data =  M*N array of image where peak was fitted.
%            result.fit = M*N array of fitted peak.
%            result.weights = M*N array of zeros marking region where fit was performed.
%
%  Typical call - 
%       j=72; res(j) = Fit_Diffraction_Spot(a,x(j), y(j), 6, 0.12); figure(1); plot(res(j).y, res(j).x,'ro');figure(2);imagesc(res(j).data.*res(j).weights); res(j)
%
%   27 June, 2006 - GEST

function result = Fit_Diffraction_Spot(imag, x, y, r_width, theta_width, X_cen, Y_cen)
      
    % Grab region of image in which data is located.
    % Grab only data for which |r- r_peak| < r_width
    %                    and   |theta-theta_peak| < theta_width
    % Do this by grabbing a smaller box and using a mask.
    clearvars ymax ymin xmax xmin xtemp ytemp theta_c radii_c radius_p theta_p
    radius_p = sqrt( (x-X_cen).*(x-X_cen) + (y-Y_cen).*(y-Y_cen) ); % Peak Radius
    theta_p = atan2( (y-Y_cen) , (x-X_cen) ); % Peak Angle
    radii_c = [ radius_p - r_width, radius_p - r_width, radius_p + r_width, radius_p + r_width];
    theta_c = [ theta_p - theta_width, theta_p + theta_width, theta_p - theta_width, theta_p + theta_width];
    xtemp = radii_c .* cos(theta_c)+X_cen;
    ytemp = radii_c .* sin(theta_c)+Y_cen;          
    xmin = max(floor(min(xtemp)),1);   %ALL THIS IS THE OLD CODE
    xmax = min(ceil(max(xtemp)),size(imag,1));
    ymin = max(floor(min(ytemp)),1);
    ymax = min(ceil(max(ytemp)),size(imag,2));
    xv = ((xmin:xmax)' - X_cen)*ones(1, ymax-ymin+1);%
    yv = ones(xmax-xmin+1,1)*((ymin:ymax)-Y_cen);
    
    % The allowed data region! 
    r = sqrt( xv.*xv + yv.*yv); % Radii at each point in allowed region
     %disp(r); disp(theta_p);%{disp(theta);
    theta = atan2(yv, xv); theta = mod( theta-theta_p+pi, 2*pi) + theta_p-pi; % Angle at each allowed point.
    data = imag(xmin:xmax,ymin:ymax); % Values of imag within allowed region.    
    weights = (r>= (radius_p-r_width)).*(r <= (radius_p+r_width)).*( theta >= (theta_p-theta_width)) .*(theta <= (theta_p+theta_width));
    
    % Estimate peak shape using Centre-of-Mass
    param = estimate_param(data.*weights,r,theta );
    
    % Now fit to Gaussian.
    param_lower = [ param(1)-r_width/2, param(2) - theta_width/2, param(3)/2, param(4)/2];
    param_upper = [ param(1)+ r_width/2, param(2) + theta_width/2, param(3)*2, param(4)*2];
    param = robust_fit(data, weights,r,theta,param, param_lower, param_upper);
    model = gaussian2d(r,theta,param);
    [err_fit, param_lin] = fit_linear_terms(data, model, weights); 
    
    % Save Results
    result.x = X_cen + param(1)*cos(param(2));
    result.y = Y_cen + param(1)*sin(param(2));
    result.radius = param(1);
    result.theta = param(2);
    result.radius_width = param(3);
    result.theta_width = param(4);
    result.amplitude = param_lin(2);
    result.background = param_lin(1);
    result.intensity = param_lin(2) * (2*pi)^(3/2) * (result.theta_width)^2 * (result.radius_width^2 + result.radius^2 ) * result.radius_width;
    result.error = err_fit;
    result.data = data;
    result.fit = param_lin(1) + param_lin(2)*model;
    result.weights = weights;
end
    
function result = robust_fit(data, weights,r,theta,param_init, param_lower, param_upper)
    % Robust fit of gaussian peak to data.
    %
    % Inputs - 
    %           data -> M*N array of values describing peak.
    %           weights -> M*N array of experimental weights for each
    %               measurement.
    %           r -> M*N array with radius of each point.
    %           theta -> M*N array with angle at each point.
    %           param_init -> Initial gaussian parameters
    %                       param_init(1) - radius of peak centre
    %                       param_init(2) - angle of peak centre
    %                       param_init(3) - radial width of peak.
    %                       param_init(4) - angular width of peak.
    %           param_lower -> Lower bound on gaussian peak parameters
    %           param_upper -> Upper bound on gaussian peak parameters
    %
    % Outputs - 
    %           results -> Fitted gaussian parameters
    %                       results(1) - radius of peak centre
    %                       results(2) - angle of peak centre
    %                       results(3) - radial width of peak.
    %                       results(4) - angular width of peak.

    % Initialize search
    result = param_init;
    step_size = (param_upper - param_lower)/10;
    Nparam = length(result);
    steps = 0;
    tol = 1e-6;
    Nmax = 30;

    % Search until Nmax search steps have been unsuccessful.    
    while(steps<Nmax)
    
    % Get current error.
    model = gaussian2d(r,theta,result);
    err_best = fit_linear_terms(data, model, weights); 
    
    flag = 0;
    
    % Try tweaking each parameter in turn to improve fit
    for j=1:Nparam
       
        % Tweak parameter
        test(j) = result(j) + step_size(j); test = min(test, param_upper); test = max(test, param_lower);
        model = gaussian2d(r,theta,test);
        err_test = fit_linear_terms(data,model,weights);
        
        % Evaluate tweak
        if (err_test<err_best) % Got better --> keep going that way
            err_best=err_test;
            result(j)= test(j);
            step_size(j)= step_size(j)*1.5;
        else
            flag = flag+1; % Got worse --> Note and try going other way
            step_size(j) = step_size(j)*-0.5;
        end
        
    end
    
    % Evaluate success/failure of last round of tweaking.
    % If all tweaks failed or all parameters are within tolerance, the step
    % failed.
    if (flag==Nparam)||(max( abs(step_size)./(param_upper-param_lower))<tol)
        steps=steps+1;
    else
            steps=0;        
    end
    
    % Optional display of progress
    %fprintf('Steps = %d, Best Error = %f, Parameters = %f %f %f %f \n', round(steps), err_best, result(1), result(2), result(3), result(4));
    end
end

function result = estimate_param(data, r,theta)
    % Estimate mean and standard deviation of peak using Centre-of-Mass
    %
    % Inputs - 
    %           data -  M*N array of intensities of peak.
    %           r -     M*N array of radii at each point.
    %           theta - M*N array of angles at each point.
    %
    % Outputs -
    %           result(1) - mean radius of peak
    %           result(2) - mean angle of peak 
    %           result(3) - standard deviation of peak radius
    %           result(4) - standard deviation of peak angle.
    
    data = data / sum(sum(data));
    rmean = sum(sum( data.*r )) ;
    tmean = sum(sum( data.*theta )) ;
    rsq = sum(sum( data.*(r-rmean).*(r-rmean) )) ;
    tsq = sum(sum( data.*(theta-tmean).*(theta-tmean) ));

    result = [ rmean, tmean, sqrt(rsq), sqrt(tsq)];
end
    
function [err, param] = fit_linear_terms( data, model, weights)
   % Fit the following equation by least squares
   %       data = param(1) + param(2)*model
   %
   % Inputs
   %        data --> M*N matrix of values
   %        model --> M*N matrix of values
   %        weights -> M*N matrix of weighting coefficients
   %
   % Outputs - 
   %        err -> standard deviation.
   %        param(1) -> constant background
   %        param(2) -> scaling factor for model.

    y  = [ sum(sum(data.*weights)); sum(sum(data.*model.*weights))];
    M = [ sum(sum(weights)), sum(sum(model.*weights));
          sum(sum(weights.*model)), sum(sum(model.*model.*weights))];
    
    param = M\y;
    tmp = data - (model * param(2) + param(1));
    err = sum(sum( tmp.*tmp.*weights))/sum(sum(weights));
    err = sqrt(err);
end
    
function result = gaussian2d(r,theta,p)
    % Compute a gaussian peak
    %   r --> array of radii
    %   theta --> array of values of theta
    %   p --> parameters describing peak.
    %       mean radius = p(1)
    %       mean angle = p(2)
    %       std dev radius = p(3)
    %       std dev angle = p(4)
    %
    %   result --> value of gaussian at all points in r and theta.
    
    result = exp(-0.5*(r-p(1)).^2./(p(3)^2)-0.5*(theta-p(2)).^2./(p(4)^2))  ;

end
    