function [pref_ori, CV] = VC_ori_tuning(theta, r, modes)
    
% 	Returns, for a dataset of responses (spikes/s) for given
% 	orientations, values for the mean vector.
% 
% 	r:	responses
% 	theta:	corresponding orientations
% 	modes:	number of modes; ie, whether uni, bi, or, tri or quadri-
% 		modal.  The default is bimodal (ie, not direction 
% 		selective).  See note below
% 	bin:	Boolean.  Whether bin compensation is requested.  This
% 		only affects the mean vector length(l) and derivative
% 		statistics.  Where the number of bins > 12, compensation
% 		is generally not necessary (Batschelet)
% 	bwdth:	radial length of arc of the bins.  necessary only if
% 		bin=True
% 	ci:	percentage for which confidence intervals will be
% 		calculated
% 
% 	Returns:	
% 
% 		pref_ori:	mean angle or angles of mean axis
% 		CV:         (0-1) - 0 means all response is in one orientation, 1
%                   means responses are equally distributed accross
%                   orientations.
% 
% 	Notes:
% 
% 	-- Modes --		
% 	Modes argument identifies how many modes the dataset is believed
% 	to possess.  This is important as the majority of circular
% 	statistics presume a unimodal dataset.  To compensate, the 
% 	angles of the dataset are multiplied by the amount of modes 
% 	to generate a unimodal dataset.  This compensation PRESUMES THAT
% 	THE DIFFERENCE IN ANGLE BETWEEN THE MODES IS EQUAL/SYMMETRICAL.
% 
% 	The angle calculated is not that of the mean vector, but that of
% 	the mean axis (or diameter) which passes through both modes of
% 	data.  Thus angle = (a, a+pi).  For this reason, a dataset with
% 	more than one mode is said to be "axial".

%   Along similar lines, if only a semi-circle of angles have been sampled,
%   0-180 degs, let's say, in the belief that the response is bimodal and
%   that the other semi-circle is redundant, the amount of modes must be
%   two for the calculations here, and also (and relatedly) to ensure that
%   this under sampling does not bias the resultant mean vector.
% 
% 	-- Conf Int --
% 	Calculation employs an approximation (see Zar) that is accurate
% 	only to one degree and that is based on the von Mises.
% 
% 	Where a NaN is returned, this is most likely due to the q value
% 	being negative such that there is no square root.  This is
% 	equates to the confidence interval being greater than 90 degrees
% 	(ie, arccos of a negative is between 90 & 270 degrees).
% 	Analytically, it can be shown that this will occur when
% 	l < (chi / 2*n)**0.5
% 
% 	References:
% 		Batschelet (1981)
% 		Zar, Biostatistical Analysis (2010)


	m = modes;
    
    %convert all to rows
    
    if ~isrow(theta)
        theta = theta';
    end
    
    if ~isrow(r)
        r = r';
    end
    
    
    
    % convert all angles to range of 0-180, as orientation is axial.
    % Also converts negative angles to positive, as mod (-x,y) will return
    % y-x when x<y.  Works for angles.
    if max(theta) > 2*pi % convert to radians if not already.
        theta = theta./180 * pi;
    end
    
     theta = mod(theta, pi);
	theta_m = mod(theta.*m, 2*pi);

    %Calculation of the mean sin and cos components of the sum of the
    %vectors, as a ratio of the sum of the lengths of the vectors.
    y = sum(r.* sin(theta_m)) / sum(abs(r)); 
    x = sum(r.* cos(theta_m)) / sum(abs(r));
    
    % If magnitude of the vector is negative, it will take the opposite
    % direction (+180 degs).  As a model of inhibition, this is taken to
    % make sense
    % the magnitude is made absolute though as the divisor as the division
    % here is employed below in order to calculate the circular variance, where the magnitude of the vector,
    % irrespective of its direction, must be positive (as it is always
    % contributing to the variance).
    
    

    % Calculation & Transformations of the mean vector (a or angle)
    % We've multiplied the vector angles by the amount of expected modes.
    % This needs to be factored out.  Also, arctan (atan) returns angles
    % only for the fourth and first quadrants.  Where the real angle
    % (before factoring out the modes) is in the third or second quads,
    % this must be accounted for by translating the angle.
    
    % First quad
    if x > 0 && y > 0
		a = 1/m * atan(y/x); % simply factor by number of modes
        
    % Second quad    
    elseif x < 0 && y > 0
		a = 1/m * (atan(y/x) + pi); % atan gives fourth quad - translate to second then factor
        
    % Third quad   
    elseif x < 0 && y < 0
		a = 1/m * (atan(y/x) + pi); % atan gives first quad ? translate to third and factor
    
    % Fourth quad
    elseif x > 0 && y < 0
		a = 1/m * (atan(y/x) + 2*pi); %atan gives fourth in negative rads ? make positive and factor
    
    % x or cos is zero
    elseif x == 0
        a = 1/m * atan(y/x);
     
    % y or sin is zero
    elseif y == 0
        if x > 0 % angle is 0 degrees/radians, atan function sufficient, no factoring, as zero
            a = atan(y/x);
        elseif x < 0 % angle is 180 degrees / pi radiand - translate atan return and factor
            a = 1/m * (atan(y/x) + pi);
        elseif x == 0 % AND x or cos is zero ... no mean vector?
            error('Cell has absolutely null mean vector.  Dunno what to do');
            
   
        end
    
    % where r is all zeros ... default to pref ori of zero, cv of 1,
    % through assigning x,y,a to all being zero
    elseif isnan(x) && isnan(y)
        x = 0;
        y = 0;
        a = 0;
    end

    
    % Length of the mean vector as a ratio of the sum of the lengths of the
    % component vectors (0 - 1).  This ratio comes from the division in the
    % calculation of `x` and `y`.  The closer this value is to one, the
    % more sharply tuned to orientation it is.
    l = (x^2 + y^2)^0.5;
    
    % Circular Variance (as in Ringach), higher this value, the more
    % `variance` there is in the angles of the component vectors. (0-1)
    CV = 1-l;

    % Preferred angles of the set of vectors.  Accounts for how many modes
    % are expected.  If the calculations are for multiple modes, then an
    % angle is returned for each of these modes.
    % This essentially involves taking the calculated preferred angle and
    % translating it around the circle `m` times, where `m` is the number
    % of modes, and translating one m-th of a circle each time.
    
    
    pref_ori = a + (0:(m-1)).* (2*pi/m);
    
    
end
