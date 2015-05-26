function [x_nozzle, y_nozzle] = nozzle_geo(y_throat, theta_max_nozzle,n)

 %y_initial = 0.1307; %meters -- radius of throat from axis of symmetry (r_th)
 %x_initial = 0;
 %n = number of points
 %theta_max_nozzle = 35;     %degrees -- maximum wall angle
 
 
 r_curve = y_throat/2; %meters -- radius of fixed wall
 
 theta_nozzle = linspace( 270, 270+theta_max_nozzle, n );
 
 x_cosd = cosd(theta_nozzle);
 y_sind = sind(theta_nozzle);
 
 x_nozzle = x_cosd * r_curve;
 y_nozzle = r_curve * ( y_sind - y_sind(1) ) + y_throat;
 
end
