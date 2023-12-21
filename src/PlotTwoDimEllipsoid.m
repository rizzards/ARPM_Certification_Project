function [ell_handle, ell_points, ax1, ax2] = PlotTwoDimEllipsoid(mu, sigma2, r, PlotAxes, PlotTangBox, color, linewidth, PlotEll, n_points)
% This function creates and plots the two dimensional ellipsoid: 
% (x - mu)' * (sigma2^(-1)) * (x - mu) = r^2 
%  INPUTS 
%   mu          : [vector] (2 x 1) 
%   sigma2      : [matrix] (2 x 2) symmetric and positive definite matrix
%   r           : [scalar] radius of the ellipsoid
%   PlotAxes    : [boolean] if true then the principal axes are plotted
%   PlotTangBox : [boolean] if true then the tangent box is plotted.
%   color       : [char] color of the line defining the ellipsoid
%   linewidth   : [scalar] width of the line defining the ellipsoid
%   PlotEll     : [boolean] if true then the ellipsoid is plotted
%   n_points    : [scalar] number of points of the ellipsoid
%  OUTPUTS
%   ell_handle  : [figure handle] ellipsoid
%   ell_points  : [matrix] (2 x n_points) points of the ellipsoid 
%   ax1, ax2    : [figure handle] principal axes

if nargin < 9 || isempty(n_points); n_points = 1000; end
if nargin < 8 || isempty(PlotEll); PlotEll = true; end
if nargin < 7 || isempty(linewidth); linewidth = 2; end
if nargin < 6 || isempty(color); color = 'k'; end
if nargin < 5 || isempty(PlotTangBox); PlotTangBox = false; end
if nargin < 4 || isempty(PlotAxes); PlotAxes = false; end
if nargin < 3 || isempty(r); r = 1; end

%% Code                  
theta = linspace(0, 2 * pi, n_points);

% compute the initial sphere
y = [r * cos(theta); r * sin(theta)];

% principal axes
y_axes1 = [-r, r; 0, 0];
y_axes2 = [0, 0; -r, r];

% spectral decomposition of sigma2
[e, Diag_lambda2] = eig(sigma2);
lambda2 = diag(Diag_lambda2);
[lambda2, order] = sort(lambda2, 'descend');
e = e(:, order);
Diag_lambda = diag(sqrt(lambda2));

% compute the ellipsoid as affine transformation of the sphere
u = e * Diag_lambda * y;
u_axes1 = e * Diag_lambda * y_axes1;
u_axes2 = e * Diag_lambda * y_axes2;
ell_points = repmat(mu, 1, n_points) + u; 
      

% plot the ellipsoid
if PlotEll
    hold on;
    ell_handle = plot(ell_points(1, :), ell_points(2, :));
    set(ell_handle, 'color', color, 'linewidth', linewidth);
    grid on;
else
    ell_handle = [];
end

% plot the tangent box
if PlotTangBox
    sigvec = sqrt(diag(sigma2));
    
    tangBox_low = [mu(1) - r * sigvec(1),  mu(1) + r * sigvec(1); mu(2) - r * sigvec(2), mu(2) - r * sigvec(2)];
    tangBox_up = [mu(1) - r * sigvec(1),  mu(1) + r * sigvec(1); mu(2) + r * sigvec(2), mu(2) + r * sigvec(2)];
    tangBox_left = [mu(1) - r * sigvec(1), mu(1) - r * sigvec(1); mu(2) - r * sigvec(2),  mu(2) + r * sigvec(2)];
    tangBox_right = [mu(1) + r * sigvec(1), mu(1) + r * sigvec(1); mu(2) - r * sigvec(2),  mu(2) + r * sigvec(2)];
            
    hold on;
    h1 = plot(tangBox_low(1,:), tangBox_low(2,:));
    hold on;
    h2 = plot(tangBox_up(1,:), tangBox_up(2,:));
    hold on;
    h3 = plot(tangBox_left(1,:), tangBox_left(2,:));
    hold on;
    h4 = plot(tangBox_right(1,:), tangBox_right(2,:));
    set([h1 h2 h3 h4], 'color', color, 'linewidth',  linewidth);
end
        
% plot the principal axes
if PlotAxes
    ax1 = plot(u_axes1(1,:) + mu(1), u_axes1(2,:) + mu(2));
    hold on;
    ax2 = plot(u_axes2(1,:) + mu(1), u_axes2(2,:) + mu(2));
    set([ax1 ax2], 'color', color, 'linewidth', linewidth);
    axis equal;
end
