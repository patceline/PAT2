%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% PAT^2 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Propagator of asteroid (and planetary) ephemerides. 
% Default bodies are: Sun, Mercury, Venus, Earth, Mars, Jupiter,
% Saturn, Uranus, Neptune, Pluto, the Moon and any asteroid. 
% If more or less bodies are to be modelled, change n (n is the number 
% of bodies in the model).
% I recommend using barycentric initial conditions (state vector and 
% mass), but body center values will also be fine.
% Initial conditions are required for all bodies: positions, velocities, 
% as well as mass. The units must be km, s, and kg. The masses of the
% barycenters of the Sun, planets and Moon are implemented, but these can
% be modified if desired. For the asteroid, you will need to provide the
% state vector and the mass.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Last modified June 19 2014,
% pat_egger@hotmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Y,runtime] = PAT2

global n beta gamma c mu
global corr year

n = 15;
beta = 1.0;
gamma = 1.0;

% re-scale the problem to avoid roundoff 
corr = 1e5; % 10^5 km
year = 365*(23*60*60+56*60+4);% 1 year in seconds
c = year*299792.458/corr;

% set up the problem (initial conditions, time span)
[mu,mass,actual_position,y0,tspan] = setup2; 

% set up the initial state vector (with new units)
y0_tilde(1:3*n,1) = y0(1:3*n,1)/corr;
y0_tilde(3*n+1:6*n,1) = year*y0(3*n+1:6*n,1)/corr;

       
% choose tolerances
Abstol = 1e-4; % 1e-6 by default
Reltol = 1e-4; % 1e-3 by default

% set up integration options
options = odeset('Mass',@mass_matrix,'MStateDependence','strong',...
    'MaxStep',1,'Abstol',Abstol,'Reltol',Reltol,'Stats','on');

%if strcmp(equation,'relativistic')

% numerical integration    
tic;
    [~,Y_tilde] = ode113(@equ_relativistic,tspan,...
                    y0_tilde,options);
    runtime = toc;
%elseif strcmp(equation,'newtonian')
  %  tic;
 %   [~,Y_tilde] = ode113(@equ_newtonian,tspan,...
                   % y0_tilde,options);
 %   runtime = toc;
%end;
% scale the ouput back to km and s       
% pos = Y_tilde(:,1:3)*corr;
% vel = Y_tilde(:,3*n+1:3*n+3)*corr/year;

% compute maximal error and error in the last time step
% plots the error in norm and component-wise
%[max_error_pos,last_error] = plot_error(pos,actual_position)

% output runtime in minutes
runtime = runtime/60;

Y = zeros(size(Y_tilde));

Y(:,1:3*n) = Y_tilde(:,1:3*n)*corr;
Y(:,3*n+1:end) = Y_tilde(:,3*n+1:end)*corr/year;

end





