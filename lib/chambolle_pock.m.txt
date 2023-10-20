% function x = chambollePock_folkert_diffraction(A, TV, b, maxit, lambda, nonnegative, x0)
function [x, num_of_tviters] = chambollePock_folkert_diffraction(A, TV, b, maxit, lambda, nonnegative, x0)
%CHAMBOLLEPOCK TV-minimization with Chambolle-Pock
%   X = CHAMBOLLEPOCK(A,TV,B) aims at solving:
%      minimize_x ||A*x-B|| + \lambda* ||TV(x)||_1
%   Either A and/or TV can be function handles that operate as:
%      A(x,1)  -  returns A*x
%      A(x,2)  -  returns A'*x
%   As an example one can use the following matrix for the TV operator:
%      n = size(A,2);
%      D = spdiags([[-ones(n-1,1);0],ones(n,1)], [0,1], n, n);
%      TV = [kron(speye(n), D); kron(D, speye(n))];
%
%   X = CHAMBOLLEPOCK(A,TV,B,MAXIT) stops the algorithm when MAXIT
%   iterations are reached, default MAXIT = 20.
%
%   X = CHAMBOLLEPOCK(A,TV,B,MAXIT,LAMBDA) pass the TV weight LAMBDA. The
%   default is LAMBDA = 10;
%
%   X = CHAMBOLLEPOCK(A,TV,B,MAXIT,LAMBDA,NONNEGATIVE) specifies if
%   nonnegativity constraints should be applied to the solution.
%
%   X = CHAMBOLLEPOCK(A,TV,B,MAXIT,LAMBDA,NONNEGATIVE,X0) pass initial
%   guess X0.
%
%   X = CHAMBOLLEPOCK(A,TV,B,MAXIT,LAMBDA,NONNEGATIVE,X0,VISUALIZE)
%   determines the level of detail in output. If VISUALIZE is 0, no output
%   is given. If it is 1, error estimates are presented. If it is larger,
%   an image of the solution is shown.

%   Copyright 2013, Folkert Bleichrodt (CWI, Amsterdam)
%
%   This code is derived from the article:
%   [1] Emil Y. Sidky, Jakob H. Jorgensen and Xiaochuan Pan, "Convex
%   optimization problem prototyping for image reconstruction in computed
%   tomography with the Chambolle-Pock algorithm", 2012
%   Download from: \url{http://arxiv.org/abs/1111.5632}


% autostop = 1;

if nargin < 4 || isempty(maxit)      , maxit  = 0;                      end
if maxit == 0                        , autostop = 1; maxit = 200; else; autostop = 0; end
if nargin < 5 || isempty(lambda)     , lambda = 10;                     end
if nargin < 6 || isempty(nonnegative), nonnegative = false;             end
if nargin < 7 || isempty(x0)         , x0 = 0*A'*b;                     end

% Check for explicit matrices
% Note it is slightly dangerours to assume that "is not a function_handle"
% implies "is numeric", however, we do this to support Spot, without doing
% so explicitly ;).
explicitA  = ~isa(A,  'function_handle');
explicitTV = ~isa(TV, 'function_handle');

% TODO: implement function handles!
% if explicitA
%     x0 = 0*A'*b;
% else
%     x0 = 0*A(b,2)';
% end


% initial vectors
u = x0;
p = zeros(size(b,1),1);
q = zeros(size(TV,1),1);

n = numel(x0);
dim = size(TV,1) / size(TV,2);

L = tvmin_powerMethod(A, TV); % computes the largest singular value L of the matrix [A,TV].
% L = 10; % for a given illumination scenario and TV operator, L is more or less constant
% for circular: 10, for spiral: 6.5
% According to theory, L should be L2 norm of A;

% parameters
tau = 1/L;
sigma = 1/L;
theta = 1;
saturation_level = .1e-3;
countit = 0;

% fprintf('Iter\n');
% fprintf('-------\n');

u_bar = u;
num_iter_autostop = 10;

for i = 1:maxit
    num_of_tviters = i;
    p = (p + sigma*(A*u_bar - b))/(1+sigma);
    q = q + sigma*TV*u_bar;
  
    q = (lambda*q)./max(lambda,abs(q));

    if nonnegative
        u_new = max(u - tau*A'*p - tau*TV'*q,0);
    else
        u_new = u - tau*A'*p - tau*TV'*q;
    end

    dist(i) = norm(u_bar(:) - u_new(:),2) / norm(u_bar(:),2) ;
    
    u_bar = u_new + theta*(u_new-u);
    u = u_new;
    
    if i>num_iter_autostop
        % autostop criterion - if relative change dynamics between
        % consecutive iterations drops below 'saturation_level' - stop
        prog(i) = abs(dist(i)-dist(i-1))/abs(median(dist(3:num_iter_autostop)));
    end
    
    if i<2*num_iter_autostop
        clc
        fprintf('OBJECT SUPPORT GENERATION [%03d]', i);
    else
        clc
        fprintf('OBJECT SUPPORT GENERATION [%03d]; CURRENT LEVEL [%.5f] >= SATURATION LEVEL [%.5f] \n', i, median(prog(end-num_iter_autostop-1:end)), saturation_level);
    end
    
%     if i>num_iter_autostop && autostop && (prog(i) < saturation_level) && (prog(i-1) < saturation_level)
% 		countit = countit+1;
%     elseif i>num_iter_autostop && autostop
%         countit = 0;
%     end
%     if countit == 10
% 		fprintf('Stop on MSE saturation.\n')
% 		i = maxit;
% 		break;
%     end
    
    if i>2*num_iter_autostop && autostop && (median(prog(end-num_iter_autostop-1:end)) < saturation_level)
% 		fprintf('Stop on MSE saturation.\n')
		i = maxit;
		break;
    end
end

x = u;

end


function L = tvmin_powerMethod(W, D)
%TVMIN_POWERMETHOD Power method for TV-min tomo system matrix
%   L = TVMIN_POWERMETHOD(W, D) computes the largest singular value L of 
%   the matrix [W,D].

niter = 5;

% random initial image
x = rand(size(W,2),1);
y = W*x;
z = D*x;

for i = 1:niter
    % power iteration
    x = W'*y + D'*z;
    % normalize
    x = x./norm(x);
    % One less matvec per iteration,
    % but double the memory usage.
    y = W*x;
    z = D*x;
    L = sqrt(y'*y + z'*z);
end

end

