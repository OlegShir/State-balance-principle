function value = filtLaplace(typeDistribution, varargin)
%%Функция возвращает гипердельтную аппроксимацию произвольной
% плотности распределения в приближённом виде.
%
%   For input varagin using:
%       - exponential distribution: lambda;
%       - gamma distribution: k, lambda;
%       - normal distribution: mu, sigma;
%       - uniform distribution: a;
%       - rayleigh distribution sigma:

% check input parameters in function
if (typeDistribution == "Exponential" && nargin == 2)||...
   (typeDistribution == "Gamma" && nargin == 3)||...
   (typeDistribution == "Normal" && nargin == 3)||...
   (typeDistribution == "Uniform" && nargin == 2)||...
   (typeDistribution == "Rayleigh" && nargin == 2)
    
    switch typeDistribution
        
        case 'Exponential'
            % C1  = 0.854; C2 = 0.146; T1 = 0.293/lambda; T2 = 1.707/lambda
            C1 = 0.854;
            C2 = 0.146;
            T1 = 0.293/varargin{1};
            T2 = 1.707/varargin{1};
            
        case 'Gamma'
            % C1  = (1 + sqrt(k + 1))/2*sqrt(k + 1); C2 = (1 - sqrt(k + 1))/2*sqrt(k + 1);
            % T1 = (k + 1 - sqrt(k + 1))/lambda; T2 = (k + 1 + sqrt(k + 1))/lambda
            C1 = (1 + sqrt(varargin{1} + 1))/2*sqrt(varargin{1} + 1);
            C2 = (1 - sqrt(varargin{1} + 1))/2*sqrt(varargin{1} + 1);
            T1 = (varargin{1} + 1 - sqrt(varargin{1} + 1))/varargin{2};
            T2 = (varargin{1} + 1 + sqrt(varargin{1} + 1))/varargin{2};
            
        case 'Normal'
            % C1  = 0.5; C2 = 0.5; T1 = mu - sigma; T2 = mu + sigma
            C1 = 0.5;
            C2 = 0.5;
            T1 = varargin{1} - varargin{2};
            T2 = varargin{1} + varargin{2};
            
        case 'Uniform'
            % C1  = 0.5; C2 = 0.5;
            % T1 = 0.211*a; T2 = 0.789*a
            C1 = 0.5;
            C2 = 0.5;
            T1 = 0.211*varargin{1};
            T2 = 0.789*varargin{1};
            
        case 'Rayleigh'
            % C1  = 0.35; C2 = 0.65; T1 = 0.773*sigma; T2 = 2.147*sigma
            C1 = 0.35;
            C2 = 0.65;
            T1 = 0.773*varargin{1};
            T2 = 2.147*varargin{1};
            
    end
    
    syms t s
    % build approximate function a(t) = C1*dir(t-T1)+C2*dir(t-T2)
    f(t) = C1 * dirac(t - T1) + C2 * dirac(t - T2);
    % get Laplace's image
    f(s) = laplace(f, t, s);
    % get Laplace's filter
    value = f(1/t);
    
else
    disp(['Введено некоректное число моментов для ', typeDistribution, ' распределения'])
    value = [];
    
end
end
