%%
function y = u_hat(uk,x)
nv = length(x);
for d = 0:100
    if factorial(nv+d)/(factorial(nv)*factorial(d)) >= length(uk)
        break;
    end
end
if d < 20 && nv < 3
    mypowers = getPowers(d,nv);
    y = zeros(size(x{1}));
    for i = 1:length(uk)
        aux = ones(size(x{1}));
        for v = 1:nv
            aux = aux.*jacobiPolynomial(mypowers(i,v),0,0,x{v});
        end
        y = y + uk(i).*aux;
    end
else
    disp(['Number of variables: ' num2str(nv) ', dimension: ' num2str(d)])
    disp('Error')
end
end
%%
% function y = u_hat(uk,x,a,b)
% if exist('a','var') == 0
%     a = -1;
% end
% if exist('b','var') == 0
%     b = 1;
% end
% if size(x,1) < size(x,2)
%     x = x';
% end
% nv = length(a);
% for d = 0:100
%     if factorial(nv+d)/(factorial(nv)*factorial(d)) >= length(uk)
%         break;
%     end
% end
% if d < 20 && nv < 3
%     mypowers = up2Power(d,nv);
%     if iscell(x)
%         y = zeros(size(x{1})); 
%     else
%         y = zeros(size(x,1),1);
%     end
%     for i = 1:length(uk)
%         aux = ones(size(y));
%         for v = 1:nv
%             if iscell(x)
%                 jp = jacobiPolynomial(mypowers(i,v),0,0,x{v},a(v),...
%                     b(v));
%             else
%                 jp = jacobiPolynomial(mypowers(i,v),0,0,x(:,v),a(v),...
%                     b(v));
%             end
%             size(jp)
%             aux = aux.*jp;
%         end
%         y = y + uk(i).*aux;
%     end
% else
%     disp(['Number of variables: ' num2str(nv) ', dimension: ' num2str(d)])
%     disp('Error')
% end
% end
%%
function pows = getPowers(deg,nvar)
pows = [];
for i = 0:deg
    pows = cat(1,pows,Powers(i,nvar));
end
end
%%
function pows = Powers(deg,nvar)
X = integer_partitions(deg,nvar);
pows = [];
for i = 1:size(X,1)
    pows = cat(1,pows,unique(perms(X(i,:)),'rows'));
end
pows = flipud(sortrows(pows));
end
%%
function S = integer_partitions(n,count)
%
% This work is licensed under a 
% Creative Commons Attribution 4.0 International License. 
% https://creativecommons.org/licenses/by/4.0/
%
% Coded by IGOR S. PERETTA - iperetta@gmail.com (2015)
%
% Arguments/inputs:
% - n: the desired nonnegative integer number to be partitioned
% - count: the maximum number of integers to build the desired number
%
% Output:
% - S: a matrix with "count" columns and as many lines as necessary to 
%      provide every possible partition of a nonnegative integer "n" 
%      considering a sum of "count" integer numbers. Related permutations 
%      are not considered here (try "help perms" for that).
%
% Examples of usage:
%         >> integer_partitions(5)
%         
%         ans =
%         
%              5     0     0     0     0
%              4     1     0     0     0
%              3     2     0     0     0
%              3     1     1     0     0
%              2     2     1     0     0
%              2     1     1     1     0
%              1     1     1     1     1
%         
%         >> integer_partitions(5,3)
%        
%         ans =
%         
%              5     0     0
%              4     1     0
%              3     2     0
%              3     1     1
%              2     2     1
%        
%         >> integer_partitions(3,6)
%         
%         ans =
%         
%              3     0     0     0     0     0
%              2     1     0     0     0     0
%              1     1     1     0     0     0
%
% Adapted from "Algorithm ZS1" in
% ANTOINE ZOGHBI and IVAN STOJMENOVIC (1998), "FAST ALGORITHMS FOR 
% GENERATING INTEGER PARTITIONS", International Journal of Computer 
% Mathematics, Volume 70, Issue 2, pages 319-332
% DOI: 10.1080/00207169808804755
%
if nargin == 1
    count = n;
end
if n < 0 || n ~= round(n)
    error('Only nonnegative integers allowed!');
elseif n == 0
    if count == 0
        S = 0;
    else
        S = zeros(1,count);
    end
else
    x = ones(1,n);
    x(1) = n;
    m = 1;
    h = 1;
    M = [x(1:m) zeros(1,n-m)];
    while x(1) ~= 1
        if x(h) == 2 
           m = m + 1;
           x(h) = 1;
           h = h - 1;
        else
           r = x(h) - 1;
           t = m - h + 1;
           x(h) = r;
           while t >= r
               h = h + 1;
               x(h) = r;
               t = t - r;
           end
           if t == 0
               m = h;
           else
               m = h + 1;
               if t > 1
                   h = h + 1;
                   x(h) = t;
               end
           end
        end
        M = cat(1,M,[x(1:m) zeros(1,n-m)]);
    end
    if count > n
        M = cat(2,M,zeros(size(M,1),count-n));
    end
    S = [];
    for i = 1:size(M,1)
        if(sum(M(i,1:count)) == n)
            S = cat(1,S,M(i,1:count));
        end
    end
end
end
%%
function y = jacobiPolynomial( n, alpha, beta, x, a, b )
if exist('a','var') == 0;
    a = -1;
end
if exist('b','var') == 0;
    b = 1;
end
if n == 0
    y = ones(size(x));
elseif n < 0
    y = zeros(size(x));
else
    if a == -1 && b == 1
        if iscell(x)
            nx = x{1};
        else
            nx = x;
        end
    else
        nx = 2*(x-a)/(b-a) - 1;
    end
    coeff = gamma(n+alpha+1)/(gamma(n+1)*gamma(alpha+1));
    expr = 0;
    for k = 0:n
        sumcoeff = pochhammer(-n, k)*pochhammer(n + alpha + beta + 1, k)...
            /pochhammer(alpha + 1, k);
        expr = expr + sumcoeff*(0.5*(1 - nx)).^k/gamma(k + 1);
    end
    y = coeff*expr;
end
end
%%
function P = pochhammer(x,n)
% If n is a negative integer, then the identity pochhammer(x, n) = 
%    1/pochhammer(x + n, -n) is used to express the result.
% The following special cases are implemented: pochhammer(x, 0) = 1, 
%   pochhammer(x, 1) = x, pochhammer(x,-1) = 1/(x - 1), 
%   pochhammer(1, n) = gamma(n + 1), pochhammer(2, n) = gamma(n + 2).
% If n is a positive integer, then expand(pochhammer(x, n)) yields the 
%   expanded polynomial x*(x + 1)*...*(x + n - 1).
% If n is not an integer, then expand(pochhammer(x, n)) yields a 
%   representation in terms of gamma.
if n == 0
    P = 1;
elseif n == 1
    P = x;
elseif n == -1
    P = 1./(x-1);
elseif x == 1
    P = gamma(n+1);
elseif x == 2
    P = gamma(n+2);
elseif n < 0
    P = 1./pochhammer(x + n,-n);
elseif n == round(n) && n > 0
    P = 1;
    for i = 0:n-1
        P = P.*(x + i);
    end
else
    P = gamma(x + n)./gamma(x);
end
end