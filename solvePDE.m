function solvePDE
% PROTOTYPE AS IT IS
clc
le = input('*[1]Load* or [2]Use matrices or [3]Enter w/ data ? ');
if isempty(le)
    le = 1;
end

mcint = @(F,a,b) prod(b-a)*(1/length(F))*sum(F);
xi_x = @(x,a,b) 2*(x-a)/(x-b) -1;
x_xi = @(xi,a,b) (b-a)*xi/2 + (b+a)/2;

if le == 1
    load('backup_lastPDE');
elseif le == 2
    nvar = 2;
    order = 2;
    a = [0 0];
    b = [1 1];
    s_x = @(x)0;
    pword = up2Power(order,nvar);
    [Q,~] = size(pword);
    kq = {@(x)0, @(x)0, @(x)0, @(x)1, @(x)0, @(x)1};
    ncond = 20;
    ordercond = 0;
    pwordcond = up2Power(ordercond,nvar);
    [W,~] = size(pwordcond);
    hc = cell(ncond,W);
    hc{1,1} = @(x) 1; 
    hc{2,1} = @(x) 1; 
    hc{3,1} = @(x) 1; 
    hc{4,1} = @(x) 1; 
    hc{5,1} = @(x) 1; 
    hc{6,1} = @(x) 1; 
    hc{7,1} = @(x) 1; 
    hc{8,1} = @(x) 1; 
    hc{9,1} = @(x) 1; 
    hc{10,1} = @(x) 1; 
    hc{11,1} = @(x) 1; 
    hc{12,1} = @(x) 1; 
    hc{13,1} = @(x) 1; 
    hc{14,1} = @(x) 1; 
    hc{15,1} = @(x) 1; 
    hc{16,1} = @(x) 1; 
    hc{17,1} = @(x) 1; 
    hc{18,1} = @(x) 1; 
    hc{19,1} = @(x) 1; 
    hc{20,1} = @(x) 1; 
    ptc = [0,0;0,0.2;0,0.4;0,0.6;0,0.8;0,1;0.2,0;0.4,0;0.6,0;0.8,0;...
        1,0;1,0.2;1,0.4;1,0.6;1,0.8;1,1;0.2,1;0.4,1;0.6,1;0.8,1];
    Vc = [0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;sin(0.2*pi);sin(0.4*pi);...
        sin(0.6*pi);sin(0.8*pi)];
    knownsolution = @(x) sin(pi*x(1,:)).*sinh(pi*x(2,:))./sinh(pi);
%     nvar = 2;
%     order = 2;
%     a = [0 0];
%     b = [1 0.01];
%     s_x = @(x) 0;
%     pword = up2Power(order,nvar);
%     [Q,~] = size(pword);
%     kq = {@(x)0, @(x)0, @(x)-1, @(x)1, @(x)0, @(x)0};
%     ncond = 12;
%     ordercond = 1;
%     pwordcond = up2Power(ordercond,nvar);
%     [W,~] = size(pwordcond);
%     hc = cell(ncond,W);
%     hc{1,1} = @(x) 0; hc{1,2} = @(x) 1; hc{1,3} = @(x) 0; 
%     hc{2,1} = @(x) 0; hc{2,2} = @(x) 1; hc{2,3} = @(x) 0;
%     hc{3,1} = @(x) 0; hc{3,2} = @(x) 1; hc{3,3} = @(x) 0;
%     hc{4,1} = @(x) 0; hc{4,2} = @(x) 1; hc{4,3} = @(x) 0;
%     hc{5,1} = @(x) 0; hc{5,2} = @(x) 1; hc{5,3} = @(x) 0;
%     hc{6,1} = @(x) 0; hc{6,2} = @(x) 1; hc{6,3} = @(x) 0;
%     hc{7,1} = @(x) 0; hc{7,2} = @(x) 1; hc{7,3} = @(x) 0;
%     hc{8,1} = @(x) 0; hc{8,2} = @(x) 1; hc{8,3} = @(x) 0;
%     hc{9,1} = @(x) 1; hc{9,2} = @(x) 0; hc{9,3} = @(x) 0;
%     hc{10,1} = @(x) 1; hc{10,2} = @(x) 0; hc{10,3} = @(x) 0;
%     hc{11,1} = @(x) 1; hc{11,2} = @(x) 0; hc{11,3} = @(x) 0;
%     hc{12,1} = @(x) 1; hc{12,2} = @(x) 0; hc{12,3} = @(x) 0;
%     ptc = [0,0;0,1/300;0,2/300;0,0.01;1,0;1,1/300;1,2/300;1,0.01;...
%         0,0;1/3,0;2/3,0;1,0];
%     Vc = [0;0;0;0;0;0;0;0;2.5;1;1;-0.5];
%     knownsolution = @(x) 1+exp(-pi^2*x(2,:)).*cos(pi*x(1,:))+...
%         0.5*exp(-(3*pi)^2*x(2,:)).*cos(3*pi*x(1,:));
else
    nvar = input('How many variables? ');
    order = input('Order of differential? ');
    uxstr = 'u(';
    a = zeros(1,nvar);
    b = zeros(1,nvar);
    for j = 1:nvar
        a(j) = input(['Which is the minimum for x_' num2str(j) '? ']);
        b(j) = input(['Which is the maximum for x_' num2str(j) '? ']);
        uxstr = [uxstr 'x_' num2str(j)];
        if j~= nvar
            uxstr = [uxstr ','];
        end
    end
    uxstr = [uxstr ')'];
    aux = input('Source function s(x) = ');
    if isa(aux,'double')
        s_x = @(x) aux;
        s_x_str = num2str(aux);
    else
        s_x = aux;
        s_x_str = func2str(aux);
        s_x_str = s_x_str(5:end);
    end
    pword = up2Power(order,nvar);
    [Q,~] = size(pword);
    kq = cell(1,Q);
    DEq = [];
    for i = 1:Q
        diff_str = [];
        for j = 1:nvar
            if pword(i,j) ~= 0
                diff_str = [diff_str 'ddx_' num2str(j) '^' ...
                    num2str(pword(i,j)) ' '];
            end
        end
        aux = input(['Coeff of [' diff_str ']' uxstr ' ? ']);
        chk0 = false;
        if isa(aux,'double')
            if aux == 0
                chk0 = true;
            end
            kq{i} = @(x) aux;
            kq_str = num2str(aux);
        else
            kq{i} = aux;
            kq_str = func2str(aux);
            kq_str = kq_str(5:end);
        end
        if ~chk0
            if(~isempty(DEq))
                DEq = [DEq ' + '];
            end
            DEq = [DEq kq_str '*[' diff_str '] ' uxstr ];
        end
    end
    DEq = [DEq ' = ' s_x_str];
    disp(DEq)
    ncond = input('How many boundary conditions? ');
    ordercond = input(['Maximum differential order for '...
        'boundary conditions? ']);
    if ordercond >= order
        ordercond = order - 1;
        warning(['New maximum order for conditions is: ' ...
            num2str(ordercond)]);
    end
    pwordcond = up2Power(ordercond,nvar);
    [W,~] = size(pwordcond);
    hc = cell(ncond,W);
    ptc = zeros(ncond,nvar);
    Vc = zeros(ncond,1);
    for c = 1:ncond
        DEq = [];
        for i = 1:size(pwordcond,1)
            diff_str = [];
            for j = 1:nvar
                if pword(i,j) ~= 0
                    diff_str = [diff_str 'ddx_' num2str(j) '^' ...
                        num2str(pword(i,j)) ' '];
                end
            end
            aux = input(['Cond#' num2str(c) '; coeff of [' diff_str ...
                ']' uxstr ' ? ']);
            chk0 = false;
            if isa(aux,'double')
                if aux == 0
                    chk0 = true;
                end
                hc{c,i} = @(x) aux;
                hc_str = num2str(aux);
            else
                hc{c,i} = aux;
                hc_str = func2str(aux);
                hc_str = hc_str(5:end);
            end
            if ~chk0
                if(~isempty(DEq))
                    DEq = [DEq ' + '];
                end
                DEq = [DEq hc_str '*[' diff_str '] ' uxstr ];
            end
        end
        ptc(c,:) = input('At which point? ');
        Vc(c) = input('Which known value? ');
        DEq = [DEq ' @( ' num2str(ptc(c,:)) ') = ' num2str(Vc(c))];
        disp(['cond#' num2str(c) ' : ' DEq]);
    end
    knownsolution = input('Known solution: ');
    save('backup_lastPDE');
end

degree = input('Which degree of approximation? ');
while true
    pwdeg = up2Power(degree,nvar);
    [N,~] = size(pwdeg);
    if N > ncond
        break;
    else
        degree = degree + 1;
    end
end
disp(['Using degree ' num2str(degree) ' for better adjustment.'])
mcnpts = input('How many points to Monte Carlo? ');
tolnum = 1e-12; %degree^2*eps;

fprintf('\nGalerkin\n')
tic
Oxxx = linspace(-1,1,floor(mcnpts^(1/nvar)));
% ========
dimension = nvar;
min_count = 1;
max_count = length(Oxxx);
counter = min_count.*ones(1,dimension);
xxx = [];
while true
    aux = zeros(nvar,1);
    for v = 1:nvar
        aux(v) = Oxxx(counter(v));
    end
    xxx = cat(2,xxx,aux);
    counter(1) = counter(1) + 1;
    for i = 1:length(counter)-1
        if counter(i) > max_count
            counter(i) = min_count;
            counter(i+1) = counter(i+1) +1;
        end
    end
    if counter(dimension) > max_count
        break;
    end
end
% ========
BAQ = zeros(1,size(pword,1));
for i = 1:length(BAQ)
    BAQ(i) = 1;
    for v = 1:nvar
        BAQ(i) = BAQ(i)*1/(b(v)-a(v))^pword(i,v);
    end
end
A = zeros(N,N);
bv = zeros(N,1);
for n = 1:N
    disp([num2str(n) ' / ' num2str(N)])
    iKPP = zeros(Q,N);
    Pn = ones(1,size(xxx,2));
    coord = zeros(size(xxx));
    for v = 1:nvar
        Pn = Pn.*jacobiPolynomial(pwdeg(n,v),0,0,xxx(v,:),-1,1);
        coord(v,:) = x_xi(xxx(v,:),a(v),b(v));
    end
    for q = 1:Q
        for m = 1:N
            Pqm = ones(1,size(xxx,2));
            for v = 1:nvar
                Pqm = Pqm.*(gamma(pwdeg(m,v)+pword(q,v)+1)/...
                    (gamma(pwdeg(m,v)+1))).*...
                    jacobiPolynomial(pwdeg(m,v)-pword(q,v), ...
                    pword(q,v),pword(q,v),xxx(v,:),-1,1);
            end
            iKPP(q,m) = mcint(kq{q}(coord).*Pn.*Pqm, -1, 1);
        end 
    end
    A(n,:) = BAQ*iKPP;
    bv(n,:) = mcint(Pn.*s_x(coord),-1,1);
end
toc
disp(['Only Galerkin: size ' num2str(size(A)) ' rank ' num2str(rank(A))...
    ' cond ' num2str(cond(A))])

fprintf('\nConditions\n')
tic
BAQ = zeros(1,size(pwordcond,1));
for i = 1:length(BAQ)
    BAQ(i) = 1;
    for v = 1:nvar
        BAQ(i) = BAQ(i)*1/(b(v)-a(v))^pwordcond(i,v);
    end
end
for c = 1:ncond
    disp([num2str(c) ' / ' num2str(ncond)])
    i = N-c+1;
    coord = ptc(c,:)';
    CD = zeros(W,N);
    for w = 1:W
        for m = 1:N
            Pwm = 1;
            for v = 1:nvar
                Pwm = Pwm.*gamma(pwdeg(m,v)+pwordcond(w,v)+1)/...
                    gamma(pwdeg(m,v)+1).* ...
                    jacobiPolynomial(pwdeg(m,v)-pwordcond(w,v),...
                    pwordcond(w,v),pwordcond(w,v),coord(v),a(v),b(v));
            end
            CD(w,m) = hc{c,w}(coord)*Pwm;
        end
    end
    A(i,:) = BAQ*CD;
    bv(i) = Vc(c);
end
toc
disp(['Galerkin + Conditions: size ' num2str(size(A)) ' rank ' ...
    num2str(rank(A)) ' cond ' num2str(cond(A))])

fprintf('\nTolerance adjustments...\n')
A(abs(A)<tolnum) = 0;
bv(abs(bv)<tolnum) = 0;
disp(['Galerkin + Conditions after tolerance: size ' ...
    num2str(size(A)) ' rank ' num2str(rank(A)) ' cond ' ...
    num2str(cond(A))])

U = A\bv;
U(abs(U)<tolnum) = 0;

if rcond(A) <= eps
    disp('Try raising polynomial degree...')
end

if nvar == 1
    ox = linspace(a(1),b(1),500)';
    oy = u_hat(U,ox,a(1),b(1));

    subplot(3,1,[1,2])
    hold on
    if(~isempty(knownsolution))
        plot(ox,knownsolution(ox),'Color',[0.75,0.75,0.75],'LineWidth',4)
    end
    plot(ox,oy,'--k','LineWidth',2)
    legend('Known solution','Approximation')
    hold off
    title(['\bfUsing degree ' num2str(degree) ...
        ' for polynomial approximation; and ' num2str(mcnpts) ...
        ' points for Monte Carlo integration'],'FontSize',14)
    set(gca,'FontSize',12);
    grid on
    subplot(3,1,3)
    if(~isempty(knownsolution))
        yyy = oy-knownsolution(ox);
        plot(ox,yyy,'-','Color',[0.5,0.5,0.5],'LineWidth',2)
    end
    title(['\bfError [ max abs / mean / std ]:   ' ...
        num2str(max(abs(yyy)),'%e') ' / ' num2str(mean(yyy),'%e') ' / '...
        num2str(std(yyy),'%e')],'FontSize',14)
    set(gca,'FontSize',12);
    grid on
elseif nvar == 2
    [ox,oy] = meshgrid(linspace(a(1),b(1),20),linspace(a(2),b(2),20));
    OX = [ox(:), oy(:)]';
    OZ = u_hat(U,OX,a,b)';
    oz = reshape(OZ,size(ox));

    subplot(3,1,[1,2])
    hold on
    if(~isempty(knownsolution))
        mesh(ox,oy,reshape(knownsolution(OX)',size(ox)),'EdgeColor',...
            'none','FaceColor',[0.75,0.75,0.75],'LineStyle','-',...
            'LineWidth',2)
    end
    mesh(ox,oy,oz,'LineStyle',':','LineWidth',1),colormap(zeros(64,3))
    legend('Known solution','Approximation')
    hold off
    title(['\bfUsing degree ' num2str(degree) ...
        ' for polynomial approximation; and ' num2str(mcnpts) ...
        ' points for Monte Carlo integration'],'FontSize',14)
    set(gca,'FontSize',12);
    grid on
    view(15,30)
    subplot(3,1,3)
    if(~isempty(knownsolution))
        yyy = (OZ-knownsolution(OX))';
        mesh(ox,oy,reshape(yyy,size(ox)),'LineStyle','-','LineWidth',1)
        colormap('gray')
    end
    title(['\bfError [ max abs / mean / std ]:   ' ...
        num2str(max(abs(yyy)),'%e') ' / ' num2str(mean(yyy),'%e') ' / '...
        num2str(std(yyy),'%e')],'FontSize',14)
    set(gca,'FontSize',12);
    grid on
    view(15,30)
else
    disp('No plot available');
end
end
%%
function y = u_hat(uk,x,a,b)
if size(x,1) < size(x,2)
    x = x';
end
nv = length(a);
for d = 0:100
    if factorial(nv+d)/(factorial(nv)*factorial(d)) >= length(uk)
        break;
    end
end
if d < 20 && nv < 3
    mypowers = up2Power(d,nv);
    y = zeros(size(x,1),1);
    for i = 1:length(uk)
        aux = ones(size(x,1),1);
        for v = 1:nv
            aux = aux.*jacobiPolynomial(mypowers(i,v),0,0,x(:,v),a(v),...
                b(v));
        end
        y = y + uk(i).*aux;
    end
else
    disp(['Number of variables: ' num2str(nv) ', dimension: ' num2str(d)])
    disp('Error')
end
end
%%
function pows = up2Power(deg,nvar)
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
    nx = 2*(x-a)/(b-a) - 1;
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
