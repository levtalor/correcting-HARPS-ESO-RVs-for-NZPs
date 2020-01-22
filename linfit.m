function [a,b,da,db,chi2,rms,res,chi2null,Fstatistic,Ftest,h1] = linfit(x,y,dy,varargin)
% [a,b,da,db,chi2,rms,res,chi2null,Fstatistic,Ftest,h1] = linfit(x,y,dy,varargin)
% fits a line to a vectors of data points and errors, ignoring NaN values.
% Optionally, makes plots.
%
% Last Modified: 20200121 LT

x = x(:);
y = y(:);
dy = dy(:);

if nargin>3
    showflag = varargin{1};
else
    showflag = 0;
end

if nargin>4
    saveflag = varargin{2};
    filename = varargin{3};
    save_fmt = varargin{4};
else
    saveflag = 0;
    filename = 'Y_vs_X';
    save_fmt = '.fig';
end

if nargin>7
    xlabel_str = varargin{5};
    ylabel_str = varargin{6};
    title_str = varargin{7};
else
    xlabel_str = 'X';
    ylabel_str = 'Y';
    title_str = 'Analytic line fit';
end

if nargin>10
    dx = varargin{8};
else
    dx = nan(size(x));
end
dx = dx(:);

% remove NaN values:
bad_ind = find(~isfinite(x) | ~isfinite(y) | ~isfinite(dy));
x(bad_ind) = [];
y(bad_ind) = [];
dy(bad_ind) = [];
dx(bad_ind) = [];

% Fit line analytically:
matrix = [ sum(x.^2./dy.^2)       sum(x./dy.^2);
           sum(x./dy.^2)          sum(1./dy.^2)];
y_sigma = y./dy.^2;
free_vector = [sum(x.*y_sigma);
               sum(1.*y_sigma)];
solution=matrix\free_vector;           
a=solution(1);b=solution(2);
err=sqrt(diag(inv(matrix)));
da=err(1);db=err(2);
% if there are errors in x, co-add them and re-fit
if nargin>10
    dy = sqrt((a*dx).^2+dy.^2);
    matrix = [ sum(x.^2./dy.^2)       sum(x./dy.^2);
               sum(x./dy.^2)          sum(1./dy.^2)];
    y_sigma = y./dy.^2;
    free_vector = [sum(x.*y_sigma);
                   sum(1.*y_sigma)];
    solution=matrix\free_vector;           
    a=solution(1);b=solution(2);
    err=sqrt(diag(inv(matrix)));
    da=err(1);db=err(2);
end

% some statistics:
res=y-(b+a*x);
rms=1.48*mad(res,1);
chi2=sum((res./dy).^2);
p1 = 1; 
p2 = 2; 
p2_p1 = p2-p1;
Nelm = length(x);
if sum(dy)>0 && isfinite(sum(dy))
    chi2null = sum(((y-nanwmean(y,dy.^-2))./dy).^2);
else
    chi2null = sum(((y-nanwmean(y,ones(size(dy)).^-2))./ones(size(dy))).^2);
end
Fstatistic = ((chi2null-chi2)/p2_p1)/(chi2/(Nelm-p2)); 
% Ftest = 1-fcdf(Fstatistic,p2_p1,(Nelm-p2));
Ftest = fcdf(Fstatistic,p2_p1,(Nelm-p2),'upper'); % to avoid roundoff err.
h1 = NaN;
% Plot
if saveflag || showflag
    h1 = figure('units','normalized','outerposition',[0.4 0.25 0.5 0.7]);
    if ~showflag
        set(gcf,'Visible', 'off'); 
    end
    hold off
    errorbar(x,y,dy,'ob');
    hold on;
    grid on;
    plot(x,b+a*x,'-k');
    if nargin>10
        herrorbar(x,y,dx,'ob');
    end
    set(gca,'Fontsize',16);
    xlabel(xlabel_str);
    ylabel(ylabel_str);
    title(title_str);
    % legend(['data: ' num2str(length(x)) ' points'],['line fit: a = ' num2str(a) '+/- ' num2str(da)],['\chi^2_{DOF} = ' num2str(chi2/(length(x)-2))]);
    legend(['F_{test} = ' num2str(Ftest) ' (' num2str(length(x)) ' points)'],['line fit: a = ' num2str(a) '+/- ' num2str(da)]);
    if saveflag
        saveas(gcf,[filename save_fmt]);
    end
end
if showflag
    fprintf('a = %f +/- %f \n',a,da);
    fprintf('b = %f +/- %f \n',b,db);
    fprintf('rms: %f chi2: %f\n',rms,chi2);
end