% Regions 

%% Finding gamma and lambda for acceptable region
tic
x0_1 = 0.79;
xf_1 = 1.00;
y0_1  = 0.5;
yf_1  = 0.63;
m1 = (yf_1-y0_1)/(xf_1-x0_1);

xl_1 = 0:0.0005:1;
yl_1 = 0:0.0005:1;
A1 = [];

for x1 = xl_1
    for y1 = yl_1
        y_upper = -m1*x1+yf_1;
        ynew_1 = y_upper-y1;
        xnew_1 = x1+x0_1;
        if (xnew_1 >= x0_1 && xnew_1 <= xf_1) && (ynew_1 >= y0_1 && ynew_1 <= yf_1)
            A1 = [A1; [xnew_1, ynew_1]]; %#ok<AGROW>
            Re = unique(A1,'rows');
        end
    end
end

lambdas  = Re(:,1);
lambfind = find(lambdas==0.805);
Re805    = Re(lambfind,:);

toc

