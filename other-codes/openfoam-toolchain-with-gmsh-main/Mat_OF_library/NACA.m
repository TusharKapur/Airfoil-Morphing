function [x_out, y_out] = NACA(thickness, maximumChamber, positionOfMaximumChamber, gurneyLength, gurneyThickness, gurneyAdditionalAngle, isClosedGurney)

t  = thickness;
m  = maximumChamber;
p  = positionOfMaximumChamber;
gl = gurneyLength;
gt = gurneyThickness;

% t=0.2; % thickness
% m=0.09; % max chamber
% p=0.4; % position of max chamber

a0 =  0.2969;
a1 = -0.1260;
a2 = -0.3516;
a3 =  0.2843;

isFiniteTE = 0;

if isFiniteTE
    a4 = -0.1015;
else
    a4 = -0.1036;
end

beta = linspace(0,pi,100)';
x = (0.5*(1-cos(beta))); % Half cosine based spacing

yt = (t/0.2)*(a0*sqrt(x)+a1*x+a2*x.^2+a3*x.^3+a4*x.^4);

xc1 = x(find(x<=p));
xc2 = x(find(x>p));
xc  = [xc1 ; xc2];

if p==0
    xu=x;
    yu=yt;

    xl=x;
    yl=-yt;
    
    zc=zeros(size(xc));
else
    yc1=(m/p^2)*(2*p*xc1-xc1.^2);
    yc2=(m/(1-p)^2)*((1-2*p)+2*p*xc2-xc2.^2);
    zc=[yc1 ; yc2];

    dyc1_dx=(m/p^2)*(2*p-2*xc1);
    dyc2_dx=(m/(1-p)^2)*(2*p-2*xc2);
    dyc_dx=[dyc1_dx ; dyc2_dx];
    theta=atan(dyc_dx);

    xu=x-yt.*sin(theta);
    yu=zc+yt.*cos(theta);

    xl=x+yt.*sin(theta);
    yl=zc-yt.*cos(theta);
end

% gurney

if gurneyLength > 0 && gurneyThickness > 0

    R = @(phi) [cos(phi) -sin(phi); sin(phi) cos(phi)];
    
    xloriginal = xl;
    yloriginal = yl;

    lastPoint = [xl(end); yl(end)];
    secondLastPoint = [xl(end-1); yl(end-1)];
    gurneyPoint = lastPoint + R(pi/2 - gurneyAdditionalAngle) * (secondLastPoint - lastPoint);
    gurneyPoint = (gurneyPoint - lastPoint) / norm(gurneyPoint - lastPoint) * gl + lastPoint;

    xl(end + 1) = gurneyPoint(1);
    yl(end + 1) = gurneyPoint(2);

    lastPoint = [xl(end); yl(end)];
    secondLastPoint = [xl(end-1); yl(end-1)];
    gurneyPoint = lastPoint + R(pi/2) * (secondLastPoint - lastPoint);
    gurneyPoint = (gurneyPoint - lastPoint) / norm(gurneyPoint - lastPoint) * gt + lastPoint;

    xl(end + 1) = gurneyPoint(1);
    yl(end + 1) = gurneyPoint(2);

    lastPoint = [xl(end); yl(end)];
    secondLastPoint = [xl(end-1); yl(end-1)];
    gurneyPoint_ = lastPoint + R(pi/2) * (secondLastPoint - lastPoint);
    gurneyPoint_ = (gurneyPoint_ - lastPoint) / norm(gurneyPoint_ - lastPoint) * gl + lastPoint;


    % find point where we hit bottom

    while 1
        p1 = [xloriginal(end); yloriginal(end)];
        p2 = [xloriginal(end-1); yloriginal(end-1)];
        p3 = gurneyPoint;
        p4 = gurneyPoint_;

    %     figure(1); clf;
    %     plot(p1(1), p1(2), 'ko'); hold on;
    %     plot(p2(1), p2(2), 'ko'); hold on;
    %     plot(p3(1), p3(2), 'ko'); hold on;
    %     plot(p4(1), p4(2), 'ko'); hold on;

        t1 = 0 : 0.1 : 1;
        t2 = 0 : 0.1 : 1;

        g1 = p1 + (p2 - p1) * t1;
        g2 = p3 + (p4 - p3) * t2;

        %plot(g1(1, :), g1(2, :), 'r-'); hold on;
        %plot(g2(1, :), g2(2, :), 'g-'); hold on;

        x1 = p1(1); y1 = p1(2);
        x2 = p2(1); y2 = p2(2);
        x3 = p3(1); y3 = p3(2);
        x4 = p4(1); y4 = p4(2);

        t = [-(x2 - x1), (x4 - x3); -(y2 - y1), (y4 - y3);] \ [x1 - x3; y1 - y3];

        xloriginal = xloriginal(1 : end-1);
        yloriginal = yloriginal(1 : end-1);

        if 0 <= t(1) && t(1) <= 1
            pc = p1 + (p2 - p1) * t(1);

            %plot(pc(1), pc(2), 'k*'); hold on;

            xloriginal(end + 1) = pc(1);
            yloriginal(end + 1) = pc(2);

            break;
        end
    end

    % add gurney flap points to lower

    xloriginal = [xloriginal; flipud(xl(end - 3 : end))];
    yloriginal = [yloriginal; flipud(yl(end - 3 : end))];

    if ~isClosedGurney
        xloriginal = xloriginal(1 : end - 1);
        yloriginal = yloriginal(1 : end - 1);
    end
    
    xl = xloriginal;
    yl = yloriginal;


%     x_ = xl;
%     y_ = yl;
%     return;
end





x_out = [flipud(xu) ; xl(2:end)];
y_out = [flipud(yu) ; yl(2:end)];


