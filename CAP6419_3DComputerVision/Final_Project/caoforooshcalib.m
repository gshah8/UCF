function K = caoforooshcalib(imname)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

syms w11 vyx vyy real

im=imread(imname);
if im(size(im,3)>1)
    im=rgb2gray(im);
end;

% Display image and get the coordinates of the follwing points in this
% order: left eye corner , right eye corner , left corner of lips, right
% corner of lips
figure; imshow(im); axis image;
[X,Y] = myginput(4);
hold on;
plot(X,Y,'r+')
hold off;

% Find the x-vanishing point
vx=cross(cross([X(1);Y(1);1],[X(2);Y(2);1]),cross([X(3);Y(3);1],[X(4);Y(4);1]));
vx=vx/vx(3);

% Find points along the direction of y-vanishing point: m1, m2, m3
m=cross(cross([X(1);Y(1);1],[X(3);Y(3);1]),cross([X(2);Y(2);1],[X(4);Y(4);1]));
m=m/m(3);
m2=cross(cross([X(1);Y(1);1],[X(4);Y(4);1]),cross([X(2);Y(2);1],[X(3);Y(3);1]));
m2=m2/m2(3);
ly=cross(m2,m);
m1=cross(cross([X(1);Y(1);1],[X(2);Y(2);1]),ly);
m1=m1/m1(3);
m3=cross(cross([X(3);Y(3);1],[X(4);Y(4);1]),ly);
m3=m3/m3(3);

% Golden ratio: assume d(m1,m3)/d(m1,m2) is equal to golden ratio
gr=(1+sqrt(5))/2;

% Use the golden ratio and the invariance of cross-ratio to find the
% y-vanishing point
eq1=((m1(2)-m3(2))/(m1(2)-m2(2)))/((vyy-m3(2))/(vyy-m2(2)))-gr;
vy(2)=double(vpa(simplify(solve(eq1,vyy))));
eq2= ly'*[vyx;vy(2);1];
vy(1)=double(vpa(simplify(solve(eq2,vyx))));
vy(3)=1;
vy=vy';

% Use orthogonality of vx and vy to find the calibration matrix K
omega=[w11, 0, 0; 0, w11, 0; 0, 0, 1];
eq=vx'*omega*vy;
w=double(vpa(simplify(solve(eq,w11))));
omega=[abs(w), 0, 0; 0, abs(w), 0; 0, 0, 1];
[K,p]=chol(inv(omega));
K=K/K(3,3);
end

