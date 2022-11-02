function [P1,P2,H,Hr,pnts] = extcaoforooshcalib(im)
% This matlab function implements a single-image variation of Cao-Foroosh 
% calibration method using mirror-symmetric objects (e.g. faces). 
% The paper to refernce is: 
% X. Cao, H. Foroosh, Camera Calibration Using Symmetric Objects, IEEE 
% Trans. on Image Processing, Volume 15, Issue 11, Pages 3614�3619, 2006.
%
% Input: imname, is the filename of an image in the current directory
% Once the image is loaded you need to click on four points (i.e. outer
% eye-corners and outer lip corners
%
% Output:   P1 and P2 are camera projection matrices for the original image 
%           and the mirrored image is the rotation matrix
%           H is the homology that maps any point to its mirror point
%           pnts is a matrix containg the set of clicked points
%           Hr is the rectifying homography
% 
% Copyright Hassan Foroosh, 2006

syms w11 vyx vyy real

%im=imread(imname);
% if im(size(im,3)>1)
%     im=rgb2gray(im);
% end;

% Display image and get the coordinates of the follwing points in this
% order: left eye outer-corner , right eye outer-corner, left corner of 
% lips, right corner of lips
% set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
imshow(im); axis image;
% [X,Y] = myginput(4);
[X,Y] = ginput();
hold on;
plot(X,Y,'r+')
hold off;

% Origin of image coordinate system is assumed to be the center pixel
h=round(size(im,1)/2);
w=round(size(im,2)/2);
X=X-w;
Y=Y-h;

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

% Using the invariance of cross-ratio, we know that the
% cross-ratio below is the same as the ratio of the world points
% corresponding to r=(X(1)-X(2))/(X(3)-X(4))
cr=((X(1)-X(2))/(X(3)-X(4)))/((vx(1)-X(1))/(vx(1)-X(4)));
% Comment the above line and uncomment the one below, if you want to use a
% variation of Cao-Froroosh method based on golden-ratio
% cr=(1+sqrt(5))/2;
Hr=rectH(X,Y,cr);

% Using similar triangles and the invariance of cross-ratio, we can find 
% the y-vanishing point
eq1=vpa(simplify(((m1(2)-m2(2))/(m2(2)-m3(2)))/((m1(2)-vyy)/(m3(2)-vyy))))-cr;
vy2=double(vpa(simplify(solve(eq1,vyy))));
% On the other hand, vy is on the line ly
eq2= ly'*[vyx;vy2;1];
vy1=double(vpa(simplify(solve(eq2,vyx))));
vy=[vy1;vy2;1];

% Use orthogonality of vx and vy to find the calibration matrix K
omega=[w11, 0, 0; 0, w11, 0; 0, 0, 1];
eq=vx'*omega*vy;
w=double(vpa(simplify(solve(eq,w11))));
% omega=Hr'*[abs(w), 0, 0; 0, abs(w), 0; 0, 0, 1]*Hr;
omega=[abs(w), 0, 0; 0, abs(w), 0; 0, 0, 1];
[K,p]=chol(inv(omega));
K=K/K(3,3);

% Find the rotation matrix from vanishing points using the fact that the
% first and second column of the camera matrix are the x and y vanishing
% points up to a scale factor
r1=inv(K)*vx/norm(inv(K)*vx);
r2=inv(K)*vy/norm(inv(K)*vy);
r3=cross(r1,r2)/norm(cross(r1,r2));
R=[r1,r2,r3];
% Make sure R is orthonormal
[U,D,V]=svd(R);
d=(D(1,1)+D(2,2)+D(3,3))/3;
D(1,1)=d;D(2,2)=d;D(3,3)=d;
R=U*D*V';
RR=sqrt(R*R'); 
R=R/RR(1,1);
z=[-1,0,0;0,1,0;0,0,1];

H=eye(3)-2*vx*ly'/(vx'*ly);

Rp=z*R;
Rp(:,3)=-Rp(:,3);
% Find the camera projection matrix P assuming that m1 is the world origin,
% because it is on the intersection of l_y (where v_y is) and the line 
% l_x = m1 x m2 (where v_x is):
% v_x^T(m1 x m2)=0 and v_y^T l_y=0 
% P1=[K*R,m2];
% P2=[K*Rp,z*hnormalise(H*m2)];
zp=[-1,0,0;0,1,0;0,0,-1];
P1=[K,hnormalise(Hr*m1)];
P2=[K*zp,z*hnormalise(Hr*m1)];

pnts=[X';Y';ones(size(X'))];


end

function Hr = rectH(X,Y,cr)

syms x3 x4 y3 real

pnts=[X';Y';ones(size(X'))];

eq1=(X(2)-X(1))-cr*(x4-x3);
eq2=x3-X(1)-X(2)+x4;
% S=solve({eq1,eq2},{x3,x4});
S=solve([eq1,eq2],[x3,x4]);

x1=[X(1:4)';Y(1:4)';ones(size(X(1:4)'))];
x2=[pnts(:,1),[pnts(1,2);pnts(2,1);1],[vpa(S.x3);pnts(2,3);1],[vpa(S.x4);pnts(2,3);1]];
Npts = 4;
A = zeros(3*Npts,9);
O = [0 0 0];
for n = 1:Npts
	X = x1(:,n)';
	x = x2(1,n); y = x2(2,n); w = x2(3,n);
	A(3*n-2,:) = [  O  -w*X  y*X];
	A(3*n-1,:) = [ w*X   O  -x*X];
	A(3*n  ,:) = [-y*X  x*X   O ];
end

[U,D,V] = svd(A,0); % 'Economy' decomposition for speed
    
% Extract homography
Hr = reshape(V(:,9),3,3)';

end
