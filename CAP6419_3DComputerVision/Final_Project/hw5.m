clear all;
close all;
clc;

% folder_path = '/home/nitish/Desktop/UCF Courses/Fall2021/CAP6419/HW5/spockxy'; 
% image_path = '/home/nitish/Desktop/UCF Courses/Fall2021/CAP6419/HW5/tpol.jpg';
image_path = '/home/nitish/Desktop/UCF Courses/Fall2021/CAP6419/HW5/stark.jpg';

% mymat = load(folder_path);

%Find camera projection matrices P1 and P2 and homographies
[P1,P2,H,Hr,pnts] = extcaoforooshcalib(image_path);


%Find rectified feature points pnts1
pnts1 = hnormalise(Hr*pnts);

z = [-1 0 0;
    0 1 0;
    0 0 1];

%Find matching feature points in mirror image
pnts2 = hnormalise(Hr*z*(hnormalise(H*pnts)));


%Apply linear triangulation to get 3d points
pnts1_t = transpose(pnts1);
pnts2_t = transpose(pnts2);

pnts1_t = pnts1_t(:,1:2);
pnts2_t = pnts2_t(:,1:2);

p3d = transpose(triangulate(pnts1_t,pnts2_t,transpose(P1),transpose(P2)));

%Plot the 3d points
figure;
scatter3(p3d(1,:),p3d(2,:),p3d(3,:));
axis vis3d;

