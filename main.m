clc;clear all;close all;
ImgID = 3;
Img = imread([num2str(ImgID),'.bmp']);
Img = (Img(:,:,1));
Img=im2uint8(Img) ;
Img_entropy=entropyfilt(Img);
Img = double(Img);
[row,col] = size(Img);
%% initial the evolving curve
phi = ones(row,col);
[nx,ny]=size(Img);
x0 = round(ny/2);
y0 = round(nx/2);
r1 = y0 - 40;
r2 = y0 + 40;
c1 = x0 - 40;
c2 = x0 + 40;
phi(r1:r2,c1:c2) = -1;
u =  phi;
figure(1);
imshow(Img,[0,255]);
hold on;
axis off;
[c, h] = contour(u, [0 0], 'g','LineWidth',2);
title('Initial contour');
%% parameter setting
sigma =3.0;
G = fspecial('gaussian', 5, sigma);
delt = 1;
epsilon=1.0; 
n = 1;
sigma1=2.0;
K=fspecial('gaussian',round(2*sigma1)+1,1); 
[Img_std] = Stdfilt(Img,3);
[Ix, Iy] = gradient(Img);
Img_delta=sqrt( Ix.^2 + Iy.^2 );
fun1 = @(x) max(x(:));
B1 = nlfilter(Img, [5 5], fun1);
fun2 = @(x) min(x(:));
B2 = nlfilter(Img, [5 5], fun2);
range=B1-B2;

switch ImgID
    case 1
        mu=400;
        eps_thr = 50;
        lambda=0.01;
    case 2
        mu=400;
        eps_thr = 50;
        lambda=0.01;
    case 3
        mu=150;
        eps_thr = 50;
        lambda=0.5;
    case 4
        mu=400;
        eps_thr = 50;
        lambda=0.01;
end
w=1./(1+lambda*range);
%% start level set evolution
while (1)
    [ux, uy] = gradient(u);
    H_u = 0.5*(1+(2/pi)*atan(u/epsilon));
    delta_H = (1/pi)*epsilon./(epsilon^2+u.^2);
    
    % compute the global term
    c1 = sum(sum(Img.*(u<0)))/(sum(sum(u<0)));
    c2 = sum(sum(Img.*(u>=0)))/(sum(sum(u>=0)));
    I=c1.*H_u+c2.*(1-H_u);
    spf_g=Img-I;
    
    % compute the local term
    [ entro_in,entro_out ] = Local_entropy( Img_entropy,K,H_u,u);
    [ std_in,std_out ] = Local_std( Img_std,K,H_u,u );
    [ gra_in,gra_out ] = Local_gradient( Img_delta,K,H_u,u );
    entropy=entro_in.*H_u+entro_out.*(1-H_u);
    spf_en=Img_entropy-entropy;
    std=std_in.*H_u+std_out.*(1-H_u);
    spf_std=Img_std-std;
    Img_gradient=gra_in.*H_u+gra_out.*(1-H_u); 
    spf_gra=Img_delta-Img_gradient;
    spf_l=spf_en+spf_std+spf_gra;
    
    % construct the complete level set function
    spf_total = w.*spf_g+(1-w).*(spf_l);
    spf_total = spf_total/(max(abs(spf_total(:))));
    u = u + delt * (( mu * spf_total.* sqrt( ux.^2 + uy.^2 )));
    u = (u >= 0) - ( u< 0);
  
    if n >= eps_thr
        break;
    end 
    u = conv2(u, G, 'same');   
    if mod(n,1)==0 
    imagesc(Img,[0 255]); colormap(gray);hold on;axis off;
    [c, h] = contour(u, [0 0], 'r','LineWidth',1);
    iterNum = [num2str(n), 'iterations'];
    title(iterNum);
    pause(0.01);
    end
    n = n + 1;
end

%% imshow the result
u = conv2(u, G, 'same');
u = (u >= 0) - ( u< 0);
bw = u;
bw(bw<0) = 0; 

figure(2);
imagesc(Img,[0 255]); colormap(gray);
hold on;
axis off;
[c, h] = contour(u, [0 0], 'r','LineWidth',2);

figure(3);
imshow(bw);


   
