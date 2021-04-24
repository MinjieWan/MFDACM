function [gra_in,gra_out ] = Local_gradient( Img_delta,K,H_u,u )
H_n=1-H_u;

g_x=Img_delta.*H_u;
g1=conv2(g_x,K,'same');
g2=conv2(H_u,K,'same');
gra_in=g1./g2;

g_y=Img_delta.*H_n;
g3=conv2(g_y,K,'same');
g4=conv2(H_n,K,'same');
gra_out=g3./g4;

end

