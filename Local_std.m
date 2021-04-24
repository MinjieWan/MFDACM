function [std_in,std_out ] = Local_std( Std,K,H_u,u )
H_n=1-H_u;

std_x=Std.*H_u;
std1=conv2(std_x,K,'same');
std2=conv2(H_u,K,'same');
std_in=std1./std2;

std_y=Std.*H_n;
std3=conv2(std_y,K,'same');
std4=conv2(H_n,K,'same');
std_out=std3./std4;

end

