function [ entro_in,entro_out ] = Local_entropy( imgn,K,H_u,u)
H_n=1-H_u;
entro_x=H_u.*imgn;
entro1=conv2(entro_x,K,'same');
entro2=conv2(H_u,K,'same');
entro_in=entro1./entro2;

entro_y=H_n.*imgn;
entro3=conv2(entro_y,K,'same');
entro4=conv2(H_n,K,'same');
entro_out=entro3./entro4;

end

