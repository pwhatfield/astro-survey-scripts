% Old script to define stellar locus

star_restriction=0*ID_galaxies;

K_FLUX_galaxies=mass_data(:,11);
J_FLUX_galaxies=mass_data(:,9);
g_FLUX_galaxies=mass_data(:,3);
i_FLUX_galaxies=mass_data(:,5);

K_MAG_galaxies=K_FLUX_galaxies*(10^23);
K_MAG_galaxies=(-1)*2.5*log10(K_MAG_galaxies)+8.9; % Convert flux to mag

J_MAG_galaxies=J_FLUX_galaxies*(10^23);
J_MAG_galaxies=(-1)*2.5*log10(J_MAG_galaxies)+8.9; % Convert flux to mag

g_MAG_galaxies=g_FLUX_galaxies*(10^23);
g_MAG_galaxies=(-1)*2.5*log10(g_MAG_galaxies)+8.9; % Convert flux to mag

i_MAG_galaxies=i_FLUX_galaxies*(10^23);
i_MAG_galaxies=(-1)*2.5*log10(i_MAG_galaxies)+8.9; % Convert flux to mag

JK_colour=J_MAG_galaxies-K_MAG_galaxies;
gi_colour=g_MAG_galaxies-i_MAG_galaxies;

for i=1:n_original;
    
    if gi_colour(i)<0.4 && JK_colour(i)<(0.12+(-0.58));
        star_restriction(i)=1;
    elseif gi_colour(i)<1.9 &&  gi_colour(i)>0.4 && JK_colour(i)<(0.12+(-0.88+0.82*gi_colour(i) - 0.21*gi_colour(i)*gi_colour(i)));
        star_restriction(i)=1;
    elseif gi_colour(i)>1.9 && JK_colour(i)<(0.12+(-0.08));
        star_restriction(i)=1;
    else
    end
        
end
