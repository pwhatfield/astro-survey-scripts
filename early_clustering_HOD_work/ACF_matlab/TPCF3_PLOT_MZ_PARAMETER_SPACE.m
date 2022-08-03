% Old script to plot redshift stellar mass parameter space

clear all
clf
close(gcf)

TPCF2_PDF_CALC

plot(Z_galaxies,MASS_galaxies,'.','MarkerSize',9)

axis([0 4 4 12])
z1=0;
z2=0.75;
z3=1;
z4=1.5;
z5=3;

m1=8.5;
m2=9;
m3=9.5;
m4=10.25;
m5=11;

line([z1,z1],[m1,m5],'Color','r','Linewidth',3)
line([z2,z2],[m1,m5],'Color','r','Linewidth',3)
line([z3,z3],[m2,m5],'Color','r','Linewidth',3)
line([z4,z4],[m3,m5],'Color','r','Linewidth',3)
line([z5,z5],[m3,m5],'Color','r','Linewidth',3)

line([z1,z2],[m1,m1],'Color','r','Linewidth',3)
line([z1,z3],[m2,m2],'Color','r','Linewidth',3)
line([z1,z5],[m3,m3],'Color','r','Linewidth',3)
line([z1,z5],[m4,m4],'Color','r','Linewidth',3)
line([z1,z5],[m5,m5],'Color','r','Linewidth',3)

xlabel('z - redshift','FontSize', 15)
ylabel('Log(M_{*}/M_{sun})','FontSize', 15)

box on


version=3;

data=[Z_galaxies,MASS_galaxies,MAG_galaxies];


save(strcat('.../TPCF3_python_data_MZ_space_v',num2str(version),'.dat'),'data','-ascii')