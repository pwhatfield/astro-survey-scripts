%Old script to generate python readable tables

%zbin=4;
%massbin=7;
version=1;

z_division=0.02;

z_interval=[0.5,0.75,1,1.25,1.7,2.3];


load(strcat('.../TPCF3_bias_all_zbin',num2str(zbin),'_massbin',num2str(massbin),'_v',num2str(version),'.mat'))


theta=distances_log(:);
tcpf_best=log10(two_point_best(:));
log_error=log10(two_point_upper(:))-log10(two_point_best(:));


theta=theta(1:1:length(theta));
tcpf_best=tcpf_best(1:1:length(tcpf_best));
log_error=log_error(1:1:length(log_error));

data=[theta,tcpf_best,log_error];


N=sum(z_distribs_1)/z_division;
N(1)=[];


distrib=N;
z_values=z_division:z_division:(z_division*length(distrib));

for i=1:length(distrib);
    if (z_values(i)<z_interval(zbin))||(z_values(i)>z_interval(zbin+1))
        distrib(i)=0;
    else
    end
end

N=[z_values',distrib']; %N=[z_values,distrib];

%save(strcat('.../TPCF3_python_data_redshift_zbin',num2str(zbin),'_massbin',num2str(massbin),'_v',num2str(version),'.dat'),'N','-ascii')
