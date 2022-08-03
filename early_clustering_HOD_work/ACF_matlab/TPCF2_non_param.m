%Very old script to calculate continuous non-binned ACF

K_separation=0.01;
K_values=-100:K_separation:100;
separation_gap=0.01;
K=normpdf(K_values,0,1);

c_1=sum((K_values(:).*K_values(:)).*(K(:)))*K_separation;
c_1=1;

c_2=sum((K(:).*K(:)))*K_separation;


range_upper=-0.5;
range_lower=-3;


RR_derriv=diff(RR)/separation_gap;
RR_derriv=diff(RR_derriv)/separation_gap;
RR_derriv_partial=RR_derriv(1+round(length(distances_log)*((range_lower+3)/(max(distances_log)-min(distances_log)))):round(length(distances_log)*((range_upper+3)/(max(distances_log)-min(distances_log)))));

A=sum(RR_derriv_partial(:).*RR_derriv_partial(:))*separation_gap;

smooth=0.15; 

non_param_D1D2=0*RR;
non_param_D1R=0*RR;
non_param_D2R=0*RR;
non_param_RR=0*RR;

X=log10(distances);

smooth_D1D2=smooth;
smooth_D1R=smooth;
smooth_D2R=smooth;
smooth_RR=smooth;


for i=1:length(RR);
    
    D1D2_count=sum(number_D1D2(1+round(length(distances_log)*((range_lower+3)/(max(distances_log)-min(distances_log)))):round(length(distances_log)*((range_upper+3)/(max(distances_log)-min(distances_log))))));
    non_param_D1D2=non_param_D1D2+D1D2(i)*normpdf(X,X(i),smooth_D1D2);
    
    D1R_count=sum(number_D1R(1+round(length(distances_log)*((range_lower+3)/(max(distances_log)-min(distances_log)))):round(length(distances_log)*((range_upper+3)/(max(distances_log)-min(distances_log))))));
    non_param_D1R=non_param_D1R+D1R(i)*normpdf(X,X(i),smooth_D1R);
    
    D2R_count=sum(number_D2R(1+round(length(distances_log)*((range_lower+3)/(max(distances_log)-min(distances_log)))):round(length(distances_log)*((range_upper+3)/(max(distances_log)-min(distances_log))))));
    non_param_D2R=non_param_D2R+D2R(i)*normpdf(X,X(i),smooth_D2R);
    
    RR_count=sum(number_RR(1+round(length(distances_log)*((range_lower+3)/(max(distances_log)-min(distances_log)))):round(length(distances_log)*((range_upper+3)/(max(distances_log)-min(distances_log))))));
    non_param_RR=non_param_RR+RR(i)*normpdf(X,X(i),smooth_RR);

    
    
end

non_param_TPCF=((non_param_D1D2-(non_param_D1R+non_param_D2R)+non_param_RR)./non_param_RR);




bin_size=0.2;
orig_bin=distances_log(2)-distances_log(1);
ratio=bin_size/orig_bin;
bin_points=-3:bin_size:0;

D1D2_bin_points=0*bin_points;
D1R_bin_points=0*bin_points;
D2R_bin_points=0*bin_points;
RR_bin_points=0*bin_points;

for i=1:(length(bin_points)-1);
    
    lowr=floor(((i-1)*ratio)+1);
    highr=floor((i*ratio));
    
    D1D2_bin_points(i)=sum(D1D2(lowr:highr));
    D1R_bin_points(i)=sum(D1R(lowr:highr));
    D2R_bin_points(i)=sum(D2R(lowr:highr));
    RR_bin_points(i)=sum(RR(lowr:highr));
    
    
end



bin_TPCF=(D1D2_bin_points-(D1R_bin_points+D2R_bin_points)+RR_bin_points)./RR_bin_points;
