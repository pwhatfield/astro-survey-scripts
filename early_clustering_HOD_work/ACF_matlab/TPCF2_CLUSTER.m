% Very old code to calculate an ACF for a set of galaxies

%%%%% Set up data to be filled
number_D1D2=zeros(1,q); % Counts number of galaxy pairs in each separation bin
number_D1R=zeros(1,q); % Counts number of dummy pairs in each separation bin
number_D2R=zeros(1,q);
number_RR=zeros(1,q); % Counts number of both pairs in each separation bin

% continuum_values=(-3:0.01:0);
% 
% D1D2_continuum=0*(-3:0.01:0);
% D1R_continuum=0*(-3:0.01:0);
% D2R_continuum=0*(-3:0.01:0);
% RR_continuum=0*(-3:0.01:0);
% 
% continuum_parameter=0.03;

%%%%% Generate randoms

% Set up dummy placeholders
RA_dummies=zeros(1,n_r); % Generate blank dummy RA
DEC_dummies=zeros(1,n_r); % Generate blank dummy dec
PROB_dummies=ones(1,n_r); % Dummy probability placeholder0
% % 
% % % Generate random numbers for random locations
% % randomvar_1=rand(1,n_r);
% % randomvar_2=rand(1,n_r);
% % 
% % 
% % % Generate dummy locations    
% % for num=1:n_r;
% %     RA_dummies(num)=RA_min+randomvar_1(num)*(RA_max-RA_min);
% %     DEC_dummies(num)=(360/(2*pi))*asin(sin(DEC_max*(2*pi/360))+randomvar_2(num)*(sin(DEC_min*(2*pi/360))-sin(DEC_max*(2*pi/360))));
% % end
% % 

RA_dummies=RA_dummies_orig(1:n_r);
DEC_dummies=DEC_dummies_orig(1:n_r);



%%%%% Calculate D1D2

extended_separation=[-0.1, separation,10];
extended_number_D1D2=[0,number_D1D2,0];

'D1D2 count start'

for i=1:n1;
    angdistance=acos(decpart1(DEC_galaxies_1(i),DEC_galaxies_2)+decpart2(DEC_galaxies_1(i),DEC_galaxies_2).*rapart(RA_galaxies_1(i),RA_galaxies_2));
    [count,bin]=histc(angdistance,extended_separation); % Count is how many in each bin, bin is which bin each element of angdistance goes in
    bin(length(bin))=q+2; % bin will be as long as angdistance, set the final element to the final bin to make sure right number generated
    totfreq=accumarray(bin,PROB_galaxies_2); % Bin and PROB same length, Each element of totfreq is the sum of the right PROBs
    extended_number_D1D2=extended_number_D1D2+PROB_galaxies_1(i)*totfreq';
    
    
end

keep=length(extended_number_D1D2);
extended_number_D1D2(keep)=[];
extended_number_D1D2(1)=[];
number_D1D2=extended_number_D1D2;

'D1D2 count done'
toc




%%%%% Calculate D1R

extended_separation=[-0.1, separation,10];
extended_number_D1R=[0,number_D1R,0];

'D1R count start'

for i=1:n1;
    angdistance=acos(decpart1(DEC_galaxies_1(i),DEC_dummies)+decpart2(DEC_galaxies_1(i),DEC_dummies).*rapart(RA_galaxies_1(i),RA_dummies));
    [count,bin]=histc(angdistance,extended_separation); % Count is how many in each bin, bin is which bin each element of angdistance goes in
    bin(length(bin))=q+2; % bin will be as long as angdistance, set the final element to the final bin to make sure right number generated
    totfreq=accumarray(bin',PROB_dummies); % Bin and PROB same length, Each element of totfreq is the sum of the right PROBs
    extended_number_D1R=extended_number_D1R+PROB_galaxies_1(i)*totfreq';
end

keep=length(extended_number_D1R);
extended_number_D1R(keep)=[];
extended_number_D1R(1)=[];
number_D1R=extended_number_D1R;

'D1R count done'
toc

if (mass_range_upper_1==mass_range_upper_2)  && (mass_range_lower_1==mass_range_lower_2)  && (z_range_upper_1==z_range_upper_2)  && (z_range_lower_1==z_range_lower_2)  && (galaxy_type_1==galaxy_type_2) && (ssfr_range_upper_1==ssfr_range_upper_2) && (ssfr_range_lower_1==ssfr_range_lower_2);

number_D2R=number_D1R;
    
else
    

%%%%% Calculate D2R

extended_separation=[-0.1, separation,10];
extended_number_D2R=[0,number_D2R,0];

'D2R count start'

for i=1:n2;
    angdistance=acos(decpart1(DEC_galaxies_2(i),DEC_dummies)+decpart2(DEC_galaxies_2(i),DEC_dummies).*rapart(RA_galaxies_2(i),RA_dummies));
    [count,bin]=histc(angdistance,extended_separation); % Count is how many in each bin, bin is which bin each element of angdistance goes in
    bin(length(bin))=q+2; % bin will be as long as angdistance, set the final element to the final bin to make sure right number generated
    totfreq=accumarray(bin',PROB_dummies); % Bin and PROB same length, Each element of totfreq is the sum of the right PROBs
    extended_number_D2R=extended_number_D2R+PROB_galaxies_2(i)*totfreq';
end

keep=length(extended_number_D2R);
extended_number_D2R(keep)=[];
extended_number_D2R(1)=[];
number_D2R=extended_number_D2R;

'D2R count done'
toc

end






% % % % % % % % 
% % % % % % % % 
% % % % % % % % 
% % % % % % % % %%%%% Calculate RR
% % % % % % % % 
% % % % % % % % extended_separation=[-0.1, separation,10];
% % % % % % % % extended_number_RR=[0,number_RR,0];
% % % % % % % % 
% % % % % % % % 'RR count start'
% % % % % % % % 
% % % % % % % % for i=1:n_r;
% % % % % % % %     angdistance=acos(decpart1(DEC_dummies(i),DEC_dummies)+decpart2(DEC_dummies(i),DEC_dummies).*rapart(RA_dummies(i),RA_dummies));
% % % % % % % %     [count,bin]=histc(angdistance,extended_separation); % Count is how many in each bin, bin is which bin each element of angdistance goes in
% % % % % % % %     bin(length(bin))=q+2; % bin will be as long as angdistance, set the final element to the final bin to make sure right number generated
% % % % % % % %     totfreq=accumarray(bin,PROB_dummies); % Bin and PROB same length, Each element of totfreq is the sum of the right PROBs
% % % % % % % %     extended_number_RR=extended_number_RR+PROB_dummies(i)*totfreq';
% % % % % % % % end
% % % % % % % % 
% % % % % % % % keep=length(extended_number_RR);
% % % % % % % % extended_number_RR(keep)=[];
% % % % % % % % extended_number_RR(1)=[];
% % % % % % % % number_RR=extended_number_RR;
% % % % % % % % 
% % % % % % % % 'RR count done'
% % % % % % % % toc
% % % % % % % % 
% % % % % % % % 







% Calculate Two-Point

galaxies_normalise_1=sum(PROB_galaxies_1);
galaxies_normalise_2=sum(PROB_galaxies_2);
dummies_normalise=sum(PROB_dummies);

D1D2=number_D1D2/(galaxies_normalise_1*galaxies_normalise_2);
D1R=number_D1R/(galaxies_normalise_1*dummies_normalise);
D2R=number_D2R/(galaxies_normalise_2*dummies_normalise);
RR=number_RR/(dummies_normalise*(dummies_normalise-1)); 

number_RR=RR*300000*300000;


load('RR_non_param')
% Final calculation

TPCF2_non_param

% two_point=(D1D2-(D1R+D2R)+RR)./RR;
two_point=non_param_TPCF;
two_point_bin=bin_TPCF;



