
%Very old code to initialise ACF calculation, and select a subset of the
%data

tic % Start timer

z_central_1=(z_range_upper_1+z_range_lower_1)/2;
mass_central_1=(mass_range_upper_1+mass_range_lower_1)/2;

z_central_2=(z_range_upper_2+z_range_lower_2)/2;
mass_central_2=(mass_range_upper_2+mass_range_lower_2)/2;


%%%%%% Loading relevant data

%%%%%%%%%%%%%%%%%%%% REDO THIS

% Load the mass and other data from CFHTLS catalogue
mass_data=load('.../VIDEO_mass.txt');
'mass loaded'

% Open the z pdf from CFHTLS catalogue
z_distribs=load('.../VIDEO_z_distribs.txt');
load('.../z_distribs_convolve');
% z_distribs=z_distribs_convolve;
'z pdfs loaded'

% Load data with SED data from CFHTLS catalogue
all_data=load('.../VIDEO_all.txt');
'all data loaded'

% Load mask as prepared by Boris
load('.../mask_restriction_TPCF')
'mask loaded'

ID_galaxies=mass_data(:,1); % IDs in the CFHTLS database
RA_galaxies=mass_data(:,26); % RAs in the CFHTLS database
DEC_galaxies=mass_data(:,27); % DECs in the CFHTLS database
FLUX_galaxies=mass_data(:,11); % FLUX in the CFHTLS database
SED_galaxies=all_data(:,39); % SED in CFHTLS database
MASS_galaxies=mass_data(:,131); % Mass in CFHTLS database
STAR_galaxies=all_data(:,41);
SFR_galaxies=mass_data(:,132);
SSFR_galaxies=mass_data(:,133);

KMAG_AUTO_galaxies=mass_data(:,133);
KMAG_APP2_galaxies=mass_data(:,133);

n_original=length(ID_galaxies); % Initial total number of galaxies


% Assign galaxies to elliptical (1), irregular/Sbc/Scd (2) or starburst (3)
TYPE_galaxies=0*ID_galaxies;
for i=1:n_original;
    if SED_galaxies(i)<21.5;
        TYPE_galaxies(i)=1;
    elseif (SED_galaxies(i)>21.5 && SED_galaxies(i)<58.5);
        TYPE_galaxies(i)=2;
    elseif SED_galaxies(i)>58.5;
        TYPE_galaxies(i)=3;
    else
    end
end


% Calculate Magnitude

MAG_galaxies=FLUX_galaxies*(10^23);
MAG_galaxies=(-1)*2.5*log10(MAG_galaxies)+8.9; % Convert flux to mag
%%%%%%%%%%% CHECK THIS


% Calculate peak z

Z_galaxies=0*ID_galaxies(:);

for i=1:n_original;
   
    [zvaluemaxpdf,zvalueindex]=max(z_distribs(i,:)); 
    Z_galaxies(i)=zvalueindex/50; 
    
end





% Save copies of originals
ID_galaxies_orig=ID_galaxies;
RA_galaxies_orig=RA_galaxies;
DEC_galaxies_orig=DEC_galaxies;
FLUX_galaxies_orig=FLUX_galaxies;
SED_galaxies_orig=SED_galaxies;
z_distribs_orig=z_distribs;
MASS_galaxies_orig=MASS_galaxies;
TYPE_galaxies_orig=TYPE_galaxies;
MAG_galaxies_orig=MAG_galaxies;
Z_galaxies_orig=Z_galaxies;
STAR_galaxies_orig=STAR_galaxies;
SFR_galaxies_orig=SFR_galaxies;
SSFR_galaxies_orig=SSFR_galaxies;

RA_min=36;
RA_max=37;
DEC_min=-5;
DEC_max=-4;


'Loaded initial data'
toc

%%% clear some stuff

load('TPCF2_initial_data')


load('dummy_mask_restriction')
load('RA_dummies_orig')
load('DEC_dummies_orig')

X=find(dummy_mask_restriction);
RA_dummies_orig(X)=[];
DEC_dummies_orig(X)=[];



if exist('conv_yn')==0;
    conv_yn=0;
else
end

if exist('conv_param')==0;
    conv_param=1;
else
end

'stop'

% load('RR_TPCF2')
load('RR_mask_full_TPCF2') 

z_distribs_conv=z_distribs*0;

for i=1:n_original
   z_distribs(i,2:302)=z_distribs(i,2:302)/sum(z_distribs(i,2:302));
   
   if conv_yn==1;
   
       for j=2:302;
           z_distribs_conv(i,2:302)=z_distribs_conv(i,2:302)+z_distribs(i,j)*normpdf(0:0.02:6,(j-1)*0.02,conv_param);
       end
   
   else
   end
   
    
end




Z_galaxies=0*ID_galaxies(:);

for i=1:n_original;
   
    [zvaluemaxpdf,zvalueindex]=max(z_distribs(i,2:302)); 
    Z_galaxies(i)=(zvalueindex-1)/50; %%%%%%%%%% DO THIS PROPERLY
    
    if delta_yn==1;
    
        z_distribs(i,2:302)=z_distribs(i,2:302)*0;
        z_distribs(i,zvalueindex)=1;
    
    else
    end
    
    
end
Z_galaxies_orig=Z_galaxies;



%%%%% Define parameters of correlation


% Define angular intervals for TPCF to be calculated
separation_gap=0.01;
separation_log=-3:0.01:0; % In log10 degrees, separation between bins
separation=(10.^separation_log)*((2*pi)/360); % In radians, separation between bins

q=length(separation)-1; % Number of bins

distances_log=zeros(1,q); % In log10 degrees, centres of bins
for num=1:q;
    distances_log(num)=(separation_log(num)+separation_log(num+1))/2; % Mid-log-bin in log10 degrees
end
distances=10.^distances_log;


% Number of bootstrap calculations performed


if exist('number_bootstrap')==1;
else
    number_bootstrap=1;
end



n_r=100000;

% 0=all 1=elliptical 2=irregular 3=starburst 4=irregular and starburst






%%%%%% Eliminate generic unwanted data

mag_restriction=0*ID_galaxies;
type_restriction=0*ID_galaxies;
star_restriction=0*ID_galaxies;


for i=1:n_original;
    
    if MAG_galaxies(i)>23.5 && FLUX_galaxies(i)~=0 && MAG_galaxies(i)~=NaN; % Magnitude faintness limit
        mag_restriction(i)=1; % Pick out ones too faint
    else
    end
end


TPCF2_STAR




% Combine all
all_restriction=(mag_restriction|mask_restriction')|star_restriction; % Pick out all ones to be deleted
X=find(all_restriction);

% Delete unwanted
ID_galaxies(X)=[];
RA_galaxies(X)=[];
DEC_galaxies(X)=[];
MAG_galaxies(X)=[];
FLUX_galaxies(X)=[];
SED_galaxies(X)=[];
z_distribs(X,:)=[];
MASS_galaxies(X)=[];
TYPE_galaxies(X)=[];
Z_galaxies(X)=[];
STAR_galaxies(X)=[];
SFR_galaxies(X)=[];
SSFR_galaxies(X)=[];

n_post=length(RA_galaxies);

'done'

random_delete_half=0;

if random_delete_half==1;
    random_restriction=randi([0 1],1,length(RA_galaxies));
else
    random_restriction=0*randi([0 1],1,length(RA_galaxies));
end


%%%%%% Isolate datasets

%%% Data set 1
ID_galaxies_1=ID_galaxies;
RA_galaxies_1=RA_galaxies;
DEC_galaxies_1=DEC_galaxies;
MAG_galaxies_1=MAG_galaxies;
FLUX_galaxies_1=FLUX_galaxies;
SED_galaxies_1=SED_galaxies;
z_distribs_1=z_distribs;
MASS_galaxies_1=MASS_galaxies;
TYPE_galaxies_1=TYPE_galaxies;
Z_galaxies_1=Z_galaxies;
STAR_galaxies_1=STAR_galaxies;
SFR_galaxies_1=SFR_galaxies;
SSFR_galaxies_1=SSFR_galaxies;

PROB_galaxies_1=zeros(1,n_post);

% Find indicies for z ranges
z_upper_index=floor((300*z_range_upper_1/6))+1;
z_lower_index=floor((300*z_range_lower_1/6))+2;

% Find probabilities for the galaxies in this z range

for i=1:n_post;
    PROB_galaxies_1(i)=sum(z_distribs_1(i,z_lower_index:z_upper_index))/sum(z_distribs_1(i,2:302)); % Calculate probability for z range
end




% Type restriction
type_restriction=zeros(1,n_post);

for i=1:n_post;
    
    if galaxy_type_1==0;
        type_restriction(i)=1;
    elseif galaxy_type_1==1 && TYPE_galaxies_1(i)==1;
        type_restriction(i)=1;
    elseif galaxy_type_1==2 && TYPE_galaxies_1(i)==2;
        type_restriction(i)=1;
    elseif galaxy_type_1==3 && TYPE_galaxies_1(i)==3;
        type_restriction(i)=1;
    elseif galaxy_type_1==4 && (TYPE_galaxies_1(i)==2 | TYPE_galaxies_1(i)==3);
        type_restriction(i)=1;
    else
    end
end

type_restriction=1-type_restriction;

% SSFR restriction
ssfr_restriction=zeros(1,n_post);

% ssfr_range_upper_1=100000;
% ssfr_range_lower_1=-100000;

for i=1:n_post;
    
    if ssfr_range_upper_1>SSFR_galaxies_1(i) && ssfr_range_lower_1<=SSFR_galaxies_1(i) ;
        ssfr_restriction(i)=1;
    else
    end
end

ssfr_restriction=1-ssfr_restriction;


% Mass restriction

mass_restriction=zeros(1,n_post);
MASS_ALT_galaxies_1=zeros(1,n_post); 
for i=1:n_post;
    
    MASS_ALT_galaxies_1(i)=MASS_galaxies_1(i)+2*log10(1+z_central_1)-2*log10(1+Z_galaxies_1(i));
    MASS_ALT_galaxies_1(i)=MASS_galaxies_1(i);
    
    if MASS_ALT_galaxies_1(i)<mass_range_upper_1 && MASS_ALT_galaxies_1(i)>mass_range_lower_1 && MASS_galaxies_1(i)~=99999;
        mass_restriction(i)=0;
    else
        mass_restriction(i)=1;
    end
end

% Combine all
all_restriction=(((mass_restriction|type_restriction)|ssfr_restriction)|random_restriction); % Pick out all ones to be deleted
X=find(all_restriction);

% Delete unwanted
ID_galaxies_1(X)=[];
RA_galaxies_1(X)=[];
DEC_galaxies_1(X)=[];
MAG_galaxies_1(X)=[];
FLUX_galaxies_1(X)=[];
SED_galaxies_1(X)=[];
z_distribs_1(X,:)=[];
MASS_galaxies_1(X)=[];
TYPE_galaxies_1(X)=[];
Z_galaxies_1(X)=[];
STAR_galaxies_1(X)=[];
SFR_galaxies_1(X)=[];
SSFR_galaxies_1(X)=[];
PROB_galaxies_1(X)=[];

% Number of galaxies first set
n1=length(ID_galaxies_1);



%%% Data set 2
ID_galaxies_2=ID_galaxies;
RA_galaxies_2=RA_galaxies;
DEC_galaxies_2=DEC_galaxies;
MAG_galaxies_2=MAG_galaxies;
FLUX_galaxies_2=FLUX_galaxies;
SED_galaxies_2=SED_galaxies;
z_distribs_2=z_distribs;
MASS_galaxies_2=MASS_galaxies;
TYPE_galaxies_2=TYPE_galaxies;
Z_galaxies_2=Z_galaxies;
STAR_galaxies_2=STAR_galaxies;
SFR_galaxies_2=SFR_galaxies;
SSFR_galaxies_2=SSFR_galaxies;

PROB_galaxies_2=zeros(1,n_post);

% Find indicies for z ranges
z_upper_index=floor((300*z_range_upper_2/6))+1;
z_lower_index=floor((300*z_range_lower_2/6))+2;

% Find probabilities for the galaxies in this z range

for i=1:n_post;
    PROB_galaxies_2(i)=sum(z_distribs_2(i,z_lower_index:z_upper_index))/sum(z_distribs_2(i,2:302)); % Calculate probability for z range
end

% Type restriction
type_restriction=zeros(1,n_post);

for i=1:n_post;
    
    if galaxy_type_2==0;
        type_restriction(i)=1;
    elseif galaxy_type_2==1 && TYPE_galaxies_2(i)==1;
        type_restriction(i)=1;
    elseif galaxy_type_2==2 && TYPE_galaxies_2(i)==2;
        type_restriction(i)=1;
    elseif galaxy_type_2==3 && TYPE_galaxies_2(i)==3;
        type_restriction(i)=1;
    elseif galaxy_type_2==4 && (TYPE_galaxies_2(i)==2 | TYPE_galaxies_2(i)==3);
        type_restriction(i)=1;
    else
    end
end

type_restriction=1-type_restriction;


% SSFR restriction
ssfr_restriction=zeros(1,n_post);

% ssfr_range_upper_2=100000;
% ssfr_range_lower_2=-100000;

for i=1:n_post;
    
    if ssfr_range_upper_2>SSFR_galaxies_2(i) && ssfr_range_lower_2<=SSFR_galaxies_2(i) ;
        ssfr_restriction(i)=1;
    else
    end
end

ssfr_restriction=1-ssfr_restriction;

% Mass restriction

mass_restriction=zeros(1,n_post);
MASS_ALT_galaxies_2=zeros(1,n_post); 
for i=1:n_post;
    
    MASS_ALT_galaxies_2(i)=MASS_galaxies_2(i)+2*log10(1+z_central_2)-2*log10(1+Z_galaxies_2(i));
    MASS_ALT_galaxies_2(i)=MASS_galaxies_2(i);
    
    if MASS_ALT_galaxies_2(i)<mass_range_upper_2 && MASS_ALT_galaxies_2(i)>mass_range_lower_2 && MASS_galaxies_2(i)~=99999;
        mass_restriction(i)=0;
    else
        mass_restriction(i)=1;
    end
end

% Combine all
all_restriction=(((mass_restriction|type_restriction)|ssfr_restriction)|random_restriction); % Pick out all ones to be deleted
X=find(all_restriction);

% Delete unwanted
ID_galaxies_2(X)=[];
RA_galaxies_2(X)=[];
DEC_galaxies_2(X)=[];
MAG_galaxies_2(X)=[];
FLUX_galaxies_2(X)=[];
SED_galaxies_2(X)=[];
z_distribs_2(X,:)=[];
MASS_galaxies_2(X)=[];
TYPE_galaxies_2(X)=[];
Z_galaxies_2(X)=[];
STAR_galaxies_2(X)=[];
SFR_galaxies_2(X)=[];
SSFR_galaxies_2(X)=[];
PROB_galaxies_2(X)=[];

% Number of galaxies first set
n2=length(ID_galaxies_2);

n1
n2

clear mass_data all_data

%%%%% Actually find the two-point correlation function

'Main calculation starting'

TPCF2_CLUSTER

D1D2_best=D1D2;
D1R_best=D1R;
D2R_best=D2R;


two_point_best=two_point;

two_point_best_bin=bin_TPCF;

RA_galaxies_1_store=RA_galaxies_1;
DEC_galaxies_1_store=DEC_galaxies_1;
PROB_galaxies_1_store=PROB_galaxies_1;
RA_galaxies_2_store=RA_galaxies_2;
DEC_galaxies_2_store=DEC_galaxies_2;
PROB_galaxies_2_store=PROB_galaxies_2;

two_point_bootstrap=zeros(number_bootstrap,q);

bin_test_size=16;

two_point_bootstrap_bin=zeros(number_bootstrap,bin_test_size);

for iii=1:number_bootstrap;
    
    iii
    
    if rand_weight==1;
    
    random_identities_1=discretesample(PROB_galaxies_1_store/sum(PROB_galaxies_1_store),n1);
    random_identities_2=discretesample(PROB_galaxies_2_store/sum(PROB_galaxies_2_store),n2);
    
    else
    
    random_identities_1=randi(n1,n1,1);
    random_identities_2=randi(n2,n2,1);
    
    end
    

    for ii=1:n1;
        RA_galaxies_1(ii)=RA_galaxies_1_store(random_identities_1(ii));
        DEC_galaxies_1(ii)=DEC_galaxies_1_store(random_identities_1(ii));
        PROB_galaxies_1(ii)=PROB_galaxies_1_store(random_identities_1(ii)); 
    end
    
    for ii=1:n2;
        RA_galaxies_2(ii)=RA_galaxies_2_store(random_identities_2(ii));
        DEC_galaxies_2(ii)=DEC_galaxies_2_store(random_identities_2(ii));
        PROB_galaxies_2(ii)=PROB_galaxies_2_store(random_identities_2(ii)); 
    end
    
    TPCF2_CLUSTER
    
    two_point_bootstrap(iii,:)=two_point;
    two_point_bootstrap_bin(iii,:)=bin_TPCF;
    
end

two_point_mean=zeros(1,q);
two_point_upper=zeros(1,q);
two_point_lower=zeros(1,q);

two_point_mean_bin=zeros(1,bin_test_size);
two_point_upper_bin=zeros(1,bin_test_size);
two_point_lower_bin=zeros(1,bin_test_size);

for i=1:q;
    
    two_point_mean(i)=mean(two_point_bootstrap(:,i));
    two_point_upper(i)=mean(two_point_bootstrap(:,i))+std(two_point_bootstrap(:,i));
    two_point_lower(i)=mean(two_point_bootstrap(:,i))-std(two_point_bootstrap(:,i));
    
end

for i=1:bin_test_size;
    
    two_point_mean_bin(i)=mean(two_point_bootstrap_bin(:,i));
    two_point_upper_bin(i)=mean(two_point_bootstrap_bin(:,i))+std(two_point_bootstrap_bin(:,i));
    two_point_lower_bin(i)=mean(two_point_bootstrap_bin(:,i))-std(two_point_bootstrap_bin(:,i));
    
end

