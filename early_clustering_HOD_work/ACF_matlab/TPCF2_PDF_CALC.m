% Old script to find weights and how much of a redshift pdf was in a given
% redshift range

clear all % Start anew!

tic % Start timer


%%%%%% Loading relevant data

%%%%%%%%%%%%%%%%%%%%

% Load the mass and other data from CFHTLS catalogue
mass_data=load('.../VIDEO_mass.txt');
'mass loaded'

% Open the z pdf from CFHTLS catalogue
z_distribs=load('.../VIDEO_z_distribs.txt');
% load('z_distribs_convolve');
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
    Z_galaxies(i)=zvalueindex/50; %%%%%%%%%% DO THIS PROPERLY
    
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

load('.../TPCF2_initial_data')


load('.../dummy_mask_restriction')
load('.../RA_dummies_orig')
load('.../DEC_dummies_orig')

X=find(dummy_mask_restriction);
RA_dummies_orig(X)=[];
DEC_dummies_orig(X)=[];




% 
% 
% 
% mass_range_upper_1=20;
% mass_range_lower_1=1;
% z_range_upper_1=4;
% z_range_lower_1=0;
% galaxy_type_1=0;
% 
% mass_range_upper_2=20;
% mass_range_lower_2=1;
% z_range_upper_2=4;
% z_range_lower_2=0;
% galaxy_type_2=0;
% 






% load('RR_TPCF2')
load('RR_mask_full_TPCF2') % 300000 dummies run, including mask


for i=1:n_original
   z_distribs(i,2:302)=z_distribs(i,2:302)/sum(z_distribs(i,2:302)); 
    
end


% Calculate peak z

% delta_yn=1;

Z_galaxies=0*ID_galaxies(:);

delta_yn=0;

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
number_bootstrap=10;
n_r=100000;

% 0=all 1=elliptical 2=irregular 3=starburst 4=irregular and starburst

% % First data set
% mass_range_upper_1=11;
% mass_range_lower_1=10;
% z_range_upper_1=4;
% z_range_lower_1=0;
% galaxy_type_1=1;
% 
% z_central_1=(z_range_upper_1+z_range_lower_1)/2;
% mass_central_1=(mass_range_upper_1+mass_range_lower_1)/2;
% 
% % Second data set
% mass_range_upper_2=10;
% mass_range_lower_2=0;
% z_range_upper_2=4;
% z_range_lower_2=0;
% galaxy_type_2=0;
% 
% z_central_2=(z_range_upper_2+z_range_lower_2)/2;
% mass_central_2=(mass_range_upper_2+mass_range_lower_2)/2;

% z_divisions=[0, 0.5, 1, 1.75, 4];







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

%%%%%%% CHECK MASK


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


% bar(binranges,bincounts,'histc')
binranges=[0,0.5,0.9,1.000000000001];
binranges_before=[0:0.01:1];
smooth=0.01;

non_param_bincounts_1=0*(0:0.01:1);
non_param_bincounts_2=0*non_param_bincounts_1;
non_param_bincounts_3=0*non_param_bincounts_1;
non_param_bincounts_4=0*non_param_bincounts_1;


z_range_upper=0.5;
z_range_lower=0;
    
z_upper_index=floor((300*z_range_upper/6))+1;
z_lower_index=floor((300*z_range_lower/6))+2;

for i=1:n_post;
    PROB_galaxies(i)=sum(z_distribs(i,z_lower_index:z_upper_index))/sum(z_distribs(i,2:302)); % Calculate probability for z range
end

bincounts_1 = histc(PROB_galaxies,binranges);
bincounts_before_1=histc(PROB_galaxies,binranges_before);

for i=1:length(non_param_bincounts_1);
   
    non_param_bincounts_1=non_param_bincounts_1+bincounts_before_1(i)*normpdf(binranges_before,binranges_before(i),smooth);
    
end




z_range_upper=1;
z_range_lower=0.5;
    
z_upper_index=floor((300*z_range_upper/6))+1;
z_lower_index=floor((300*z_range_lower/6))+2;

for i=1:n_post;
    PROB_galaxies(i)=sum(z_distribs(i,z_lower_index:z_upper_index))/sum(z_distribs(i,2:302)); % Calculate probability for z range
end

bincounts_2 = histc(PROB_galaxies,binranges);
bincounts_before_2=histc(PROB_galaxies,binranges_before);

for i=1:length(non_param_bincounts_2);
   
    non_param_bincounts_2=non_param_bincounts_2+bincounts_before_2(i)*normpdf(binranges_before,binranges_before(i),smooth);
    
end


z_range_upper=1.75;
z_range_lower=1;
    
z_upper_index=floor((300*z_range_upper/6))+1;
z_lower_index=floor((300*z_range_lower/6))+2;

for i=1:n_post;
    PROB_galaxies(i)=sum(z_distribs(i,z_lower_index:z_upper_index))/sum(z_distribs(i,2:302)); % Calculate probability for z range
end

bincounts_3 = histc(PROB_galaxies,binranges);
bincounts_before_3=histc(PROB_galaxies,binranges_before);

for i=1:length(non_param_bincounts_3);
   
    non_param_bincounts_3=non_param_bincounts_3+bincounts_before_3(i)*normpdf(binranges_before,binranges_before(i),smooth);
    
end


z_range_upper=4;
z_range_lower=1.75;
    
z_upper_index=floor((300*z_range_upper/6))+1;
z_lower_index=floor((300*z_range_lower/6))+2;

for i=1:n_post;
    PROB_galaxies(i)=sum(z_distribs(i,z_lower_index:z_upper_index))/sum(z_distribs(i,2:302)); % Calculate probability for z range
end

bincounts_4 = histc(PROB_galaxies,binranges);

bincounts_4 = histc(PROB_galaxies,binranges);
bincounts_before_4=histc(PROB_galaxies,binranges_before);

for i=1:length(non_param_bincounts_4);
   
    non_param_bincounts_4=non_param_bincounts_4+bincounts_before_4(i)*normpdf(binranges_before,binranges_before(i),smooth);
    
end



