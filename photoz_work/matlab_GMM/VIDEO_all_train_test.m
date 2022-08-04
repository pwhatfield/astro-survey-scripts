% All GPz
clear all

% Script that implements something close to 'GMM-All' described in https://arxiv.org/abs/2009.01952

%Number of galaxies
N_COSMOS=995049;
N_XMM=1674689;


filters=10;

% Load data
dataPath = '.../X_train.txt';
X_train = csvread(dataPath);
dataPath = '.../X_test.txt';
X_test = csvread(dataPath);
dataPath = '.../X_train_error.txt';
X_train_error = csvread(dataPath);
dataPath = '.../X_test_error.txt';
X_test_error = csvread(dataPath);
dataPath = '.../Y_train.txt';
Y_train = csvread(dataPath);
dataPath = '.../Y_test.txt';
Y_test = csvread(dataPath);
dataPath = '.../Y_spec.txt';
Y_spec = csvread(dataPath);
dataPath = '.../Y_normal.txt';
Y_normal = csvread(dataPath);
dataPath = '.../sigma_normal.txt';
sigma_normal = csvread(dataPath);

% Number of times to repeat after resampling the data based on error bars
% on photometry
number_resample=100;


%%%%%%%%%%%%%% SETUP %%%%%%%%%%%%%% 


%%%% GPzV2


rng(1); % fix random seed
addpath ../GPz/ % path to GPz
addpath(genpath('../minFunc_2012/'))       % path to minfunc
 
%%%%%%%%%%%%%% Model options %%%%%%%%%%%%%%%%
 
m = 500 %160;   %m=100 works for original                             % number of basis functions to use [required]
 
method = 'VC';                          % select a method, options = GL, VL, GD, VD, GC and VC [required]
 
joint = true;                           % jointly learn a prior linear mean-function [default=true]
 
heteroscedastic = true;                 % learn a heteroscedastic noise process, set to false if only interested in point estimates [default=true]
 
normalize = true;                       % pre-process the input by subtracting the means and dividing by the standard deviations [default=true]

maxIter = 500;                          % maximum number of iterations [default=200]
maxAttempts = 50;                       % maximum iterations to attempt if there is no progress on the validation set [default=infinity]
 
 
trainSplit = 0.5;                       % percentage of data to use for training
validSplit = 0.5;                       % percentage of data to use for validation
testSplit  = 0.0;                       % percentage of data to use for testing

inputNoise = false;                      % false = use mag errors as additional inputs, true = use mag errors as additional input noise

csl_method = 'normal';                  % cost-sensitive learning option: [default='normal']
                                        %       'balanced':     to weigh
                                        %       rare samples more heavily during training
                                        %       'normalized':   assigns an error cost for each sample = 1/(z+1)
                                        %       'normal':       no weights assigned, all samples are equally important
 
binWidth = 0.1;                         % the width of the bin for 'balanced' cost-sensitive learning [default=range(output)/100]
%%%%%%%%%%%%%% Start of script %%%%%%%%%%%%%%%% 

CSL_yn=0;

font_size=18;
line_width=4;
 
outPath = [];                           % if set to a path, the output will be saved to a csv file.
      


% read data from file
%load_data_intro

n=length(Y_train); % Number of training spectroscopic sources


% select training, validation and testing sets from the data
[training,validation,testing] = sample(n,trainSplit,validSplit,testSplit); 
 
% get the weights for cost-sensitive learning

%Normal
omega = getOmega(Y_train,csl_method,binWidth); % Original omega values for

% % GCSL
dataPath = '.../omega_GCSL.txt';  % Weights for the galaxies                                   
omega_weights = csvread(dataPath);

% Load best mixtures
dataPath = '.../best_mixture_train.txt';  % Load data that describes what mixtures the galaxies have been divided into by the GMM                                   
train_best_mixture = csvread(dataPath);

dataPath = '.../best_mixture_test.txt';  % Load data that describes what mixtures the galaxies have been divided into by the GMM                                                                      
test_best_mixture = csvread(dataPath);


% initialize the model
number_mixtures=20;

mu_individual_resample=zeros(1,length(Y_test));
sigma_individual_resample=zeros(1,length(Y_test));
nu_individual_resample=zeros(1,length(Y_test));
beta_i_individual_resample=zeros(1,length(Y_test));
gamma_individual_resample=zeros(1,length(Y_test));

mu_storage=zeros(number_resample,length(Y_test));

% Set up empty arrays
mu_all=zeros(1,length(Y_test));
sigma_all=zeros(1,length(Y_test));
nu_all=zeros(1,length(Y_test));
beta_i_all=zeros(1,length(Y_test));
gamma_all=zeros(1,length(Y_test));

mean_all=zeros(1,length(Y_test));

data_size=size(X_train);
num_train_galaxies=data_size(1);

data_size=size(X_test);
num_test_galaxies=data_size(1);

% Set up redshift arrays
redshift_values=-1:0.05:9;
redshift_distribution_array=zeros(num_test_galaxies,length(redshift_values));

redshift_ilbert=[0:0.04:6,6.1:0.1:10];
redshift_distribution_ilbert_array=zeros(num_test_galaxies,length(redshift_ilbert));


for q=1:number_resample; % Cycle over resamplings
    
    % Do the resampling
    for j=1:num_train_galaxies;
        for k=1:filters;
            X_train_resampled(j,k)=X_train(j,k)+normrnd(0,X_train_error(j,k));
        end
    end
    
    for j=1:num_test_galaxies;
        for k=1:filters;
            X_test_resampled(j,k)=X_test(j,k)+normrnd(0,X_test_error(j,k));
        end
    end
    


    for i=1:number_mixtures; % Cycle over the mixtures



        I = find(train_best_mixture == i);
        J = find(test_best_mixture == i);
        
        % Determine number of basis functions assigned to the mixture
        m_mixture=ceil(m*length(J)/length(Y_test));
        m_mixture=max(m_mixture,30);

        if length(I)>20;

            [training,validation,testing] = sample(length(I),trainSplit,validSplit,testSplit); % Divide the mixture into training and validation etc.

            
            % Weighted validation, weights determined by training and
            % testing colour-mag distribution
            omega_divide=omega_weights(I);
            for j=1:length(training);

                random_number=rand;

                if random_number>1/(1+omega_divide(j));
                    training(j)=0;
                    validation(j)=1;
                else
                    training(j)=1;
                    validation(j)=0;
                end

            end
            
            num_train=sum(training);
            num_validation=sum(validation);



            % Do the actually training GPz
            model = init(X_train_resampled(I,:),Y_train(I),method,m_mixture,'omega',omega(I),'training',training,'heteroscedastic',heteroscedastic,'joint',joint,'normalize',normalize);
            % train the model
            model = train(model,X_train_resampled(I,:),Y_train(I),'omega',omega(I),'training',training,'validation',validation,'maxIter',maxIter,'maxAttempts',maxAttempts); 

            % use the model to generate predictions for the test set
            [mu_individual_resample(J),sigma_individual_resample(J),nu_individual_resample(J),beta_i_individual_resample(J),gamma_individual_resample(J)] = predict(X_test_resampled(J,:),model);

        elseif length(I)>0;

            for k=1:length(J);
                % Treatment for if a mixture is too small to run GPz on
                % (normally when there are just a tiny number of sources
                % far from anything else
                mu_individual_resample(J(k))=mean(Y_train(I));
                sigma_individual_resample(J(k))=std(Y_train(I))^2;

            end




        else

        end

    end
    
    % Store the means from the individual resamples
    mu_storage(q,:)=mu_individual_resample(:);
    
    
    
    % Build up the non-Gaussian distributions from the multiple resamples
    for j=1:num_test_galaxies;
        redshift_distribution_array(j,:)=redshift_distribution_array(j,:)+normpdf(redshift_values,mu_individual_resample(j),sigma_individual_resample(j)^0.5);
    end
    
    for j=1:num_test_galaxies;
        redshift_distribution_ilbert_array(j,:)=redshift_distribution_ilbert_array(j,:)+normpdf(redshift_ilbert,mu_individual_resample(j),sigma_individual_resample(j)^0.5)/number_resample;
    end
    


end

'A'


for i=1:num_test_galaxies;
    
    [A,I]=max(redshift_distribution_array(i,:));
    mu_all(i)=redshift_values(I);
    
    mean_all(i)=sum((redshift_distribution_array(i,:).*redshift_values)/sum(redshift_distribution_array(i,:)));
    
    mu_all(i)=mean_all(i);
    
    sigma_all(i)=sum((redshift_distribution_array(i,:).*redshift_values.*redshift_values)/sum(redshift_distribution_array(i,:)))-(mean_all(i)^2);
end

'B'



mu_all=min(mu_all,9);
mu_all=mu_all';
sigma_all=sigma_all';

nan_indicies=isnan(sigma_all);
sigma_all(nan_indicies)=5;


for i=1:num_test_galaxies;
    
    if nan_indicies(i)==1;
   
        redshift_distribution_ilbert_array(i,:)=normpdf(redshift_ilbert,mu_all(i),sigma_all(i)^0.5);
    else
    end
    
end




'nan indicies sum'
sum(nan_indicies)


 
% mu     = the best point estimate
% nu     = variance due to data density
% beta_i = variance due to output noise
% gamma  = variance due to input noise
% sigma  = nu+beta_i+gamma
 
% compute any metrics here, e.g. RMSE
 % reduce the sample for efficient plotting
[x_all,y_all,color_all,counts_all]=reduce(Y_test,mu_all,sigma_all,200);
 
figure;scatter(x_all,y_all,12,log(color_all),'s','filled');title('Uncertainty');xlabel('Spectroscopic Redshift');ylabel('Photometric Redshift');colormap jet;line([0,5],[0,5],'Linestyle','--');
figure;scatter(x_all,y_all,12,log(counts_all),'s','filled');title('Density');xlabel('Spectroscopic Redshift');ylabel('Photometric Redshift');colormap jet;line([0,5],[0,5],'Linestyle','--');



Y_all=mu_all;
Y_template=Y_test;


'C'

null_Y_spec=Y_spec==-1000;
Y_spec(null_Y_spec)=-99;

% Array of all the relevant redshifts and uncertainties
redshifts=[Y_spec,Y_template,Y_normal,sigma_normal.^0.5,Y_all,sigma_all.^0.5];

redshifts_COSMOS=redshifts(1:N_COSMOS,:);
redshifts_XMM=redshifts(1+N_COSMOS:N_COSMOS+N_XMM,:);

redshifts_distribution_COSMOS=redshift_distribution_array(1:N_COSMOS,:);
redshifts_distribution_XMM=redshift_distribution_array(1+N_COSMOS:N_COSMOS+N_XMM,:);

redshifts_distribution_COSMOS_ilbert=redshift_distribution_ilbert_array(1:N_COSMOS,:);
redshifts_distribution_XMM_ilbert=redshift_distribution_ilbert_array(1+N_COSMOS:N_COSMOS+N_XMM,:);

'D'

csvwrite('.../GPz_prediction_data_COSMOS_all_MASTER.txt',redshifts_COSMOS)
csvwrite('.../GPz_prediction_data_XMM_all_MASTER.txt',redshifts_XMM)

'E'

csvwrite('.../GPz_prediction_data_ilbert_pdfs_COSMOS_all_MASTER.txt',redshifts_distribution_COSMOS_ilbert)
csvwrite('.../GPz_prediction_data_ilbert_pdfs_XMM_all_MASTER.txt',redshifts_distribution_XMM_ilbert)

'job all sorted'

%%%%%%%%%%%%%% Display Results %%%%%%%%%%%%%%%% 
 
% compute metrics 
 
%root mean squared error, i.e. sqrt(mean(errors^2))
rmse_all = sqrt(metrics(Y_test,mu_all,sigma_all,@(Y_test,mu_all,beta_i_all) (Y_test-mu_all).^2)); 

% mean log likelihood mean(-0.5*errors^2/sigma -0.5*log(sigma)-0.5*log(2*pi))
mll_all = metrics(Y_test,mu_all,sigma_all,@(y,mu_all,sigma_all) -0.5*(y-mu_all).^2./sigma_all - 0.5*log(sigma_all)-0.5*log(2*pi));
 
% fraction of data where |z_spec-z_phot|/(1+z_spec)<0.15
fr15_all = metrics(Y_test,mu_all,sigma_all,@(y,mu_all,sigma_all) 100*(abs(y-mu_all)./(y+1)<0.15));
 
% fraction of data where |z_spec-z_phot|/(1+z_spec)<0.05
fr05_all = metrics(Y_test,mu_all,sigma_all,@(y,mu_all,sigma_all) 100*(abs(y-mu_all)./(y+1)<0.05));
 
% bias, i.e. mean(errors)
bias_all = metrics(Y_test,mu_all,sigma_all,@(y,mu_all,sigma_all) y-mu_all);
 
% print metrics for the entire data
fprintf('RMSE\t\tMLL\t\tFR15\t\tFR05\t\tBIAS\n')
fprintf('%f\t%f\t%f\t%f\t%f\n',rmse_all(end),mll_all(end),fr15_all(end),fr05_all(end),bias_all(end))
%  

