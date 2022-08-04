% Construct GCSL weights

%flux_u 1
%flux_HSC-G 2
%flux_HSC-R 3
%flux_HSC-I 4
%flux_HSC-Z 5
%flux_HSC-Y 6
    
%flux_Y 7
%flux_J 8
%flux_H 9
%flux_Ks 10

%%%%%%%%%%%%%% CONFIGURATION %%%%%%%%%%%%%%%%
rng(1);                                 % fix random seed
 
addpath ../VB_GMM/                         % path to the VB library
 
K = 20;                                % number of mixtures to use [required]
 

maxIter = 100;                          % maximum number of iterations
tol = 1e-10;                            % early stoping criteria if no significant improvement is gained
init= 'conditional';                    % how to initialize ('rand', 'conditional', 'kmeans')
cov = 'full';                           % type of covariance ('full', 'diag')
display=1;                              % display progress
 
 
trainSplit = 1.0;                       % percentage of data to use for training
validSplit = 0.0;                       % percentage of data to use for validation
testSplit  = 0.0;                       % percentage of data to use for testing

 
bins = 100;                            % number of samples to plot the CDFs and PDFs


% set training options
options.cyc=maxIter;
options.tol=tol;
options.init=init;
options.cov=cov;
options.display=display;

%%%%%%


% read data from file
%load_data_intro

% Bespoke colour spaces:

colour_space_train=0*X_train;

% Train
colour_space_train(:,1)=-2.5*((X_train(:,2)-X_train(:,4))); %g-i
colour_space_train(:,2)=-2.5*((X_train(:,8)-X_train(:,10))); %J-K
colour_space_train(:,3)=-2.5*((X_train(:,2)-X_train(:,1))); %g-u
colour_space_train(:,4)=-2.5*((X_train(:,2)-X_train(:,3))); %g-r
colour_space_train(:,5)=-2.5*((X_train(:,2)-X_train(:,5))); %g-z
colour_space_train(:,6)=-2.5*((X_train(:,2)-X_train(:,6))); %g-Y1
colour_space_train(:,7)=-2.5*((X_train(:,2))); %g
colour_space_train(:,8)=-2.5*((X_train(:,8)-X_train(:,7))); %J-Y2
colour_space_train(:,9)=-2.5*((X_train(:,8)-X_train(:,9))); %J-H
colour_space_train(:,10)=-2.5*((X_train(:,8))); %J


colour_space_test=0*X_test;

% Train
colour_space_test(:,1)=-2.5*((X_test(:,2)-X_test(:,4))); %g-i
colour_space_test(:,2)=-2.5*((X_test(:,8)-X_test(:,10))); %J-K
colour_space_test(:,3)=-2.5*((X_test(:,2)-X_test(:,1))); %g-u
colour_space_test(:,4)=-2.5*((X_test(:,2)-X_test(:,3))); %g-r
colour_space_test(:,5)=-2.5*((X_test(:,2)-X_test(:,5))); %g-z
colour_space_test(:,6)=-2.5*((X_test(:,2)-X_test(:,6))); %g-Y1
colour_space_test(:,7)=-2.5*((X_test(:,2))); %g
colour_space_test(:,8)=-2.5*((X_test(:,8)-X_test(:,7))); %J-Y2
colour_space_test(:,9)=-2.5*((X_test(:,8)-X_test(:,9))); %J-H
colour_space_test(:,10)=-2.5*((X_test(:,8))); %J



% GCSL
filters_used=1:10;

model_train = gmmvar_missing(colour_space_train(:,filters_used),K, options);
model_test = gmmvar_missing(colour_space_test(:,filters_used),K, options);

P_train = getP_new_GMM(colour_space_train(:,filters_used), model_train);
P_test = getP_new_GMM(colour_space_train(:,filters_used), model_test);

%omega = P_test./P_train; 
omega_GMM=(P_test+0.01)./(P_train+0.01);

%omega_GMM=min(omega_GMM,100);
%omega_GMM=max(omega_GMM,0.1);


csvwrite('.../omega_GCSL.txt',omega_GMM)


%%%%%
%%% Find best groups



% Train
n_train=length(X_train(:,1));
best_mixture_train=zeros(1,n_train);

mixture_probability=zeros(n_train,model_train.K);

for i=1:n_train;
    
    
    
    for k=1:model_train.K
        
        mk = model_train.mus(k,:);
        Ck = model_train.Sigmas(:,:,k);
        Ck=0.5*(Ck+Ck');
        
        mixture_probability(i,k)=mvnpdf(colour_space_train(i,filters_used),mk,Ck)*model_train.weights(k);
    
    end
    
    [val, idx] = max(mixture_probability(i,:));
    best_mixture_train(i)=idx;
    
    mixture_probability(i,:)=mixture_probability(i,:)/sum(mixture_probability(i,:));

end

mixture_probability_train=mixture_probability;


% Test
n_test=length(X_test(:,1));
best_mixture_test=zeros(1,n_test);

mixture_probability=zeros(n_test,model_train.K);

for i=1:n_test;
    
    
    for k=1:model_train.K
        
        mk = model_train.mus(k,:);
        Ck = model_train.Sigmas(:,:,k);
        Ck=0.5*(Ck+Ck');
        
        mixture_probability(i,k)=mvnpdf(colour_space_test(i,filters_used),mk,Ck)*model_train.weights(k);
    
    end
    
    [val, idx] = max(mixture_probability(i,:));
    best_mixture_test(i)=idx;
    
    mixture_probability(i,:)=mixture_probability(i,:)/sum(mixture_probability(i,:));

end

mixture_probability_test=mixture_probability;


csvwrite('.../best_mixture_train.txt',best_mixture_train)
csvwrite('.../best_mixture_test.txt',best_mixture_test)

csvwrite('.../mixture_probabilities_train.txt',mixture_probability_train)
csvwrite('.../mixture_probabilities_test.txt',mixture_probability_test)
