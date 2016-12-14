%
% Script demonstrating the use of the function
% function_ClipPoisGaus_denoising    [version 2.11 released 8 December 2008]
%
% This script demonstrates the algorithm described in the publication
% Foi, A., "Practical denoising of clipped or overexposed noisy images", Proc. 16th European Signal Process. Conf., EUSIPCO 2008, Lausanne, Switzerland, August 2008.
%
% Required functions:   BM3D.m and related dll/mex files,  available from http://www.cs.tut.fi/~foi/GCF-BM3D,  or an alternative denoising filter.
%                       function_ClipPoisGaus_stdEst2D  and  ml_fun_reg,  available from http://www.cs.tut.fi/~foi/sensornoise.html
% Optional function:    function_Errors    -  computes error criteria,  available from http://www.cs.tut.fi/~lasip/2D/
%
% Alessandro Foi,  firstname.lastname@tut.fi  www.cs.tut.fi/~foi  Tampere University of Technology  2008
% ------------------------------------------------------------------------------------------------------

clear all
close all

disp(' ')
if 0   %% RAW DATA
    ztilde=double(imread('image_fuji_staircase1224x922.tiff'));   % this is a subsampled green channel raw-data from a Fujifilm FinePix S9600 camera
    % ztilde=double(imread('image_fuji_bottle1224x922.tiff'));      % this is a subsampled blue channel raw-data from a Fujifilm FinePix S9600 camera
    ztilde=ztilde/15840;             %  or divide by max(ztilde(:)), if that can be taken as the clipping value
    clipping_below=1;  %%%% on off   %  RAW-DATA IS ASSUMED TO BE CLIPPED FROM ABOVE AND BELOW
    clipping_above=1;  %%%% on off
    prior_density=0;                 %  type of prior density to use for ML    (0)
    %                                %    0: zero_infty uniform prior density (R+);  (default, use this for raw-data)
    %                                %    1: zero_one uniform prior density [0,1];
    %                                %    2: -infty_infty uniform prior density (R);
    add_noise=0;                     %  no need to add noise on already noisy data
elseif 1   %% SYNTHETIC NOISE ON TEST IMAGES
    add_noise=1;                     %  add noise and clip noisy data
    %% CHOOSE IMAGE AND NOISE SETTINGS
    if 0
        y=im2double(imread('image_man1024.tiff'));
        chi=30;        a=1/chi;
        b=0.1^2;
        clipping_below=1;  %%%% on off   %  RAW-DATA IS ASSUMED TO BE CLIPPED FROM ABOVE AND BELOW
        clipping_above=1;  %%%% on off
        prior_density=1;                 %  type of prior density to use for ML    (0)
        %                                %    0: zero_infty uniform prior density (R+);  (default, use this for raw-data)
        %                                %    1: zero_one uniform prior density [0,1];
        %                                %    2: -infty_infty uniform prior density (R);
    elseif 0
        y=im2double(imread('image_testpat1024.tiff'));
        a=0;
        b=0.2^2;
        clipping_below=1;  %%%% on off   %  RAW-DATA IS ASSUMED TO BE CLIPPED FROM ABOVE AND BELOW
        clipping_above=1;  %%%% on off
        prior_density=1;                 %  type of prior density to use for ML    (0)
        %                                %    0: zero_infty uniform prior density (R+);  (default, use this for raw-data)
        %                                %    1: zero_one uniform prior density [0,1];
        %                                %    2: -infty_infty uniform prior density (R);
    elseif 0
        % y=im2double(imread('image_cameraman256.png'));
        y=im2double(imread('image_lena512.png'));
        chi=40;        a=1/chi;
        b=0;
        clipping_below=1;  %%%% on off   %  RAW-DATA IS ASSUMED TO BE CLIPPED FROM ABOVE AND BELOW
        clipping_above=0;  %%%% on off
        prior_density=1;                 %  type of prior density to use for ML    (0)
        %                                %    0: zero_infty uniform prior density (R+);  (default, use this for raw-data)
        %                                %    1: zero_one uniform prior density [0,1];
        %                                %    2: -infty_infty uniform prior density (R);
    else %% some smooth/analytic signal
        sizez=512; % do not exceed 1500, please!
        [XMATRIX,YMATRIX]=meshgrid(1:sizez,1:sizez);
        XMATRIX=XMATRIX-min(XMATRIX(:));  YMATRIX=YMATRIX-min(YMATRIX(:));
        XMATRIX=XMATRIX/max(XMATRIX(:));  YMATRIX=YMATRIX/max(YMATRIX(:));
        y=sin(2*pi*XMATRIX);
        y=0.7*y+0.5;
        a=0;
        b=0.15^2;
        clipping_below=1;  %%%% on off   %  RAW-DATA IS ASSUMED TO BE CLIPPED FROM ABOVE AND BELOW
        clipping_above=1;  %%%% on off
        prior_density=2;                 %  type of prior density to use for ML    (0)
        %                                %    0: zero_infty uniform prior density (R+);  (default, use this for raw-data)
        %                                %    1: zero_one uniform prior density [0,1];
        %                                %    2: -infty_infty uniform prior density (R);
    end
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ESTIMATE NOISE PARAMETERS or USE KNOWN PARAMETERS
estimate_parameters=1;            %  estimate parameters from image  (otherwise use given fitparamstrue=[a,b,c,...])

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% TEXT and FIGURES  (set both to zero for silent processing)
text_verbosity_main=3;     % print to screen nothing (0), or some information about the ongoing processing (1-3)   [default 3]
figure_verbosity_main=1;   % show figures?  0: no figures,  1: show some figures  [default 1]



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% NOISE ESTIMATION ALGORITHM MAIN PARAMETERS   (needed in case estimate_parameters==1)
%% (these parameters are used by the function function_ClipPoisGaus_stdEst2D)

polyorder=1;                     %  polynomial order of the root-polynomial model to be estimated [default 1, i.e. linear]
%                                %   Note: a large order results in overfitting and difficult and slow convergence of the ML optimization, therefore it is recommended to increase speed_factor when polyorder>1.
%                                %         Further, higher order models typically result in problems of invertibility of the stabilizing, debiasing, and declipping transformations.
variance_power=1;                %  power of the variance [default 1, i.e. affine/linear variance]
%                                %  the standard-deviation function has the form \sigma=sqrt(a*y^polyorder+b*y^(polyorder-1)+c*y^(polyorder-2)+...).^variance_power, where y is the unclipped noise-free signal.
%                                %   The usual Poissonian-Gaussian model has the form \sigma=sqrt(a*y+b).
median_est=1;                    %  0: sample standard deviation;  1: MAD   (1)
LS_median_size=1;                %  size of median filter for outlier removal in LS fitting (enhances robustness for initialization of ML) 0: disabled  [default 1 = auto]
tau_threshold_initial=1;         %  (initial) scaling factor for the tau threshold used to define the set of smoothness   [default 1]

% prior_density=0;               %  type of prior density to use for ML    (0)
%                                %    0: zero_infty uniform prior density (R+);  (default, use this for raw-data)
%                                %    1: zero_one uniform prior density [0,1];
%                                %    2: -infty_infty uniform prior density (R);

level_set_density_factor=1;      %   density of the slices in for the expectations   [default 1 ( = 600 slices)]   (if image is really small it should be reduced, say, to 0.5 or less)
integral_resolution_factor=1;    %   integral resolution (sampling) for the finite sums used for evaluatiing the ML-integral   [default 1]
speed_factor=1;                  %   factor controlling simultaneously density and integral resolution factors  [default 1] (use large number, e.g. 1000, for very fast algorithm)

text_verbosity=1;                %  how much info to print to screen 0: none, 1: little, 2: a lot
figure_verbosity=3;              %  show/keep figures?        [default 1]
%                                     0: no figures
%                                     1: only figure with final ML result is shown and kept
%                                     2: few figures are shown during processing but none kept
%                                     3: few figures are shown but only figure with final ML result is kept
%                                     4: show and keep all figures

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% STRING FOR CALLING A DENOISING ALGORITHM (MAIN PARAMETERS)

use_bm3d=1;                      % use BM3D denoising algorithm  [default 1]
bm3d_profile='np';               % complexity profile of BM3D denoising algorithm  [default 'np']
%                                    'np'  normal profile
%                                    'lc'  low-complexity profile (faster, slightly lower quality)
bm3d_string=['[PSNR, Efztilde] = BM3D(1,fztilde,sigma*255,''',bm3d_profile,''');'];    %   default:  bm3d_string=['[PSNR, Efztilde] = BM3D(1,fztilde, sigma*255,''',bm3d_profile,''');']

other_string='[Efztilde]= other_denoising_algorithm(fztilde,sigma);' ;  % if use_bm3d==0




%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%  END OF DEMONSTRATION SCRIPT SETTIGS  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



if add_noise
    %% ADD NOISE & CLIP
    if 0
        disp(' !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! fixed pseudonoise !!!!!!!!!!!!!!!!!!');
        disp(' ');
        randn('seed', 0);
        rand('seed', 0);
        randn('state', 0);
        rand('state', 0);
    end

    if a==0   % no Poissonian component
        z=y;
    else      % Poissonian component
        chi=1/a;
        z=poissrnd(max(0,chi*y))/chi+min(y,0);      %%% NOTE!!!!  Whenever y is not positive, it is not possible have poissonian noise for the negative samples!!!
    end
    if b<0
        disp('  !!! The parameter b has to be non-negative !!!   (setting b=0)')
        b=0;
    end
    z=z+sqrt(b)*randn(size(y));   % Gaussian component

    if 0 %% REVERSE  (THIS IS USEFUL TO MAKE TESTS ABOUT ABSORPION IMAGING)
        y=1-y;
        z=1-z;
        b=a+b;
        a=-a;
        clipping_aboveb=clipping_above;
        clipping_above=clipping_below;
        clipping_below=clipping_aboveb;
        clear clipping_aboveb
    end

    % CLIPPING
    if clipping_above
        z=min(z,1);
    end
    if clipping_below
        z=max(0,z);
    end
    ztilde=z;
    clear z
    check_errors=1;                  %  compute PSNR, etc. against original image
    fitparamstrue=[a b];
else
    check_errors=0;                  %  original image is not available and estimation errors cannot be checked
end
clear XMATRIX YMATRIX


if figure_verbosity_main
    figure, imshow(ztilde),title('ztilde   noisy raw observation'), colorbar
end
if check_errors
    if exist('function_Errors','file')==2||exist('function_Errors','file')==6
        disp('Error criteria for observed noisy image \tilde{z} (or z)')
        function_Errors(y,ztilde,1+text_verbosity_main>0);
        disp(' ');
    end
end


if estimate_parameters==0&&~exist('fitparamstrue','var')
    disp(' ! standard-deviation function parameters are unknown and will be estimated.');
    estimate_parameters=1;
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%  PREPARING PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if use_bm3d
    denoise_string=bm3d_string;
else
    denoise_string=other_string;
end
if estimate_parameters
    noise_est_vector=[polyorder,variance_power,median_est,prior_density,LS_median_size,tau_threshold_initial,level_set_density_factor,integral_resolution_factor,speed_factor,text_verbosity,figure_verbosity];
    variance_power=[];
else
    noise_est_vector=fitparamstrue;
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%  CALLING MAIN FUNCTION  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[y_hat,Eztilde,fitparams]=function_ClipPoisGaus_denoising(ztilde,clipping_below,clipping_above,variance_power,noise_est_vector,denoise_string,[],[],[],1,[],text_verbosity_main,figure_verbosity_main);
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% D(ztilde)
if check_errors
    if exist('function_Errors','file')==2||exist('function_Errors','file')==6
        disp('Error criteria for denoised image D(ztilde)')
        function_Errors(y,Eztilde,1+text_verbosity_main>0);
        disp(' ');
    end
end
if figure_verbosity_main
    figure,imshow(Eztilde),  title('DENOISED  D(ztilde)'), colorbar
end

%% \hat{y}
if check_errors
    if exist('function_Errors','file')==2||exist('function_Errors','file')==6
        disp('Error criteria for denoised and declipped image \hat{y}')
        function_Errors(y,y_hat,1+text_verbosity_main>0);
        disp(' ');
    end
end
if figure_verbosity_main
    figure,imshow(y_hat), title('DENOISED and DECLIPPED \\hat\{y\}'), colorbar
    figure,imshow(y_hat,[]), title('DENOISED and DECLIPPED \\hat\{y\}  [stretched]'), colorbar
    pause(eps)
end

%% cross-sections
if figure_verbosity_main
    if check_errors
        linenn=round(1024/6);
        linenn=min(linenn,size(ztilde,1));
        figure,plot([y(linenn,:);ztilde(linenn,:);Eztilde(linenn,:);y_hat(linenn,:)]')
        legend('y','ztilde', 'D(ztilde)','\\hat\{y\}')
    else
        linenn=230;
        linenn=min(linenn,size(ztilde,1));
        figure,plot([ztilde(linenn,:);Eztilde(linenn,:);y_hat(linenn,:)]')
        legend('ztilde', 'D(ztilde)','\\hat\{y\}')
    end
    xlim([1,size(ztilde,2)]);
    legend('location','best');
    grid on
    pause(eps)
end

