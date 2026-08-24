function [thetas,the_t,spread_matrix]=GetVelocityRadon(data,windowsize)
%
%function applies a sliding time window over the RBCV line scan data,
%within each window a radon transform is applied to find the angle of the
%lines. The windows overlap to provide smoothing and better accuracy when
%taking the 'average angle' of the RBC shadows within each window.
%this code is adapted directly from Drew, 2010
%adapted by Kira, Dori and Orla, 2017, updated by Kira May 2018
%
%INPUTS
%data - the matrix of time X space data 
%windowsize - number of lines to use in estimating velocity.  must be a
%multiple of 4
%
%OUTPUTS
%thetas - the time varying angle of the space-time image
%the_t - time points of the angle estimates (in lines) 
%spreadmatrix - matix of variances as a function of angles at each time 
%point
%
%if more information is required regarding how this function works, please
%refer to the detailed methods paper referenced below: 
%Drew et al. (2010). J Comput Neurosci. 29, 5-11. 

%define step size for sliding window
%this is because the data will be averaged across the window, then step 1/4
%of the window and average again - this means data is overlapped in
%averaging by multiple windows, and provides a smoothing effect
%divide the no of pixels in the linescan by 4
stepsize=round(.25*windowsize); %pixels
nlines=size(data,1); %time, i.e. no Lines
npoints=size(data,2); %space, i.e. Pixels
%calculate how many steps can be taken with the step size specified and
%size of the data - this will be used for looping the data 
nsteps=floor(nlines/stepsize)-3;
%find the edges
%angle of line will not be greater than 180 degrees
%check each angle (between 0-180), to find slope of line
angles=(0:179); 
%once found crude angle, use this for FINE estimates of angle 
%(in 0.25 increments)
angles_fine=-2:.25:2; 

%create empty matrices to fill
%no steps (i.e. time windows) * how many angles to assess
spread_matrix=zeros(nsteps,length(angles));
%no steps (i.e. time windows) * how many fine angles to assess
spread_matrix_fine=zeros(nsteps,length(angles_fine));
%no steps (i.e. time window)- as will get precise angle for each step 
%(in degrees)
thetas=zeros(nsteps,1);

%hold matrix will be multiplied by the mean, for normalising later
hold_matrix=ones(windowsize,npoints);
%NaNs for every step 
the_t=NaN*ones(nsteps,1); 

%loop for every step 
for k=1:nsteps

    %some sort of time vector? (used to create time later) - use to get
    %sample rate
    the_t(k)=1+(k-1)*stepsize+windowsize/2; 
    
    % take info from Continuous data
    % loop through data in steps of step size * 4 
    % this will overlap by one time box on each loop 
    % 0:T/4 time box = stepsize+windowsize/2 - i.e. first box
    % 4 steps later (i.e. 4 time boxes on, 4*stepsize)
    % - (k-1)*stepsize+windowsize
    data_hold=data(1+(k-1)*stepsize:(k-1)*stepsize+windowsize,:);
    
    % normalise to remove inhomogeneities in illumination in space and to
    % remove the baseline intensity
    data_hold=data_hold-mean(data_hold(:))*hold_matrix;%subtract the mean
    
    %send time boxed continuous data into radon transform, with the 180
    %possible angles 
    % The Radon transform is the projection of the image intensity along a 
    % radial line oriented at a specific angle.
    % THETA is a vector, then R is a matrix in which each column is the 
    % Radon transform for one of the angles in THETA. If you omit THETA, it
    % defaults to 0:179.
    radon_hold=radon(data_hold,angles);%radon transform
    
    %find the variance of the radon tranform values    
    spread_matrix(k,:)=var(radon_hold);%take variance
    [~, the_theta]=max(spread_matrix(k,:));%find max variance
    thetas(k)=angles(the_theta);     
    %re-do radon with finer increments around first estimate of the maximum
    radon_hold_fine=radon(data_hold,thetas(k)+angles_fine);
    spread_matrix_fine(k,:)=var(radon_hold_fine);
    clear the_theta; 
    [~, the_theta]=max(spread_matrix_fine(k,:));
    thetas(k)=thetas(k)+angles_fine(the_theta);
    
    %figures 
    %if debugging the code, can use to check what is going on
    %turned off from automatically plotting
    if false
        if rem(k, 1000) == 0 %plot every 1000 steps
            
            %not sure this is the right way to plot the line/angle
            %calculated for the time window
            figure(k)
            %make triangle to plot onto line scan 
            %(i.e. can see if calc is right)
            ange = -1*(thetas(k)-90);
            a = tan(ange);
            x = [1:npoints];
            y = x*a;
            y = y - min(y)+1;

            lineim = zeros(size(data_hold));
            for i = 1:npoints
                lineim(round(y(i)), npoints+1-i) = 1;
            end
            
            lineim = lineim(1:size(data_hold,1), 1:size(data_hold,2));
            
            subplot(2,1,1);
            imagesc(imfuse(data_hold, lineim, 'ColorChannels',[0 2 1]));
            colormap gray;
            title('data');
            subplot(2,1,2);
            imagesc(radon_hold);
            colormap gray;
            title('radon transform');
            
        end %end of plotting every 1000 steps
    end %end of suppressing plots
    
end %end of steps loop 

thetas=-1*(thetas-90); %rotate theta calculated 

end %end of function

