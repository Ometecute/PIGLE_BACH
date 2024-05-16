%%MD SIMULATION BELOW%%
%this section should have been completed in the previous tutorial 
%------------------------------%
clear;
eta=1; %1/picoseconds
mass=1; %amu
Ts=200; %Kelvin
Kb=0.8314; %A^2 amu ps^-2 K^-1
dt=.01; %picoseconds
dt2=dt^2; %this will be handy in the equation of motion to speed it up
nstep=100000; %number of timesteps in the simulation
ndim=2; %number of dimensions

scaling=sqrt(1/dt*eta*Kb*Ts*mass)/mass; %there is an additional 1/m included becaute we are going to use units of accelleration.

imp=sqrt(2)*scaling.*randn(ndim,nstep); %sqrt 2 because the impulse scaling is the STD not FWHM for randn. 

v(1:ndim,1:nstep)=ones.*nan;
r(1:ndim,1:nstep)=ones.*nan;

r(:,1)=zeros; % origin
v(:,1)=zeros; % start at rest

for ind=2:nstep  
    r(:,ind)=r(:,ind-1)+v(:,ind-1)*dt+dt2/2.*(-1*eta.*v(:,ind-1)+imp(:,ind-1));
    v(:,ind)=(r(:,ind)-r(:,ind-1))/dt;
end 

plot(r(1,:),r(2,:))
KE=mass*(v(1,:).^2+v(2,:).^2);
disp('T_act='); Ta=mean(KE)/Kb; disp(Ta);
disp('T_set=');disp(Ts);
%------------------------------%

%%
%%where ISF calculation tutorial starts
%------------------------------%

%linspace creates evenly spaced points between 0 and 5
KVec=linspace(0,5,30);

%using the formula given in the tutorial
%notice the transpose sign, 'linspace' returns a row vector
AKt=exp(-1i*(KVec'*r(1,:)+KVec'*r(2,:)));

%The first entry of fft is the function to be Fourier transformed
%The second entry specifies the number of points used to take discrete Fourier transform (DFT). We let MatLab run automatically and did not specify any value.
%The third and last entry specifies the direction along which dimension to take the Fourier transform from. In this case, 2 tells MatLab to take each row as input.
AKw=fft(AKt,[],2);

%We wish to investigate A of the first numK values of Delta K...
numK=10;
%...of the numTime time steps
numTime=100;

IKt=ifft(AKw.*conj(AKw),[],2);

%Using figure environment
figure;
%Suppress figure refreshing
hold on
%Set labels
xlabel('Time/ps')
ylabel('$\Re(A(\Delta\bm{K},t))$','Interpreter','latex')
%Go through numK's Delta K
for i = 1:numK
%'real' applies on a complex valued array returns the real part of each entry
%Use linspace for evenly spaced points across the numTime time range
    plot(linspace(0, dt*numTime, numTime), real(AKt(i, 1:numTime)))
end

%Start a figure afresh
figure;
hold on
xlabel('Time/ps')
ylabel('$\Im(A(\Delta\bm{K},t)$','Interpreter','latex')
for i = 1:numK
%Seek the imaginary part using the 'imag' function.
    plot(linspace(0, dt*numTime, numTime), imag(AKt(i, 1:numTime)))
end

%Start afresh
figure;
hold on
%Set label
xlabel('$\Re(A(\Delta K,t))$','Interpreter','latex')
ylabel('$\Im(A(\Delta K,t))$','Interpreter','latex')
for i = 1:numK
%Extract the each parts into arrays to plot against each other
    plot(real(AKt(i,1:numTime)),imag(AKt(i,1:numTime)))
end
axis equal

%------optional------%
%plots the argument and modulus of the function 
figure;
hold on
xlabel('Time/ps')
ylabel('$|A(\Delta K,t)|$','Interpreter','latex')
for i = 1:numK
    plot(linspace(0,dt*numTime,numTime),abs(AKt(i,1:numTime)))
end
axis equal

figure;
hold on
xlabel('Time/ps')
ylabel('$Arg(\Delta K,t)$','Interpreter','latex')
for i = 1:numK
    plot(linspace(0,dt*numTime,numTime),angle(AKt(i,1:numTime)))
end
axis square
%--------------------%

%Equations based on the analytic solution provided in the tutorial
%The '@' operator is a short hand that creates local functions
phi      = @(t) 1-exp(-eta*abs(t));
chi      = @(DK) DK*sqrt(Kb*Ts/mass)/eta;
isf      = @(t,DK) exp(-chi(DK)^2*(eta*abs(t)-phi(t)));
%isf_fast is the ballistic part of the motion that uses approximation derived in the tutorial
isf_fast = @(t,DK) exp(-chi(DK)^2*eta^2*t.^2/2);

for ind =2:length(KVec)
    figure

    plot(linspace(-dt*nstep/2,dt*nstep/2,nstep),fftshift(real(IKt(ind,:))./nstep),'Color','#800000');
    %The 'prepareCurveData' function comprises of many steps
    %   prepares rows, or even matrices, into columns
    %   takes the real part of each entry
    %   convert the data type of all entries into double
    [xData, yData] = prepareCurveData(linspace(-dt*nstep/2,dt*nstep/2,nstep), fftshift(real(IKt(ind,:))./nstep));
    %Plot the data and specify the shape of the marker and the colour accordingly
    plot(xData,yData,'x','Color','#48399e'); hold on;
    
    %The plot has been shifted so it centres at the origin
    %Hence we uses linspace from -100 to 100
    %The function 'isf' defined above takes two inputs
    %   time range - -100 to 100 as centralised
    %   KVec, which is the momentum transfer that has been defined before
    %Specify the plot as red solid line with the attributes below    
    %plot(linspace(-100,100,500),isf(linspace(-100,100,500),KVec(ind)),'-','Color','#800000','LineWidth',3);
    
    xlabel('Time (centralised)')
    ylabel('ISF')
    xlim([-100 100]);

    %The tilde operator cleverly negates the range of points being selected
    %Hence, excludePoints is now all points outside the range of [-.01,.01]
    excludedPoints=~excludedata(xData,yData,'Domain',[-.01 .01]);
    % Set up fittype and options.
    ft = fittype( 'a*exp(-b*abs(x))', 'independent', 'x', 'dependent', 'y' );
    %Specify the fitting method, i.e. least square method
    %The user usually does not need to specify the method and MatLab chooses automatically
    %These lines are being included to show the versatility of the function and the many attribute that can be exploit when needed
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    %Do not show the intermediate steps
    opts.Display = 'Off';
    %The coefficients has to be non negative
    opts.Lower = [0 0];
    %Initialise the first pair of coefficients
    opts.StartPoint = [0.5 20];
    %Tell MatLab to ignore the points outside range of consideration
    opts.Exclude = excludedPoints;
    
    % Fit model to data.
    [fitresult, gof] = fit( xData, yData, ft, opts );
    %Store the coefficients, b, of varying Delta K, into an array called DKs
    DKs(ind)=fitresult.b;
    
    %Specify the fitting type as a Gaussian function
    %Different modes of Gaussian curve is available to MatLab
    %The '1' in 'gauss1' tells MatLab that this is an average Gaussian distribution with only one term
    %If instead it was 'gauss2', MatLab will fit the data to the sum of two different Gaussian distributions
    f = fit(xData,yData,'gauss1');
    %Store the standard deviations, sigma's, in an array, derived from c1
    sigmas(ind) = f.c1/sqrt(2);
    
end

%Plotting standard deviation, sigma, against momentum transfer
%Specify the marker's shape, size, and colour
figure
plot(KVec(2:end),DKs(2:end),'x','Color','#800000','MarkerSize',15)
xlabel('$\Delta \mathbf{K}$','Interpreter','latex')
ylabel('$\alpha$','Interpreter','latex')

%Plotting standard deviation, sigma, against momentum transfer
%Specify the marker's shape, size, and colour
figure
plot(KVec(2:end),sigmas(2:end),'x','Color','#800000','MarkerSize',15)
xlabel('$\Delta \mathbf{K}$','Interpreter','latex')
ylabel('$\sigma$','Interpreter','latex')
