%% Inspect the results

% The stimulus is 100 pixels (in both height and weight), and this corresponds to
% 10 degrees of visual angle.  To convert from pixels to degreees, we multiply
% by 10/100.
results = results_quick;
cfactor = 10/100;

% Visualize the location of each voxel's pRF
figure; hold on;
set(gcf,'Units','points','Position',[100 100 400 400]);
cmap = jet(size(results.ang,1));
for p=1:size(results.ang,1)
  xpos = results.ecc(p) * cos(results.ang(p)/180*pi) * cfactor;
  ypos = results.ecc(p) * sin(results.ang(p)/180*pi) * cfactor;
  ang = results.ang(p)/180*pi;
  sd = results.rfsize(p) * cfactor;
  h = drawellipse(xpos,ypos,ang,2*sd,2*sd);  % circle at +/- 2 pRF sizes
  set(h,'Color',cmap(p,:),'LineWidth',2);
  set(scatter(xpos,ypos,'r.'),'CData',cmap(p,:));
end
drawrectangle(0,0,10,10,'k-');  % square indicating stimulus extent
axis([-10 10 -10 10]);
straightline(0,'h','k-');       % line indicating horizontal meridian
straightline(0,'v','k-');       % line indicating vertical meridian
axis square;
set(gca,'XTick',-10:2:10,'YTick',-10:2:10);
xlabel('X-position (deg)');
ylabel('Y-position (deg)');
saveas(gca, '/Users/tong_processor/Desktop/Loic/retinotopy/loic_retino/daveTopup_ling/sub-F019/analysis/prf_location_rel_stimuli','png')
%%

%% Perform some setup
load('/Users/tong_processor/Desktop/Loic/retinotopy/loic_retino/daveTopup_ling/sub-F019/analysis/results_quick_06162023.mat');
load('/Users/tong_processor/Desktop/Loic/retinotopy/loic_retino/daveTopup_ling/sub-F019/code/pipeline/kayWedgeRing.mat');
stimulus=finalstimuli;
data=finaldata;
results=results_quick;
% Define some variables
res = [100 100];                    % row x column resolution of the stimuli
resmx = 100;                        % maximum resolution (along any dimension)
hrf = results.options.hrf;          % HRF that was used in the model
degs = results.options.maxpolydeg;  % vector of maximum polynomial degrees used in the model

% Pre-compute cache for faster execution
[d,xx,yy] = makegaussian2d(resmx,2,2,2,2);


% Prepare the stimuli for use in the model
stimulusPP = {};
for p=1:length(stimulus)
  stimulusPP{p} = squish(stimulus{p},2)';  % this flattens the image so that the dimensionality is now frames x pixels
  stimulusPP{p} = [stimulusPP{p} p*ones(size(stimulusPP{p},1),1)];  % this adds a dummy column to indicate run breaks
end

% Define the model function.  This function takes parameters and stimuli as input and
% returns a predicted time-series as output.  Specifically, the variable <pp> is a vector
% of parameter values (1 x 5) and the variable <dd> is a matrix with the stimuli (frames x pixels).
% Although it looks complex, what the function does is pretty straightforward: construct a
% 2D Gaussian, crop it to <res>, compute the dot-product between the stimuli and the
% Gaussian, raise the result to an exponent, and then convolve the result with the HRF,
% taking care to not bleed over run boundaries.
modelfun = @(pp,dd) conv2run(posrect(pp(4)) * (dd*[vflatten(placematrix(zeros(res),makegaussian2d(resmx,pp(1),pp(2),abs(pp(3)),abs(pp(3)),xx,yy,0,0) / (2*pi*abs(pp(3))^2))); 0]) .^ posrect(pp(5)),hrf,dd(:,prod(res)+1));

% Construct projection matrices that fit and remove the polynomials.
% Note that a separate projection matrix is constructed for each run.
polymatrix = {};
for p=1:length(degs)
  polymatrix{p} = projectionmatrix(constructpolynomialmatrix(size(data{p},4),0:degs(p)));
end

%% Inspect the data and the model fit

% Which voxel should we inspect?  Let's inspect the second voxel.
vx = 60;
vy = 100; %AP
vz = 38/2;

% For each run, collect the data and the model fit.  We project out polynomials
% from both the data and the model fit.  This deals with the problem of
% slow trends in the data.
datats = {};
modelts = {};
for p=1:length(data)
  datats{p} =  polymatrix{p}*squeeze(data{p}(vx,vy,vz,:));
  modelts{p} = polymatrix{p}*modelfun(results.params(1,:,vx*vy*vz),stimulusPP{p});
end

% Visualize the results
figure; hold on;
set(gcf,'Units','points','Position',[100 100 1000 100]);
plot(cat(1,datats{:}),'r-');
plot(cat(1,modelts{:}),'b-');
straightline(300*(1:4)+.5,'v','g-');
xlabel('Time (s)');
ylabel('BOLD signal');
ax = axis;
axis([.5 1200+.5 ax(3:4)]);
title('Time-series data');
