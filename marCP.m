
function recon = marCP( sinogram, thetas,nDetectors, dSize, cx, cy, ...
  Nx, Ny, dx, dy, window )

  % Sigma and tau are parameters
  nrmK = sqrt(6);
  sigma = 1/nrmK;
  tau = 1/nrmK;

  scale = [3 3 3 4 4]; % choose size of shearlets at each scale
  qmf1 = MakeONFilter('Symmlet',4); % choose a wavelet filter
  qmf2 = MakeONFilter('Symmlet',4); % choose a wavelet filter

  ndir = 0; % number of directions (2^(ndir+2)+2 directions)      
  alpha = (2^(ndir+2)+2);

  % Define function handles -------------------------------------------
  
  applyK = @(x) ctBackProject( sinogram, thetas, dSize, cx, cy, Nx, Ny, ...
    dx, dy, 'iso' );
  applyKTrans = @(y) ctRadon( img, delta, nDetectors, dSize, thetas, ...
    'iso' );

  computeYnp1 = @(yn,xnBar,sigma) ...
    computeYnp1_deblurImageTightFrame_pockChambolle_shearlet( ...
    yn, xnBar, sigma, applyK );

  computeCost = @(x) computeCost_inpaintTightFrame_arrayVersion( ...
    x, rho, R, img, applyK );

  % -------------------------------------------------------------------


  % Perform Metal Artifact Reduction

  reconSino = chamboullePock( sinogram, applyK, applyKTrans, computeYnp1, ...
    sigma, tau, computeCost );

  recon = ctIRadon( reconSino, thetas, dSize, cx, cy, Nx, Ny, ...
    dx, dy, window );


  %-------------------------------------------------------------------
  
  function out = chamboullePock( img, applyK, applyKTrans, computeYnp1, ...
    sigma, tau, computeCost )

    p = inputParser;
    p.addRequired('img');
    p.addRequired('applyK');
    p.addRequired('applyKTrans');
    p.addRequired('computeYnp1');
    p.addOptional('computeCost','none');

    p.parse( img, applyK, applyKTrans, computeYnp1, computeCost );
    inputs = p.Results;
    computeCost = inputs.computeCost;

    xn = img;
    xnBar = xn;
    yn = applyK(xn);

    theta = 1;

    costs = zeros(1,maxIter+1);
    for n = 1:maxIter

      %%% Compute ykp1.
      ynp1 = computeYnp1(yn,xnBar,sigma); %here
      xnp1 = computeXnp1_inpaintImageTightFrame_pockChambolle( ...
        xn, ynp1, tau, applyKTrans, R, img, rho );
      xnp1Bar = xnp1 + theta*(xnp1 - xn);

      yn = ynp1;
      xn = xnp1;
      xnBar = xnp1Bar;    

      if nargin > 4 costs(n+1) = computeCost(xn); end;
      if mod(n,showTrigger) == 0
        if createFigFlag == 1
          figure('Name','xk inside CP')
          hxn = gcf;
          createFigFlag = 0;
        end

        showx(xn,hxn);               
      end    

    end
    
    out = xn;
  end

end