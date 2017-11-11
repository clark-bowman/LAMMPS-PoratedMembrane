%% Pore radius calculation from VMD rendering (front)
radius = zeros(100000, 1);
for i = 1:100000
    % Loop through frames...
    if exist(strcat('analysis.',sprintf('%05d',i-1),'.bmp'),'file') ~= 2
        break
    end
    A=imread(strcat('analysis.',sprintf('%05d',i-1),'.bmp'));
    
    % Read binary images...
    B=(A(:,:,1)==0);
    [x,y]=size(B);
    
    % Restrict to center region of image...
    B=double(B(floor(x/2)-600:floor(x/2)+600,floor(y/2)-600:floor(y/2)+600));
    
    % Smoothing by convolution, then thresholding...
    s = ones(3);
    C = (conv2(B,s) > 5);
    Ct = sum(sum(C));
    
    % If a hole exists, calculate radius in terms of integral approx.
    if (Ct > 1)
        Cx = sum(C)*(1:1203)'/Ct;
        Cy = (1:1203)*sum(C,2)/Ct;
        [Cx2,Cy2] = meshgrid(((1:1203)-Cx).^2,((1:1203)-Cy).^2);
        Gr2 = sum(sum(sqrt(Cx2+Cy2).*C))/Ct;
        radius(i) = 1.5*Gr2;
    end
end
radius = radius(1:i-1);
%% Membrane width calculation from VMD rendering (side)
width = zeros(20, 1);
for i = 1:20
    % Loop through frames...
    if exist(strcat('side.',sprintf('%05d',i-1),'.bmp'),'file') ~= 2
        break
    end
    A=imread(strcat('side.',sprintf('%05d',i-1),'.bmp'));
    
    % Read binary images...
    B=(A(:,:,1)==0);
    [x,y]=size(B);
    
    % Restrict to center region of image...
    B=double(B(floor(x/2)-600:floor(x/2)+600,floor(y/2)-600:floor(y/2)+600));
    
    % Smoothing by convolution, then thresholding...
    s = ones(3);
    C = (conv2(B,s) < 5);
    
    % Loop over rows and save instantaneous width...
    width(i) = 9999;
    for j = 2:1202
        tempoi = find(C(j,2:1202));
        boxx = max(tempoi) - min(tempoi);
        if boxx < width(i)
            width(i) = boxx;
        end
    end
end

% Calculate mean membrane width.
width = mean(width(1:i-1));