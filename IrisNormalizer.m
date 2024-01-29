function IrisNormalizer(inputImageURL, minRadiusPupil, maxRadiusPupil,minRadiusIris,maxRadiusIris, outputDirectory)
    % Task 1: Read the input image
    originalImage = imread(inputImageURL);

    % Task 2: Convert the color image to an 8-bit gray-level image
    grayImage = rgb2gray(originalImage);
    
    
    % Task 3: Find pupil and iris boundaries with Circular Hough Transform
    [xc1,yc1, radiiP] = CircularHoughTransform(grayImage, minRadiusPupil , maxRadiusPupil);

    
    [xc2,yc2, radiiI] = CircularHoughTransform(grayImage, minRadiusIris, maxRadiusIris);

    xc = (xc1 + xc2)/2;
    yc = (yc1 + yc2)/2;
    
    red = originalImage(:,:,1);
    green = originalImage(:,:,2);
    blue = originalImage(:,:,3);
   
    % Task 4: Daugman Rubber Sheet Model
    unwrappedGray = DaugmanRubberSheetModel(grayImage, xc,yc, radiiP,radiiI);

    unwrappedRed = DaugmanRubberSheetModel(red, xc,yc, radiiP,radiiI);
    unwrappedGreen = DaugmanRubberSheetModel(green, xc,yc, radiiP,radiiI);
    unwrappedBlue = DaugmanRubberSheetModel(blue, xc,yc, radiiP,radiiI);
    unwrappedColor = cat(3,unwrappedRed,unwrappedGreen,unwrappedBlue);
    figure,imshow(unwrappedColor)

    % Task 5: Output unwrapped images
    outputGrayFilename = fullfile(outputDirectory, 'unwrapped_gray.pgm');
    outputColorFilename = fullfile(outputDirectory, 'unwrapped_color.ppm');

    imwrite(unwrappedGray, outputGrayFilename);
    imwrite(unwrappedColor, outputColorFilename);
end
function [xc,yc, radii] = CircularHoughTransform(grayImage, minRadius, maxRadius)
    % Step 1: Initialization
    [rows, cols] = size(grayImage);
    accumulator =zeros(rows, cols,maxRadius-minRadius+1);
    
    % Step 2: Voting
    
    gaussianFilter = fspecial('gaussian',35, 25);
    blurImage = imfilter(grayImage, gaussianFilter,'symmetric');


    edgeImage = edge(blurImage, 'canny',[]);
    
    
    
    [edgeRows, edgeCols] = find(edgeImage);
    

    for i = 1:numel(edgeRows)
        
        for radius = minRadius:maxRadius
            for theta = 0:360
                xc = round(edgeCols(i) - radius * cos(theta));
                yc = round(edgeRows(i) - radius * sin(theta));
                
                if xc >= 1 && xc <= cols && yc >= 1 && yc <= rows
                    accumulator(yc, xc,radius - minRadius + 1) = accumulator(yc, xc,radius - minRadius + 1) + 1;
                end
            end
        end
    end
    
     
    figure, imshow(grayImage)
    % Obtain the coordinates of the maximum value of the accumulator
    [y,x,r]=ind2sub(size(accumulator),find(accumulator==max(max(max(accumulator)))));
    r=r+minRadius;
    
    xc = x;
    yc =y;
    radii = r;

    disp(x)
    disp(y)
    disp(r)
    % Draw over the picture the generated circle
    hold on
    th = 0:pi/50:2*pi;
    xunit = r * cos(th) + x;
    yunit = r * sin(th) + y;
    plot(xunit, yunit);


    
end



  


function [unwrapped_image] = DaugmanRubberSheetModel(img, xc,yc,radiiP,radiiI)
 % Define polar coordinates
   
   
    disp(xc)
    disp(yc)
    r =  0:1/180:1;

    
    sAngles = 360;
    angles =(0:pi/sAngles:pi-pi/sAngles) + pi/(2*sAngles);
    nAngles = length(angles);
    
% Calculate pupil points and iris points that are on the same line
    x1 = ones(size(angles))*xc;
    y1 = ones(size(angles))*yc;
    x2 = xc + 10*sin(angles);
    y2 = yc + 10*cos(angles);
    dx = x2 - x1;
    dy = y2 - y1;
    slope = dy./dx;
    intercept = yc - xc .* slope;
    
    xout = zeros(nAngles,2);
    yout = zeros(nAngles,2);
    for i = 1:nAngles
            % Coefficients
        A = 1 + slope(i)^2;
        B = 2 * (slope(i) * intercept(i) - slope(i) * yc - xc);
        C = xc^2 + yc^2 + intercept(i)^2 - 2 * intercept(i) * yc - radiiP^2;
    
        % Solve quadratic equation
        discriminant = B^2 - 4 * A * C;
    
        if discriminant < 0
            % No intersection
            xOut = [];
            yOut = [];
        elseif discriminant == 0
            % One intersection point
            xOut = -B / (2 * A);
            yOut = slope(i) * xOut + intercept(i);
        else
            % Two intersection points
            xout1 = (-B + sqrt(discriminant)) / (2 * A);
            xout2 = (-B - sqrt(discriminant)) / (2 * A);
    
            yout1 = slope(i) * xout1 + intercept(i);
            yout2 = slope(i) * xout2 + intercept(i);
    
            xOut = [xout1, xout2];
            yOut = [yout1, yout2];
        end
        xout(i,:) = xOut;
        yout(i,:) = yOut;
    end


    % Get samples on limbus boundary
    xRightIris = yc + radiiI * cos(angles);
    yRightIris = xc + radiiI * sin(angles);
    xLeftIris = yc - radiiI * cos(angles);
    yLeftIris = xc - radiiI * sin(angles);
    
    % Get samples in radius direction
    xrt = (1-r)' * xout(:,1)' + r' * yRightIris;
    yrt = (1-r)' * yout(:,1)' + r' * xRightIris;
    xlt = (1-r)' * xout(:,2)' + r' * yLeftIris;
    ylt = (1-r)' * yout(:,2)' + r' * xLeftIris;

  
   
    % Create Normalized Iris Image
    
    image = uint8(reshape(interp2(double(img),[xrt(:);xlt(:)],[yrt(:);ylt(:)]),length(r), 2*length(angles)));
    
   
  
   
        
    
    
    


    unwrapped_image = image;
   

  
   

    
   
  
end
