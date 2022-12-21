clear all; close all; clc;

%%%%  Project: Barcode Scan
%%%%  Student ID: 27952827

%% === STEP 1 - IMPORTING THE IMAGE ===
% I = imread('IMG_20200317_151526.jpg');  % bar code 1
I = imread('IMG_20200317_151537.jpg');  % bar code 2
% I = imread('IMG_20200317_151803.jpg');  % bar code 3
J = imresize(I, 0.35);

%% === STEP 2 - RGB>GRAYSCALE CONVERSION ===
R = J(:,:,1); G = J(:,:,2); B = J(:,:,3);
grayscale = (0.3*R) + (0.59*G) + (0.114*B);    % Weighted Method

clear R G B J I;      % clearing variables that are not needed anymore
%% === STEP 3 - IMAGE BINARIZATION ===
T = 130;        % Threshold 

black = (grayscale < T);        % Splitting black & white bits
white = (grayscale >= T);

barCodeImg = zeros( size(grayscale) );     % Pre-allocation
barCodeImg(black) = 0; barCodeImg(white) = 1; % Binarizing based on T value

clear black white grayscale;
%% === STEP 4 - DETECTING THE CENTER-POINT ===
rowMid = ceil( length( barCodeImg(:,1) ) /2 );  % Index for the row at center
whiteBars = 0; whitePixels = zeros(15,1);   % Pre-allocation
blackBars = 0; blackPixels = zeros(15,1);

% Extracting 3 rows from the center. 
% I chose 3 instead of 1 to take the average to account for incorrect bits
midSection = barCodeImg( rowMid-1:rowMid+1, : ); 
c = 1;      % counter for going through columns 1-by-1
i = 1; j = 1;   % counters for indexing black and white pixels

while blackBars < 15    % while the counter hasn't reached the guard bar
    window = midSection(:,c);   % Extracting one column
    
    % Upon encountering first black bar
    if sum( window ) ~= 3       % if it's not all white bits
        while sum( window ) < 2     % checks and counts black bits 1-by-1
            blackPixels(i) = blackPixels(i) + 1;
            
            c = c + 1;
            window = midSection(:,c);
        end
        blackBars = blackBars + 1;
        i = i + 1;
    end
    
    % For subsequent white bars after first black bar
    if ( sum( window ) >= 2 ) && blackBars && blackBars < 15 
        while sum( window ) >= 2    % checks and counts white bits 1-by-1
            whitePixels(j) = whitePixels(j) + 1;
            
            c = c + 1;
            window = midSection(:,c);
        end
        whiteBars = whiteBars + 1;
        j = j + 1;
    end
    c = c + 1;
end
c = c - 2; % readjusting column no. to end of the first black guard bar

% Determining coordinates of barcode center-point
barWidth = min(blackPixels, whitePixels); % finding no. of pixels in one bar
barWidth = mode(barWidth)+1; % using min & mode to find the average no.
x_point = c + 0.5*barWidth; % centerpoint = 1st black guard bar + half a bar

clear i j;
%% === STEP 5 - SCANNING FOR T-VALUES ===
cR = c + 3*barWidth; cL = c - 2*barWidth; % column index for right/left-side scan

whiteBars = 0; whitePixels = zeros(2,1);    % Pre-allocation
blackBars = 0; blackPixels = zeros(2,1);

T = zeros(12,1); T1 = T; T2 = T; T3 = T; T4 = T;    % Pre-allocation

for i = 0:5
    % Right-side Scanning
    for n = 1:2     % scanning 2 black and 2 white bars (1 digit equivalent)
        window = midSection(:,cR+1);    % extracting one column
        
        % Checking and correcting cR's position if necessary
        while (sum( window ) == 3) && blackBars == 0 % if first bar is not black
            cR = cR+1;
            window = midSection(:,cR+1);
        end
        
        if sum( window ) ~= 3           % if it's a black bar
            while sum( window ) < 2     % counting the black bits
                blackPixels(n) = blackPixels(n) + 1;

                cR = cR + 1;
                window = midSection(:,cR+1);
            end
            blackBars = blackBars + 1;
        end

        if sum( window ) >= 2           % if it's a white bar
            while sum( window ) >= 2    % counting the white bits
                whitePixels(n) = whitePixels(n) + 1;

                cR = cR + 1;
                window = midSection(:,cR+1);
            end
            whiteBars = whiteBars + 1;
        end
    end
    % Finding T values for the right side, Tn(7:12) 
    T(i+7) = sum( blackPixels + whitePixels ) / barWidth;
    
    T1(i+7) = round( ( blackPixels(1) + whitePixels(1) ) / barWidth );
    T2(i+7) = round( ( whitePixels(1) + blackPixels(2) ) / barWidth );
    T3(i+7) = round( ( blackPixels(2) + whitePixels(2) ) / barWidth );
    T4(i+7) = round( whitePixels(2) / barWidth );
    
    blackPixels = [0; 0]; whitePixels = [0; 0]; % resetting for re-use
    
    % Left-side Scanning
    for n = 1:2
        window = midSection(:,cL);
        
        if sum( window ) ~= 3           % if it's a black bar
            while sum( window ) < 2     % counting the black bits
                blackPixels(n) = blackPixels(n) + 1;

                cL = cL - 1;
                window = midSection(:,cL);
            end
            blackBars = blackBars + 1;
        end

        if sum( window ) >= 2           % if it's a white bar 
            while sum( window ) >= 2    % counting the white bits
                whitePixels(n) = whitePixels(n) + 1;

                cL = cL - 1;
                window = midSection(:,cL);
            end
            whiteBars = whiteBars + 1;
        end
    end
    % Finding T values for the left side, Tn(1:6) 
    T(6-i) = sum( blackPixels + whitePixels ) / barWidth;
    
    T1(6-i) = round( ( blackPixels(1) + whitePixels(1) ) / barWidth );
    T2(6-i) = round( ( whitePixels(1) + blackPixels(2) ) / barWidth );
    T3(6-i) = round( ( blackPixels(2) + whitePixels(2) ) / barWidth );
    T4(6-i) = round( whitePixels(2) / barWidth );
    
    blackPixels = [0; 0]; whitePixels = [0; 0]; % resetting for re-use
end
clear blackPixels whitePixels window;

% T1
T1 = T1./T;
T1( T1 <= (2.5/7) ) = 2; T1( T1 <= (3.5/7) ) = 3; 
T1( T1 <= (4.5/7) ) = 4; T1( T1 < 1 ) = 5;

% T2
T2 = T2./T;
T2( T2 <= (2.5/7) ) = 2; T2( T2 <= (3.5/7) ) = 3; 
T2( T2 <= (4.5/7) ) = 4; T2( T2 < 1 ) = 5;

% T4
T4 = T4./T;
T4( T4 <= (1.5/7) ) = 1; T4( T4 <= (2.5/7) ) = 2;
T4( T4 <= (3.5/7) ) = 3; T4( T4 <= (4.5/7) ) = 4; T4( T4 < 1 ) = 5;

ISBN = zeros(12,1); parityDigits = zeros(12,1); % Pre-allocation

%% === STEP 6 - DECODING BARCODE NUMBER === 

% Decoding with given T-table
for n = 1:12    % for each digit
    switch(1)           % check which condition is true
       case all( [T1(n) T2(n)] == [ 2 2 ] )
          number = 6; parity = 1;
       case all( [T1(n) T2(n)] == [ 2 3 ] )
          number = 0; parity = 0;       
       case all( [T1(n) T2(n)] == [ 2 4 ] )
          number = 4; parity = 1;       
       case all( [T1(n) T2(n)] == [ 2 5 ] )
          number = 3; parity = 0;
          
       case all( [T1(n) T2(n)] == [ 3 2 ] )
          number = 9; parity = 0;
       case all( [T1(n) T2(n)] == [ 3 3 ] ) && T4(n) == 2
          number = 2; parity = 1;  
       case all( [T1(n) T2(n)] == [ 3 3 ] ) && T4(n) == 3
          number = 8; parity = 1;    
       case all( [T1(n) T2(n)] == [ 3 4 ] ) && T4(n) == 2
          number = 1; parity = 0;  
       case all( [T1(n) T2(n)] == [ 3 4 ] ) && T4(n) == 1
          number = 7; parity = 0;   
       case all( [T1(n) T2(n)] == [ 3 5 ] )
          number = 5; parity = 1; 
          
       case all( [T1(n) T2(n)] == [ 4 2 ] )
          number = 9; parity = 1;  
       case all( [T1(n) T2(n)] == [ 4 3 ] ) && T4(n) == 2
          number = 2; parity = 0;
       case all( [T1(n) T2(n)] == [ 4 3 ] ) && T4(n) == 1
          number = 8; parity = 0;
       case all( [T1(n) T2(n)] == [ 4 4 ] ) && T4(n) == 1
          number = 1; parity = 1;
       case all( [T1(n) T2(n)] == [ 4 4 ] ) && T4(n) == 2
          number = 7; parity = 1;
       case all( [T1(n) T2(n)] == [ 4 5 ] )
          number = 5; parity = 0;
          
       case all( [T1(n) T2(n)] == [ 5 2 ] )
          number = 6; parity = 0;
       case all( [T1(n) T2(n)] == [ 5 3 ] )
          number = 0; parity = 1;
       case all( [T1(n) T2(n)] == [ 5 4 ] )
          number = 4; parity = 0; 
       case all( [T1(n) T2(n)] == [ 5 5 ] )
          number = 3; parity = 1;
    end
    
    ISBN(n) = number;           % decoded ISBN without country code
    parityDigits(n) = parity;   % parity digits for country code decoding
end
clear T T1 T2 T3 T4 i n parity;

% Plotting center-point and confirming the points where the decoding stops
figure(1)
imshow(barCodeImg);
axis on
hold on;
plot( x_point, rowMid , 'r+','MarkerSize', 10,'LineWidth',2 );
plot( cL, rowMid , 'g+','MarkerSize', 10,'LineWidth',1.5 );
plot( cR, rowMid , 'g+','MarkerSize', 10,'LineWidth',1.5 );

%% === STEP 7 - FINDING COUNTRY CODE ===

% Decoding based on given parity values, E = 1 and O = 0
switch (1)          % check which condition is true
    case isequal( parityDigits(1:6), ~parityDigits(7:12) ) % if OOOOOOEEEEEE
        % leave it as it is because it is UPC
    case isequal( parityDigits(1:6), [ 0 0 0 0 0 0 ]' )
        ISBN = [0; ISBN(:)];    % insert country code into ISBN
    case isequal( parityDigits(1:6), [ 0 0 1 0 1 1 ]' )
        ISBN = [1; ISBN(:)];
    case isequal( parityDigits(1:6), [ 0 0 1 1 0 1 ]' )
        ISBN = [2; ISBN(:)];
    case isequal( parityDigits(1:6), [ 0 0 1 1 1 0 ]' )
        ISBN = [3; ISBN(:)];
    case isequal( parityDigits(1:6), [ 0 1 0 0 1 1 ]' )
        ISBN = [4; ISBN(:)];
    case isequal( parityDigits(1:6), [ 0 1 1 0 0 1 ]' )
        ISBN = [5; ISBN(:)];
    case isequal( parityDigits(1:6), [ 0 1 1 1 0 0 ]' )
        ISBN = [6; ISBN(:)];
    case isequal( parityDigits(1:6), [ 0 1 0 1 0 1 ]' )
        ISBN = [7; ISBN(:)];
    case isequal( parityDigits(1:6), [ 0 1 0 1 1 0 ]' )
        ISBN = [8; ISBN(:)];
    case isequal( parityDigits(1:6), [ 0 1 1 0 1 0 ]' )
        ISBN = [9; ISBN(:)];
end

% Displaying scanned code
fprintf('The scanned barcode is ');
fprintf('%d',ISBN);
fprintf('\n');

%% === STEP 8 - MODULO CHECK ===
if length(ISBN) == 13       % For EAN-13 barcodes
    identity = [1 3 1 3 1 3 1 3 1 3 1 3]'; % even x 3, odd x 1
else                        % For UPC barcodes
    identity = [1 3 1 3 1 3 1 3 1 3 1]';
end

checksum = sum( identity.*ISBN(1:end-1) );
modulo = 10 - mod(checksum,10); % checks if barcode was scanned correctly

% Printing outcome
if modulo == ISBN(end)
    fprintf('Modulo is %d; Scanned barcode is correct.\n', modulo);
else                        
    fprintf('Modulo is %d. Scanned barcode is not correct.\n', modulo);
end

% print('-bestfit','barcode_3','-dpdf')
% ---------------------------
% END OF CODE