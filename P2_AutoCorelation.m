[X,F_s] = audioread('C:\Users\Hosein\Desktop\voice.wav','native');
X = double(X);
N=40;
frame_size = N*F_s/1000;
X_size =  length(X);
overlap = 0.5;
num_of_window = ceil( (X_size - frame_size)/(frame_size*overlap) ) -1;



candidate_3lvl_cp = find_pitch_freq(X,1,num_of_window,frame_size,overlap,N);
candidate_cp = find_pitch_freq(X,2,num_of_window,frame_size,overlap,N);
[ pointsX_3level , pointsY_3level ] = calcute_fc(candidate_3lvl_cp,num_of_window,1);
[ pointsX_cp , pointsY_cp ] = calcute_fc(candidate_cp,num_of_window,0);
[all_pointX_cp , all_pointY_cp]= all_points(candidate_cp);
[all_pointX_3lvl , all_pointY_3lvl]= all_points(candidate_3lvl_cp);

%3 level center Cliping
subplot(4,1,4)
plot(all_pointX_3lvl,all_pointY_3lvl,'.','MarkerEdgeColor','r')
ylabel("Pitch Candidate"+sprintf('\n')+"3level Center Clipping (Hz)",'Fontsize',10);
xlim([0 num_of_window])
ylim([50 200])


%pitch with 3level 
subplot(4,1,3)
plot(pointsX_3level,pointsY_3level,'b','LineWidth',1.5)
ylabel("Pitch Frequency"+sprintf('\n')+"3level Center Clipping (Hz)",'Fontsize',9);
xlim([0 num_of_window])
ylim([50 200])

%center Cliping 
subplot(4,1,2)
plot(all_pointX_cp,all_pointY_cp,'.','MarkerEdgeColor','g')
ylabel("Pitch Candidate"+sprintf('\n')+"Center Clipping (Hz)",'Fontsize',10);
xlim([0 num_of_window])
ylim([50 200])

%pitch with cp
subplot(4,1,1)
plot(pointsX_cp,pointsY_cp,'b','LineWidth',1.5)
ylabel("Pitch Frequency"+sprintf('\n')+"Center Clipping(Hz)",'Fontsize',9);
xlim([0 num_of_window])
ylim([50 200])
title("AutoCorrelation Function , Frame Length="+string(N)+"ms , Number Of Windows="+string(num_of_window));



function [ pointX , pointY ] = all_points(array)
    pointX = NaN(1,5*length(array));
    pointY = NaN(1,5*length(array));
    for i=1:length(array)
       for j=1:5
           if array(i,j)>75 && array(i,j)<200
               pointX((i-1)*5+j) = i;
               pointY((i-1)*5+j) = array(i,j);
           end
           
       end
    end
end

function [ pointsX , pointsY ] = calcute_fc(final_candidate,num_of_window,islvl)
    pointsX = NaN(1,num_of_window);
    pointsY = NaN(1,num_of_window);
    for i=1:num_of_window
       for j=1:5
          if  final_candidate(i,j) < 75 ||  final_candidate(i,j) >200 
              final_candidate(i,j)  = nan;
          end
       end
    end
    
    final_candidate = smooth2a(final_candidate,13-islvl*2,7-islvl*5);
    
    for i=1:num_of_window
        pointsX(i) =  i ;
    end
    
    for i=1:num_of_window
        for j=1:5
            if final_candidate(i,j) > 75 && final_candidate(i,j) <200
                pointsY(i) = final_candidate(i,j); 
                break;
            end
        end   
    end
end

function candidate = find_pitch_freq(X,type_of_centerclipping,num_of_window,frame_size,overlap,N)
    rect_windows = rectwin(frame_size);
    acf = [];
    Energy = [];
    for i=1:num_of_window
        y= repmat(rect_windows,1,1) .* X((i-1)*overlap*frame_size+1:(i-1)*overlap*frame_size + frame_size);
        e=energy(y);
        Energy = [ Energy , e];
        z=zero_crossing(y);
        if type_of_centerclipping==1
            c_c = center_cliping_3level(y);   
        elseif type_of_centerclipping==2
            c_c = center_cliping(y);
        end
        ac = autocorelation(c_c);
        if e < 500000
            acf = [acf ;  zeros(1,frame_size) ];
        else
            acf = [acf ; ac ]; 
        end
    end

    for i=1:num_of_window

        locs = find_locs(acf(i,:));

        for k=1:length(locs)
            candidate(i,k) = 1/((locs(k))*N/(frame_size*1000));  
        end
        if length(locs)<=1
           candidate(i,1)=0; 
           candidate(i,2)=0;
           candidate(i,3)=0;
           candidate(i,4)=0;
           candidate(i,5)=0;
        end
    end
end



function locs = find_locs(y)
    p = smooth(smooth(smooth(y,10),10)) ;
    [pks,locs] = findpeaks(p);
    new_locs = [];
    minimum = min(5,length(locs));
    
    for j=1:minimum
        maximum = max(pks);
        index = find(pks==maximum);
        if length(index)==1
            new_locs = [ new_locs , locs(index) ];
            pks = pks(pks~=maximum);
            locs = locs(locs~=locs(index)); 
        elseif length(index)==2
            x = locs(index);
            new_locs = [ new_locs , x(1) ,x(2)  ];
            pks = pks(pks~=maximum);
            locs = locs(locs~=x(1));
            locs = locs(locs~=x(2));
        end   
    end
    locs = new_locs;   
end


function output = center_cliping(y)
    output = zeros(1,length(y));
    maxvalue = max(abs(y));
    coef = 0.2;
    for i=1:length(y)
       if y(i) >= coef* maxvalue
           output(1,i) = y(i) -  coef* maxvalue;
       elseif y(i) <= -coef* maxvalue
           output(1,i) = y(i) +  coef* maxvalue;
       end
    end
end


function output = center_cliping_3level(y)
    output = zeros(1,length(y));
    maxvalue = max(abs(y));
    coef = 0.2;
    for i=1:length(y)
       if y(i) >= coef* maxvalue
           output(1,i) = 1;
       elseif y(i) <= -coef* maxvalue
           output(1,i) = -1;
       end
    end
end

function e = energy(y)
    e=0;
    for j=1:length(y)
        e = e + y(j)*y(j);
    end
end

function z = zero_crossing(y)
    z=0;
    for j=2:length(y)
        z = z + abs( sign( y(j) ) - sign( y(j-1) ))/2 ;
    end
    z = z/(2*length(y));
end

function auto_cor = autocorelation(y)
    auto_cor = zeros(1,length(y));
    for k=0:length(y)-1
        sum=0;
        for j=k:length(y)-1
            sum =  sum + y(j+1)*y(j+1-k);
        end
        auto_cor(k+1) = sum;
    end
end

function matrixOut = smooth2a(matrixIn,Nr,Nc)
    if nargin < 2, error('Not enough input arguments!'), end
    N(1) = Nr; 
    if nargin < 3, N(2) = N(1); else N(2) = Nc; end
    if length(N(1)) ~= 1, error('Nr must be a scalar!'), end
    if length(N(2)) ~= 1, error('Nc must be a scalar!'), end
    [row,col] = size(matrixIn);
    eL = spdiags(ones(row,2*N(1)+1),(-N(1):N(1)),row,row);
    eR = spdiags(ones(col,2*N(2)+1),(-N(2):N(2)),col,col);
    A = isnan(matrixIn);
    matrixIn(A) = 0;
    nrmlize = eL*(~A)*eR;
    nrmlize(A) = NaN;
    matrixOut = eL*matrixIn*eR;
    matrixOut = matrixOut./nrmlize;
end
