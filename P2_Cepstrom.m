[X,F_s] = audioread('C:\Users\Hosein\Desktop\voice.wav','native');
X = double(X);
N=40;
frame_size = N*F_s/1000;
X_size =  length(X);
overlap = 0.5;
num_of_window = ceil( (X_size - frame_size)/(frame_size*overlap) ) -1;
rect_windows = rectwin(frame_size);

candidate_3lvl = cepstrom_analysis(X,rect_windows,frame_size,N,overlap,num_of_window,1);
candidate_cp = cepstrom_analysis(X,rect_windows,frame_size,N,overlap,num_of_window,2);

[ pX_3lvl ,pY_3lvl ] = calcute_fc(candidate_3lvl,num_of_window);
[pointX_3lvl , pointY_3lvl ] = all_points(candidate_3lvl,F_s);
[ pX_cp ,pY_cp ] = calcute_fc(candidate_cp,num_of_window);
[pointX_cp , pointY_cp ] = all_points(candidate_cp,F_s);

subplot(4,1,4)
plot(pointX_3lvl,pointY_3lvl,'.r')
xlim([0 num_of_window])
ylabel("Pitch Candidates"+sprintf('\n')+"3level Center Clipping (Hz)",'Fontsize',9);

subplot(4,1,3)
plot(pX_3lvl,pY_3lvl,'LineWidth',1.5)
ylabel("Pitch Frequency"+sprintf('\n')+"3level Center Clipping (Hz)",'Fontsize',9);
xlim([0 num_of_window])
ylim([50,250])

subplot(4,1,2)
plot(pointX_cp,pointY_cp,'.g')
xlim([0 num_of_window])
ylabel("Pitch Candidates"+sprintf('\n')+"Center Clipping (Hz)",'Fontsize',9);

subplot(4,1,1)
plot(pX_cp,pY_cp,'LineWidth',1.5)
ylabel("Pitch Frequency"+sprintf('\n')+"Center Clipping (Hz)",'Fontsize',9);
xlim([0 num_of_window])
ylim([50,250])
title("Cepstral Analysis , Frame Length="+string(N)+"ms , Number Of Windows="+string(num_of_window));

function candidate = cepstrom_analysis(X,rect_windows,frame_size,N,overlap,num_of_window,type)
    candidate = zeros(num_of_window,5);
    for i=1:num_of_window
        y= repmat(rect_windows,1,1) .* X((i-1)*overlap*frame_size+1:(i-1)*overlap*frame_size + frame_size);
        e = energy(y);
        if type==1
            y = center_cliping(y);
        else
           y =  center_cliping_3level(y);
        end
        if e < 1000000
            continue
        else
            cepstrom=real(ifft(log(abs(fft(y)))));
            n_ceps=length(cepstrom);
            cepstrom=cepstrom(1:n_ceps/2);
            [ pks , locs ] = findpeaks(cepstrom);
            ter = 5;
            for j=1:ter
                [ maximum , max_index ] = max(pks);
                pks = pks(pks~=maximum);
                if  length(max_index)>1
                    candidate(i,j) = locs(max_index(1));
                    candidate(i,j+1) = locs(max_index(2));
                    locs = locs(locs~=locs(max_index(1)));
                    locs = locs(locs~=locs(max_index(1)));
                    j = j+1;
                elseif length(max_index)==1
                    candidate(i,j) = locs(max_index);
                    locs = locs(locs~=locs(max_index));
                end
            end      
        end   
    end
end


function [ pointX , pointY ] = all_points(array,F_s)
    pointX = NaN(1,5*length(array));
    pointY = NaN(1,5*length(array));
    for i=1:length(array)
       for j=1:5
           if array(i,j)>75 
               pointX((i-1)*5+j) = i;
               pointY((i-1)*5+j) = F_s / array(i,j);
           end
           
       end
    end
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

function [ pointsX , pointsY ] = calcute_fc(final_candidate,num_of_window)
    pointsX = NaN(1,num_of_window);
    pointsY = NaN(1,num_of_window);
    for i=1:num_of_window
       for j=1:5
          if  final_candidate(i,j) < 75 ||  final_candidate(i,j) >200 
              final_candidate(i,j)  = nan;
          end
       end
    end
    
    final_candidate = smooth2a(final_candidate-15,20,20);
    
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
