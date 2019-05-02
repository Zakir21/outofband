function [B,b] = Freq_Diff_BF_Sph(geometry, SIGNALS, rec_weights)

lowFID = geometry.lowFID; % Low end of bandwidth
highFID = geometry.highFID; % High end of bandwidth
df_i = geometry.df_i; % Difference frequency (number of bins)
df = geometry.df; % Difference frequency
xRange = geometry.xRange; % Test coordinates (x)
yRange = geometry.yRange; % Test coordinates (y)
recArray = geometry.recArray; % Array geometry (x,y columns)

ii = 0;
for fi = lowFID:highFID-df_i
    ii = ii + 1;
    for n = 1:numRec
        RangeTemp = sqrt((xRange-recArray(n,1)).^2+(yRange-recArray(n,2)).^2);
        tau = RangeTemp/c;
        
        % 1/r is optional
        A = exp(1i*2*pi*df*tau)./RangeTemp;
        
        % Option of normalizing outputs by frequency.
%         y(n) = rec_weights(n)*conj(SIGNALS(n,fi))*SIGNALS(n,fi+df_i)/( norm(SIGNALS(:,fi)/numRec)*norm(SIGNALS(:,fi+df_i)) );
        y(n) = rec_weights(n)*conj(SIGNALS(n,fi))*SIGNALS(n,fi+df_i);
        
        xhat(:,:,n) = A*y(n);
    end
    b(:,:,ii) = sum(xhat,3);
    B(:,:,ii) = abs( sum(xhat,3) ).^2;
end
