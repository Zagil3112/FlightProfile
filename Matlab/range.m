function R = range(lat,long)

s(1)=0;
n=length(lat)-1;

for i=1:n
	lat1=lat(i); % Latitude 1
	lat2=lat(i+1); % Latitude 2
	long1=long(i); % Longitude 1
	long2=long(i+1); % Longitude 2
	s(i+1)=distVincenty(lat1,lat2,long1,long2);
end

