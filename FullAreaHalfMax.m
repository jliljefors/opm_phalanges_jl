function FAHM = FullAreaHalfMax(sourcedistribution,latency,i_phalange)
%UNTITLED Calculates the full area at half max amplitude
%   Detailed explanation goes here
[~,i_latency] = min(abs(sourcedistribution.time-latency));

% Half max level
half_max = max(sourcedistribution.avg{i_phalange}.pow(:,i_latency))/2;

% Find triangles that have at least one point with amplitude >= half max
i_vertices = find(sourcedistribution.avg{i_phalange}.pow(:,i_latency)>=half_max);
[triangles,~] = find(ismember(sourcedistribution.tri,i_vertices)); 
triangles = sourcedistribution.tri(triangles,:);

% Sum area of triangles and divide by 3 (since its a triangle per point). 
FAHM = sum(calculateTriangleAreas(sourcedistribution.pos, triangles))/3;
end