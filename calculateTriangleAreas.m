function A = calculateTriangleAreas(vertices, triangles)
    % vertices: Mx3 array of vertex points
    % triangles: Nx3 array of triangles (indices into the vertices array)
    % areas: Nx1 array of triangle areas
    
    A = zeros(size(triangles, 1), 1);
    for i = 1:size(triangles, 1)
        idx = triangles(i, :);
        v1 = vertices(idx(1), :);
        v2 = vertices(idx(2), :);
        v3 = vertices(idx(3), :);
        A(i) = 0.5 * norm(cross(v2 - v1, v3 - v1));
    end
end