local CEIL = math.ceil

local function draw_flat_top_triangle(v0,v1,v2,points,caller)
    local v0x,v0y,v0z = v0.x,v0.y,v0.z
    local v1x,v1y,v1z = v1.x,v1.y,v1.z
    local v2x,v2y,v2z = v2.x,v2.y,v2.z
    local LIST,LUT = points.list,points.LUT
    local m0 = (v2x - v0x) / (v2y - v0y)
    local m1 = (v2x - v1x) / (v2y - v1y)
    local y_start = CEIL(v0y - 0.5)
    local y_end   = CEIL(v2y - 0.5) - 1
    local n = #LIST
    for y=y_start,y_end do 
        local px0 = m0 * (y + 0.5 - v0y) + v0x
        local px1 = m1 * (y + 0.5 - v1y) + v1x
        local x_start = CEIL(px0 - 0.5)
        local x_end   = CEIL(px1 - 0.5)
        for x=x_start,x_end do
            n = n + 1
            LIST[n] =  {x=x,y=y}
            if not LUT[y] then LUT[y] = {} end
            LUT[y][x] = true
            if caller then caller(x,y) end
        end
    end
end

local function draw_flat_bottom_triangle(v0,v1,v2,points,caller)
    local v0x,v0y,v0z = v0.x,v0.y,v0.z
    local v1x,v1y,v1z = v1.x,v1.y,v1.z
    local v2x,v2y,v2z = v2.x,v2.y,v2.z
    local LIST,LUT = points.list,points.LUT
    local m0 = (v1x - v0x) / (v1y - v0y)
    local m1 = (v2x - v0x) / (v2y - v0y)
    local y_start = CEIL(v0y - 0.5)
    local y_end   = CEIL(v2y - 0.5) - 1
    local n = #LIST
    for y=y_start,y_end do 
        local px0 = m0 * (y + 0.5 - v0y) + v0x
        local px1 = m1 * (y + 0.5 - v0y) + v0x
        local x_start = CEIL(px0 - 0.5)
        local x_end   = CEIL(px1 - 0.5)
        for x=x_start,x_end do
            n = n + 1
            LIST[n] =  {x=x,y=y}
            if not LUT[y] then LUT[y] = {} end
            LUT[y][x] = true
            if caller then caller(x,y) end
        end
    end
end

local function get_triangle_points(vector0,vector1,vector2,caller)
    local points = {list={},LUT={}}
    if vector1.y < vector0.y then vector0,vector1 = vector1,vector0 end
    if vector2.y < vector1.y then vector1,vector2 = vector2,vector1 end
    if vector1.y < vector0.y then vector0,vector1 = vector1,vector0 end
    if vector0.y == vector1.y then
        if vector1.x < vector0.x then vector0,vector1 = vector1,vector0 end
        draw_flat_top_triangle(vector0,vector1,vector2,points,caller)
    elseif vector1.y == vector2.y then
        if vector2.x < vector1.x then vector1,vector2 = vector2,vector1 end
        draw_flat_bottom_triangle(vector0,vector1,vector2,points,caller)
    else
        local alpha_split = (vector1.y-vector0.y) / (vector2.y-vector0.y)
        local split_vertex = { 
            x = vector0.x + ((vector2.x - vector0.x) * alpha_split),      
            y = vector0.y + ((vector2.y - vector0.y) * alpha_split),
        }
        if vector1.x < split_vertex.x then
            draw_flat_bottom_triangle(vector0,vector1,split_vertex,points,caller)
            draw_flat_top_triangle   (vector1,split_vertex,vector2,points,caller)
        else
            draw_flat_bottom_triangle(vector0,split_vertex,vector1,points,caller)
            draw_flat_top_triangle   (split_vertex,vector1,vector2,points,caller)
        end
    end
    return points.list
end

return {
    create=get_triangle_points
}