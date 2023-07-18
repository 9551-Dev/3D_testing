local s = 0.5

local mat = require("mat")

local function u(u)
    return (u+1)/3
end
local function v(v)
    return v/4
end

local objects = {
    {
        verts={
            -s,-s,-s, u(1) ,v(0),
            s,-s,-s,  u(0) ,v(0),
            -s,s,-s,  u(1) ,v(1),
            s,s,-s,   u(0) ,v(1),
            -s,-s,s,  u(1) ,v(3),
            s,-s,s,   u(0) ,v(3),
            -s,s,s,   u(1) ,v(2),
            s,s,s,    u(0) ,v(2),
            -s,-s,-s, u(1) ,v(4),
            s,-s,-s,  u(0) ,v(4),
            -s,-s,-s, u(2) ,v(1),
            -s,-s,s,  u(2) ,v(2),
            s,-s,-s,  u(-1),v(1),
            s,-s,s,   u(-1),v(2),
        },
        indices={
            1,3,2,   3,4,2,
            5,9,6,   6,9,10,
            3,7,4,   4,7,8,
            5,6,8,   5,8,7,
            3,11,12, 3,12,7,
            13,4,8,  13,8,14
        }
    },
}

--[[local objects = {
    {
        verts={
            -s,-s,-s, 0,1,
            s,-s,-s,  1,1,
            -s,s,-s,  0,0,
            s,s,-s,   1,0,
            -s,-s,s,  1,1,
            s,-s,s,   0,1,
            -s,s,s,   1,0,
            s,s,s,    0,0
        },
        indices={
            1,3,2, 3,4,2,
            2,4,6, 4,8,6,
            3,7,4, 4,7,8,
            5,6,8, 5,8,7,
            1,5,3, 3,5,7,
            1,2,5, 2,6,5
        }
    },
}]]

local res = 5

local worg,horg = love.graphics.getDimensions()
love.window.setMode(worg,horg,{resizable=true})
local w = worg / res
local h = horg / res

local const = 10

local function getT(a,b,p)
    local v1 = {a[1]-b[1], a[2]-b[2]}
    local v2 = {a[1]-p[1], a[2]-p[2]}

    return (v1[1]*v2[1] + v1[2]*v2[2]) / (v1[1]*v1[1] + v1[2]*v1[2])
end

local function intY(a,b,y)
    local oy = (b[2] - a[2])
    if oy == 0 then oy = 5e-50 end

    return a[1] + ((y - a[2]) * (b[1] - a[1])) / oy
end

local function transform_screen_space(p,w,h)
    local inverse_z = 1/p[4]
    return {
        ( p[1]*inverse_z+1)*w/2,
        (-p[2]*inverse_z+1)*h/2,
        inverse_z,
        p[4],
        p[5]*inverse_z,
        p[6]*inverse_z
    }
end

local function cull(a,b,c)
    local i1,i2,i3 = a[1],a[2],a[4]

    local a1 = b[1]-i1
    local a2 = b[2]-i2
    local a3 = b[4]-i3

    local b1 = c[1]-i1
    local b2 = c[2]-i2
    local b3 = c[4]-i3

    return (a2*b3 - a3*b2)*i1+
        (a3*b1 - a1*b3)   *i2+
        (a1*b2 - a2*b1)   *i3
end

local function makePerspective(w, h, nearPlane, farPlane, fov)
    local aspectRatio = w/h
    local fov_rad = math.rad(fov)
    local tan_half_fov = math.tan(fov_rad * 0.5)

    local A = 1 / (aspectRatio * tan_half_fov)
    local B = 1 / tan_half_fov
    local C = (farPlane + nearPlane) / (farPlane - nearPlane)
    local D = 1
    local E = -(2 * farPlane * nearPlane) / (farPlane - nearPlane)

    return mat.new(4, 4,
        A, 0, 0, 0,
        0, B, 0, 0,
        0, 0, C, E,
        0, 0, D, 0
    )
end

local sqrt = math.sqrt

local function cross(x1,y1,z1,x2,y2,z2)
    return y1*z2-z1*y2,
        z1*x2-x1*z2,
        x1*y2-y1*x2
end

local function dot(x1,y1,z1,x2,y2,z2)
    return x1*x2 + y1*y2 + z1*z2
end

local function normalize(x,y,z)
    local len = sqrt(x*x+y*y+z*z)
    return x/len,y/len,z/len
end

local function normalize_vector(v)
    local x,y,z = normalize(v.x,v.y,v.z)
    return {
        x=x,
        y=y,
        z=z
    }
end

local function get_look_vector(yaw,pitch)
    return normalize_vector{
        x = math.sin(yaw)*math.cos(pitch),
        y = -math.sin(pitch),
        z = math.cos(yaw)*math.cos(pitch)
    }
end

local function get_move_vector(yaw)
    return normalize_vector{
        x = math.sin(yaw),
        y = 0,
        z = math.cos(yaw)
    }
end

local function new_vector(x,y,z)
    return {
        x = x or 0,
        y = y or 0,
        z = z or 0
    }
end

local function add_vector(a,b)
    return {
        x = a.x + b.x,
        y = a.y + b.y,
        z = a.z + b.z
    }
end

local function sub_vector(a,b)
    return {
        x = a.x - b.x,
        y = a.y - b.y,
        z = a.z - b.z
    }
end

local function scale_vector(a,scalar)
    return {
        x = a.x * scalar,
        y = a.y * scalar,
        z = a.z * scalar
    }
end

local function makeLookat(fromx,fromy,fromz,atx,aty,atz)
    local nx,ny,nz = normalize(atx-fromx,aty-fromy,atz-fromz)
    local ux,uy,uz = normalize(cross(0,1,0,nx,ny,nz))
    local vx,vy,vz = normalize(cross(nx,ny,nz,ux,uy,uz))

    return mat.new(4,4,
        ux,uy,uz,-dot(fromx,fromy,fromz,ux,uy,uz),
        vx,vy,vz,-dot(fromx,fromy,fromz,vx,vy,vz),
        nx,ny,nz,-dot(fromx,fromy,fromz,nx,ny,nz),
        0,0,0,1
    )
end

local RAD,SIN,COS = math.rad,math.sin,math.cos
local function create_rotation_eulers(eulers)
    local x = RAD(eulers.x)
    local y = RAD(eulers.y)
    local z = RAD(eulers.z)
    local sx = SIN(x)
    local sy = SIN(y)
    local sz = SIN(z)

    local cx = COS(x)
    local cy = COS(y)
    local cz = COS(z)

    return mat.new(4,4,
        cy*cz,-cy*sz,sy,0,
        (sx*sy*cz) +(cx*sz), (-sx*sy*sz)+(cx*cz),-sx*cy,0,
        (-cx*sy*cz)+(sx*sz), (cx*sy*sz) +(sx*cz),cx*cy,0,
        0,0,0,1
    )
end

local per

local CEIL,MAX,MIN,FLOOR = math.ceil,math.max,math.min,math.floor

local function lerp(v1,v2,t)
    return (1 - t) * v1 + t * v2
end

local function get_bary_coord(x,y,p1,p2,p3)
    local div = ((p2[2]-p3[2])*(p1[1]-p3[1]) + (p3[1]-p2[1])*(p1[2]-p3[2]))
    local ba = ((p2[2]-p3[2])*(x-p3[1]) + (p3[1]-p2[1])*(y-p3[2])) / div
    local bb = ((p3[2]-p1[2])*(x-p3[1]) + (p1[1]-p3[1])*(y-p3[2])) / div

    return {ba,bb,1-ba-bb}
end

local function interpolate_uv(bary_a,bary_b,bary_c,uv1,uv2,uv3)
    return {
        uv1[5] * bary_a + uv2[5] * bary_b + uv3[5] * bary_c,
        uv1[6] * bary_a + uv2[6] * bary_b + uv3[6] * bary_c
    }
end

local function draw_flat_top_triangle(v0,v1,v2,tex,origin,fragment)
    local v0x,v0y = v0[1],v0[2]
    local v1x,v1y = v1[1],v1[2]
    local v2x,v2y = v2[1],v2[2]
    local m0 = (v2x - v0x) / (v2y - v0y)
    local m1 = (v2x - v1x) / (v2y - v1y)
    local y_start = CEIL(v0y - 0.5)
    local y_end   = CEIL(v2y - 0.5) - 1

    for y=y_start,y_end do
        local px0 = m0 * (y + 0.5 - v0y) + v0x
        local px1 = m1 * (y + 0.5 - v1y) + v1x
        local x_start = CEIL(px0 - 0.5)
        local x_end   = CEIL(px1 - 0.5)

        local sx_start = intY(origin[1],origin[2],y)
        local sx_end   = intY(origin[1],origin[3],y)
        local t1 = getT(origin[1],origin[2],{sx_start,y})
        local t2 = getT(origin[1],origin[3],{sx_end,y})
        local w1 = lerp(origin[1][3],origin[2][3],t1)
        local w2 = lerp(origin[1][3],origin[3][3],t2)
        local z1 = lerp(origin[1][4],origin[2][4],t1)
        local z2 = lerp(origin[1][4],origin[3][4],t2)

        for x=x_start,x_end do
            local bary = get_bary_coord(x,y,v0,v1,v2)
            local tpos = interpolate_uv(bary,v0,v1,v2)

            local div = sx_end - sx_start
            local t3 = (x - sx_start) / ((div == 0) and 1 or div)


            local z = 1/lerp(w1,w2,t3)

            fragment(x,y,lerp(z1,z2,t3),
                MAX(0,MIN(CEIL(tpos[1]*z*(tex.w+res)),tex.w)-1),
                MAX(0,MIN(CEIL(tpos[2]*z*(tex.h+res)),tex.h)-1),{t3=origin[1],r=bary[1],g=bary[2],b=bary[3]}
            )
        end
    end
end

local function draw_flat_bottom_triangle(v0,v1,v2,tex,origin,fragment)
    local v0x,v0y = v0[1],v0[2]
    local v1x,v1y = v1[1],v1[2]
    local v2x,v2y = v2[1],v2[2]
    local m0 = (v1x - v0x) / (v1y - v0y)
    local m1 = (v2x - v0x) / (v2y - v0y)
    local y_start = CEIL(v0y - 0.5)
    local y_end   = CEIL(v2y - 0.5) - 1

    for y=y_start,y_end do
        local px0 = m0 * (y + 0.5 - v0y) + v0x
        local px1 = m1 * (y + 0.5 - v0y) + v0x
        local x_start = CEIL(px0 - 0.5)
        local x_end   = CEIL(px1 - 0.5)

        local sx_start = intY(origin[1],origin[2],y)
        local sx_end   = intY(origin[1],origin[3],y)
        local t1 = getT(origin[1],origin[2],{sx_start,y})
        local t2 = getT(origin[1],origin[3],{sx_end,y})
        local w1 = lerp(origin[1][3],origin[2][3],t1)
        local w2 = lerp(origin[1][3],origin[3][3],t2)
        local z1 = lerp(origin[1][4],origin[2][4],t1)
        local z2 = lerp(origin[1][4],origin[3][4],t2)

        for x=x_start,x_end do
            local bary =  get_bary_coord(x,y,v0,v1,v2)
            local tpos =  interpolate_uv(bary,v0,v1,v2)

            local div = sx_end - sx_start
            local t3 = (x - sx_start) / ((div == 0) and 1 or div)

            local z = 1/lerp(w1,w2,t3)

            fragment(x,y,lerp(z1,z2,t3),
                MAX(0,MIN(CEIL(tpos[1]*z*(tex.w+res)),tex.w)-1),
                MAX(0,MIN(CEIL(tpos[2]*z*(tex.h+res)),tex.h)-1),{t3=m0,r=bary[1],g=bary[2],b=bary[3]}
            )
        end
    end
end

local function slope(x1,y1,x2,y2)
    return (y2-y1)/(x2-x1)
end

local function bary(x,y,p1,p2,p3)
    local p23y_delta,p13x_delta,p32x_delta = p2[2]-p3[2],p1[1]-p3[1],p3[1]-p2[1]

    local xp3_delta,yp3_delta = x-p3[1],y-p3[2]

    local div = (p23y_delta*p13x_delta + p32x_delta*(p1[2]-p3[2]))

    local dot_a = (p23y_delta*xp3_delta    + p32x_delta*yp3_delta) / div
    local dot_b = ((p3[2]-p1[2])*xp3_delta + p13x_delta*yp3_delta) / div

    return dot_a,dot_b,1-dot_a-dot_b
end

local function raster_triangle(p1,p2,p3,tex,frag)
    local ori = {p1,p2,p3}
    if p2[2] < p1[2] then p1,p2 = p2,p1 end
    if p3[2] < p2[2] then p2,p3 = p3,p2 end
    if p2[2] < p1[2] then p1,p2 = p2,p1 end
    if p1[2] == p2[2] then
        if p2[1] < p1[1] then p1,p2 = p2,p1 end
        draw_flat_top_triangle(p1,p2,p3,tex,ori,frag)
    elseif p2[2] == p3[2] then
        if p3[1] < p2[1] then p2,p3 = p3,p2 end
        draw_flat_bottom_triangle(p1,p2,p3,tex,ori,frag)
    else
        local alpha_split = (p2[2]-p1[2]) / (p3[2]-p1[2])
        local split_vertex = {
            lerp(p1[1],p3[1],alpha_split),
            lerp(p1[2],p3[2],alpha_split),
            lerp(p1[3],p3[3],alpha_split),
            lerp(p1[4],p3[4],alpha_split),
            lerp(p1[5],p3[5],alpha_split),
            lerp(p1[6],p3[6],alpha_split)
        }

        if p2[1] < split_vertex[1] then
            draw_flat_bottom_triangle(p1,p2,split_vertex,tex,ori,frag)
            draw_flat_top_triangle   (p2,split_vertex,p3,tex,ori,frag)
        else
            draw_flat_bottom_triangle(p1,split_vertex,p2,tex,ori,frag)
            draw_flat_top_triangle   (split_vertex,p2,p3,tex,ori,frag)
        end
    end
end


local function raster_triangle(p1,p2,p3,tex,frag)
    if p1[2] > p3[2] then p1,p3 = p3,p1 end
    if p1[2] > p2[2] then p1,p2 = p2,p1 end
    if p2[2] > p3[2] then p2,p3 = p3,p2 end

    local split_alpha = (p2[2]-p1[2])/(p3[2]-p1[2])

    local split_point =  {
        lerp(p1[1],p3[1],split_alpha),
        lerp(p1[2],p3[2],split_alpha),
        lerp(p1[3],p3[3],split_alpha),
        lerp(p1[4],p3[4],split_alpha),
        lerp(p1[5],p3[5],split_alpha),
        lerp(p1[6],p3[6],split_alpha)
    }

    local left_point,right_point = p2,split_point
    if left_point[1] > right_point[1] then
        left_point,right_point = right_point,left_point
    end

    local delta_left_top  = 1/slope(p1[1],p1[2],left_point[1], left_point[2])
    local delta_right_top = 1/slope(p1[1],p1[2],right_point[1],right_point[2])

    local delta_left_bottom  = 1/slope(p3[1],p3[2],left_point[1], left_point[2])
    local delta_right_bottom = 1/slope(p3[1],p3[2],right_point[1],right_point[2])

    local x_left,x_right = p1[1],p1[1]
    for y=math.ceil(p1[2]),math.ceil(left_point[2])-1 do

        for x=math.floor(x_left-0.5),math.ceil(x_right-0.5) do
            local bary_a,bary_b,bary_c = bary(math.ceil(x),math.ceil(y),p1,left_point,right_point)
            local depth = p1[3]*bary_a+left_point[3]*bary_b+right_point[3]*bary_c

            local z = 1 / depth

            local tpos = interpolate_uv(bary_a,bary_b,bary_c,p1,left_point,right_point)

            frag(x,y,depth,
                MAX(0,MIN(CEIL(tpos[1]*z*(tex.w+res)),tex.w)-1),
                MAX(0,MIN(CEIL(tpos[2]*z*(tex.h+res)),tex.h)-1),{r=1,g=1,b=1}

            )
        end

        x_left,x_right = x_left+delta_left_top,x_right+delta_right_top
    end

    x_left,x_right = left_point[1],right_point[1]
    for y=math.ceil(left_point[2]),math.ceil(p3[2]) do

        for x=math.ceil(x_left-0.5),math.ceil(x_right-0.5) do
            local bary_a,bary_b,bary_c = bary(math.ceil(x),math.ceil(y),left_point,right_point,p3)
            local depth = left_point[3]*bary_a+right_point[3]*bary_b+p3[3]*bary_c

            local tpos =  interpolate_uv(bary_a,bary_b,bary_c,left_point,right_point,p3)

            local z = 1 / depth

            frag(x,y,depth,
                MAX(0,MIN(CEIL(tpos[1]*z*(tex.w+res)),tex.w)-1),
                MAX(0,MIN(CEIL(tpos[2]*z*(tex.h+res)),tex.h)-1),{r=1,g=1,b=1}

            )
        end

        x_left,x_right = x_left+delta_left_bottom,x_right+delta_right_bottom
    end
end

local function raster_triangle(p1,p2,p3,tex,frag)
    if p1[2] > p3[2] then p1,p3 = p3,p1 end
    if p1[2] > p2[2] then p1,p2 = p2,p1 end
    if p2[2] > p3[2] then p2,p3 = p3,p2 end

    local split_alpha = (p2[2]-p1[2])/(p3[2]-p1[2])

    local split_point =  {
        lerp(p1[1],p3[1],split_alpha),
        lerp(p1[2],p3[2],split_alpha),
        lerp(p1[3],p3[3],split_alpha),
        lerp(p1[4],p3[4],split_alpha),
        lerp(p1[5],p3[5],split_alpha),
        lerp(p1[6],p3[6],split_alpha)
    }

    local left_point,right_point = p2,split_point
    if left_point[1] > right_point[1] then
        left_point,right_point = right_point,left_point
    end

    local row_top_min = math.max(0,math.floor(p1[2]+0.5))
    local mid_y_floored = math.floor(p2[2]+0.5)
    local row_top_max = mid_y_floored-1
    local row_bottom_min = math.max(mid_y_floored,0)
    local row_bottom_max = math.ceil(p3[2]-0.5)

    local top_delta_y = left_point[2]-p1[2]
    local top_left_gradient = (left_point[1]-p1[1])/top_delta_y
    local top_right_gradient = (right_point[1]-p1[1])/top_delta_y

    local top_projection = row_top_min + 0.5 - p1[2]
    local top_left_x = p1[1] + top_left_gradient * top_projection
    local top_right_x = p1[1] + top_right_gradient * top_projection

    for y=row_top_min,row_top_max do
        local collumn_min = math.ceil(top_left_x-0.5)
        local collumn_max = math.ceil(top_right_x-0.5)-1

        for x=collumn_min,collumn_max do
            local bary_a,bary_b,bary_c = bary(x,y,left_point,right_point,p3)
            local depth = left_point[3]*bary_a+right_point[3]*bary_b+p3[3]*bary_c

            local tpos =  interpolate_uv(bary_a,bary_b,bary_c,left_point,right_point,p3)

            local z = 1 / depth

            frag(x,y,depth,
                MAX(0,MIN(CEIL(tpos[1]*z*(tex.w+res)),tex.w)-1),
                MAX(0,MIN(CEIL(tpos[2]*z*(tex.h+res)),tex.h)-1),{r=1,g=1,b=1}

            )
        end

        top_left_x = top_left_x + top_left_gradient
        top_right_x = top_right_x + top_right_gradient
    end

    local bottom_delta_y = p3[2] - left_point[2]
    local bottom_left_gradient = (p3[1]-left_point[1]) / bottom_delta_y
    local bottom_right_gradient = (p3[1]-right_point[1]) / bottom_delta_y

    local bottom_projection = row_bottom_min + 0.5 - left_point[2]
    local bottom_left_x = left_point[1] + bottom_left_gradient * bottom_projection
    local bottom_right_x = right_point[1] + bottom_right_gradient * bottom_projection

    for y=row_bottom_min,row_bottom_max do
        local collumn_min = math.ceil(bottom_left_x - 0.5)
        local collumn_max = math.ceil(bottom_right_x - 0.5) - 1

        for x=collumn_min,collumn_max do
            local bary_a,bary_b,bary_c = bary(x,y,left_point,right_point,p3)
            local depth = left_point[3]*bary_a+right_point[3]*bary_b+p3[3]*bary_c

            local tpos =  interpolate_uv(bary_a,bary_b,bary_c,left_point,right_point,p3)

            local z = 1 / depth

            frag(x,y,depth,
                MAX(0,MIN(CEIL(tpos[1]*z*(tex.w+res)),tex.w)-1),
                MAX(0,MIN(CEIL(tpos[2]*z*(tex.h+res)),tex.h)-1),{r=1,g=1,b=1}

            )
        end

        bottom_left_x = bottom_left_x + bottom_left_gradient
        bottom_right_x = bottom_right_x + bottom_right_gradient
    end
end

local function round_x(x)
    return math.ceil(x-0.5)
end
local function round_y(y)
    return math.floor(y+0.5)
end

local function round_xy(x,y)
    return round_x(x),round_y(y)
end

local function round_vertex(vertex)
    vertex[1],vertex[2] = round_xy(vertex[1],vertex[2])
end

local function bary(x,y,p1,p2,p3)
    --[[round_vertex(p1)
    round_vertex(p2)
    round_vertex(p3)]]

    local p23y_delta,p13x_delta,p32x_delta = p2[2]-p3[2],p1[1]-p3[1],p3[1]-p2[1]

    local xp3_delta,yp3_delta = x-p3[1],y-p3[2]

    local div = (p23y_delta*p13x_delta + p32x_delta*(p1[2]-p3[2]))

    local dot_a = (p23y_delta*xp3_delta    + p32x_delta*yp3_delta) / div
    local dot_b = ((p3[2]-p1[2])*xp3_delta + p13x_delta*yp3_delta) / div

    return dot_a,dot_b,1-dot_a-dot_b
end

local function raster_triangle(p1,p2,p3,tex,frag)
    if p1[2] > p3[2] then p1,p3 = p3,p1 end
    if p1[2] > p2[2] then p1,p2 = p2,p1 end
    if p2[2] > p3[2] then p2,p3 = p3,p2 end

    local split_alpha = (p2[2]-p1[2])/(p3[2]-p1[2])

    local split_point =  {
        lerp(p1[1],p3[1],split_alpha),
        lerp(p1[2],p3[2],split_alpha),
        lerp(p1[3],p3[3],split_alpha),
        lerp(p1[4],p3[4],split_alpha),
        lerp(p1[5],p3[5],split_alpha),
        lerp(p1[6],p3[6],split_alpha)
    }

    local left_point,right_point = p2,split_point
    if left_point[1] > right_point[1] then
        left_point,right_point = right_point,left_point
    end

    local delta_left_top  = 1/slope(p1[1],p1[2],left_point[1], left_point[2])
    local delta_right_top = 1/slope(p1[1],p1[2],right_point[1],right_point[2])

    local delta_left_bottom  = 1/slope(p3[1],p3[2],left_point[1], left_point[2])
    local delta_right_bottom = 1/slope(p3[1],p3[2],right_point[1],right_point[2])

    -- flat bottom
    local top_projection = math.floor(p1[2]+0.5) + 0.5 - p1[2]
    local bottom_projection = math.floor(p2[2]+0.5) + 0.5 - left_point[2]

    local x_left,x_right = p1[1] + delta_left_top * top_projection,p1[1] + delta_right_top * top_projection
    if delta_left_top then
    for y=math.floor(p1[2]+0.5),math.floor(p2[2]+0.5)-1 do

        for x=math.ceil(x_left-0.5),math.ceil(x_right-0.5)-1 do

            local bary_a,bary_b,bary_c = bary(x,y,p1,left_point,right_point)
            local depth = p1[3]*bary_a+left_point[3]*bary_b+right_point[3]*bary_c

            local z = 1 / depth

            local tpos = interpolate_uv(bary_a,bary_b,bary_c,p1,left_point,right_point)

            frag(x,y,depth,
                MAX(0,MIN(CEIL(tpos[1]*z*(tex.w+res)),tex.w)-1),
                MAX(0,MIN(CEIL(tpos[2]*z*(tex.h+res)),tex.h)-1),{r=depth,g=depth,b=depth}

            )
        end

        x_left,x_right = x_left+delta_left_top,x_right+delta_right_top
    end
    end


    -- flat top
    if delta_left_bottom == delta_left_bottom then
        x_left,x_right = left_point[1] + delta_left_bottom * bottom_projection,right_point[1] + delta_right_bottom * bottom_projection
        for y=math.floor(p2[2]+0.5),math.ceil(p3[2]-0.5) do

            for x=math.ceil(x_left-0.5),math.ceil(x_right-0.5)-1 do

                local bary_a,bary_b,bary_c = bary(x,y,left_point,right_point,p3)
                local depth = left_point[3]*bary_a+right_point[3]*bary_b+p3[3]*bary_c

                local tpos =  interpolate_uv(bary_a,bary_b,bary_c,left_point,right_point,p3)

                local z = 1 / depth

                frag(x,y,depth,
                    MAX(0,MIN(CEIL(tpos[1]*z*(tex.w+res)),tex.w)-1),
                    MAX(0,MIN(CEIL(tpos[2]*z*(tex.h+res)),tex.h)-1),{r=depth,g=depth,b=depth}

                )
            end

            x_left,x_right = x_left+delta_left_bottom,x_right+delta_right_bottom
        end
    end
end

local tex
function love.load()
    tex = love.image.newImageData("dice_skin.png")
end

local function printt(t)
    local str = ""
    for k,v in pairs(t) do str = str .. tostring(v) .. "\n" end
    error(str)
end

local pitch,yaw = 0,0
local camera_pos = new_vector(0,0,-2)

local scale = mat.new(4,4,
    3,0,0,0,
    0,3,0,0,
    0,0,3,0,
    0,0,0,3
)

local die = false
local selected_pixels = {}

local function render_pixel(screen,x,y,z,tx,ty,debug)
    local r,g,b,a = debug.r,debug.g,debug.b,1
    local _r,_g,_b,_a = tex:getPixel(tx,ty)
    --local r,g,b = r*_r,g*_g,b*_b
        --local r,g,b,a = debug.u,0,debug.v,1

    if a > 0.5 then

        if not screen[y] then screen[y] = {} end
        if screen[y] and screen[y][x] and screen[y][x].w and screen[y][x].w < z then
            if selected_pixels[x] and selected_pixels[x][y] then
                screen[y][x] = {w=0,1,0,0,1}
            else
                screen[y][x] = {w=z,r,g,b,1}
            end
        elseif not screen[y][x] then
            screen[y][x] = {w=z,r,g,b,1}
        end

    end
end

local function interpolate_vertex(v1,v2,alpha)
    return {
        (1 - alpha) * v1[1] + alpha * v2[1],
        (1 - alpha) * v1[2] + alpha * v2[2],
        (1 - alpha) * v1[3] + alpha * v2[3],
        (1 - alpha) * v1[4] + alpha * v2[4],
        (1 - alpha) * v1[5] + alpha * v2[5],
        (1 - alpha) * v1[6] + alpha * v2[6],
    }
end

local function clip_1_vertex(tris,v1,v2,v3)
    local alpha1 = (-v1[3]) / (v2[3]-v1[3])
    local alpha2 = (-v1[3]) / (v3[3]-v1[3])

    local v10 = interpolate_vertex(v1,v2,alpha1)
    local v01 = interpolate_vertex(v1,v3,alpha2)

    tris[#tris+1] = {v10,v2,v3,split=true}
    tris[#tris+1] = {v3,v01,v10,split=true}
end

local function clip_2_vertices(tris,v1,v2,v3)
    local alpha1 = (-v1[3]) / (v3[3]-v1[3])
    local alpha2 = (-v2[3]) / (v3[3]-v2[3])

    local v10 = interpolate_vertex(v1,v3,alpha1)
    local v01 = interpolate_vertex(v2,v3,alpha2)

    tris[#tris+1] = {v3,v01,v10,split=true}
end

local function handle_triangle(tri_list,a,b,c)
    local v1x,v1y,v1z,v1w = a[1],a[2],a[3],a[4]
    local v2x,v2y,v2z,v2w = b[1],b[2],b[3],b[4]
    local v3x,v3y,v3z,v3w = c[1],c[2],c[3],c[4]

    if v1x >  v1w and v2x >  v2w and v3x >  v3w then return end

    if v1x < -v1w and v2x < -v2w and v3x < -v3w then return end

    if v1y >  v1w and v2y >  v2w and v3y >  v3w then return end

    if v1y < -v1w and v2y < -v2w and v3y < -v3w then return end

    if v1z >  v1w and v2z > -v2w and v3z >  v3w then return end

    if v1z < 0    and v2z < 0    and v3z < 0    then return end

    if v1z < 0 then
        if v2z < 0 then
            clip_2_vertices(tri_list,a,b,c)
        elseif v3z < 0 then

            clip_2_vertices(tri_list,a,c,b)
        else
            clip_1_vertex(tri_list,a,b,c)
        end
    elseif v2z < 0 then
        if v3z < 0 then
            clip_2_vertices(tri_list,b,c,a)
        else
            clip_1_vertex(tri_list,b,a,c)
        end
    elseif v3z < 0 then
        clip_1_vertex(tri_list,c,a,b)
    else
        tri_list[#tri_list+1] = {a,b,c}
    end
end

function love.draw()
    local per = makePerspective(w,h,const/100,10,70)
    local screen = {}

    local fx,fy,fz = camera_pos.x,camera_pos.y,camera_pos.z
    local look_vector = add_vector(camera_pos,get_look_vector(math.rad(yaw),math.rad(pitch)))
    local ax,ay,az = look_vector.x,look_vector.y,look_vector.z

    local lookat = makeLookat(fx,fy,fz,ax,ay,az)

    for k,v in pairs(objects) do
        local vertices = {}

        local vert = v.verts
        local tri  = v.indices

        local n = 0
        for i=1,#vert,5 do
            n = n + 1

            local scaled_vertice = mat.vector(
                vert[i],vert[i+1],vert[i+2],1
            )*scale

            local rotated_vertice = scaled_vertice*create_rotation_eulers{
                x=0,
                y=0,
                z=0
            }

            local model = rotated_vertice*lookat

            local v = model*per

            v[5] = vert[i+3]
            v[6] = vert[i+4]

            vertices[n] = v
        end

        local triangles = {}

        for i=1,#tri,3 do
            n = n + 1

            handle_triangle(triangles,vertices[tri[i]],vertices[tri[i+1]],vertices[tri[i+2]])
        end

        for i=1,#triangles do
            local t = triangles[i]
            local a,b,c = t[1],t[2],t[3]
            --if cull(a,b,c) > 0 then
                raster_triangle(
                    transform_screen_space(a,w,h),
                    transform_screen_space(b,w,h),
                    transform_screen_space(c,w,h),
                    {tex=tex,w=tex:getWidth(),h=tex:getHeight()},
                    function(x,y,z,tx,ty,debug)
                        render_pixel(screen,x,y,z,tx,ty,debug)
                    end
                )
            --end
        end
    end
    local sc = {}
    local n = 0
    for k,v in pairs(screen) do
        for _k,_v in pairs(v) do
            n = n + 1
            sc[n] = {_k*res,k*res,_v[1],_v[2],_v[3],_v[4]}
        end
    end
    love.graphics.setPointSize(res)
    love.graphics.points(sc)

    if die then error("finished frame") end
end

local pitch_lim = {-89,89}
local speed = 0.05

function love.update(dt)
    if love.keyboard.isDown("up") then
        local new_pitch = pitch - 100*dt
        if new_pitch > pitch_lim[1] then pitch = new_pitch end
    end
    if love.keyboard.isDown("down") then
        local new_pitch = pitch + 100*dt
        if new_pitch < pitch_lim[2] then pitch = new_pitch end
    end
    if love.keyboard.isDown("left") then
        yaw = yaw - 100*dt
    end
    if love.keyboard.isDown("right") then
        yaw = yaw + 100*dt
    end
    if love.keyboard.isDown("a") then
        local move_vector = get_move_vector(math.rad(yaw-90))
        camera_pos = add_vector(camera_pos,scale_vector(move_vector,speed))
    end
    if love.keyboard.isDown("d") then
        local move_vector = get_move_vector(math.rad(yaw+90))
        camera_pos = add_vector(camera_pos,scale_vector(move_vector,speed))
    end
    if love.keyboard.isDown("w") then
        local move_vector = get_move_vector(math.rad(yaw))
        camera_pos = add_vector(camera_pos,scale_vector(move_vector,speed))
    end
    if love.keyboard.isDown("s") then
        local move_vector = get_move_vector(math.rad(yaw))
        camera_pos = sub_vector(camera_pos,scale_vector(move_vector,speed))
    end
    if love.keyboard.isDown("lctrl") then
        camera_pos = add_vector(camera_pos,scale_vector(new_vector(0,-1,0),speed))
    end
    if love.keyboard.isDown("lshift") then
        camera_pos = add_vector(camera_pos,scale_vector(new_vector(0,1,0),speed))
    end
    if love.keyboard.isDown("space") then
        die = true
    end
end

function love.mousepressed(x,y,b)
    local x,y = math.floor(x/res),math.floor(y/res)
    if not selected_pixels[x] then selected_pixels[x] = {} end
    selected_pixels[x][y] = b ~= 2
end

function love.wheelmoved(dx,dy)
    const = const + dy
end

function love.resize()
    local res = 5
    worg,horg = love.graphics.getDimensions()
    w = worg / res
    h = horg / res
end