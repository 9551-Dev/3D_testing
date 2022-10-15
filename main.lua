local s = 0.5

local mat = require("mat")

local function u(u)
    return (u+1)/3
end
local function v(v)
    return v/4
end

--[[local objects = {
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
}]]

local objects = {
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
}

local res = 10

local worg,horg = love.graphics.getDimensions()
love.window.setMode(worg,horg,{resizable=true})
local w = worg / res
local h = horg / res

local const = 0

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

local function makePerspective(w,h,n,f,fov)
    local aspectRatio = 1/w * h
    fov = math.rad(fov)

    return mat.new(4,4,
        aspectRatio/math.tan(fov*0.5),0,0,0,
        0,1/(math.tan(fov*0.5)),0,0,
        0,0,-f/(f-n),-1,
        0,0,-f*n/(f-n),1
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

local function interpolate_uv(bary,uv1,uv2,uv3)
    return {
        uv1[5] * bary[1] + uv2[5] * bary[2] + uv3[5] * bary[3],
        uv1[6] * bary[1] + uv2[6] * bary[2] + uv3[6] * bary[3]
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
                MAX(0,MIN(CEIL(tpos[1]*z*tex.w),tex.w-1)),
                MAX(0,MIN(CEIL(tpos[2]*z*tex.h),tex.h-1)),{t3=t3,u=tpos[1],v=tpos[2]}
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
                MAX(0,MIN(CEIL(tpos[1]*z*tex.w),tex.w-1)),
                MAX(0,MIN(CEIL(tpos[2]*z*tex.h),tex.h-1)),{t3=t3,u=tpos[1],v=tpos[2]}
            )
        end
    end
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

local tex
function love.load()
    tex = love.image.newImageData("jesus.png")
    per = makePerspective(w,h,-0.5,90,90)
end

local function printt(t)
    local str = ""
    for k,v in pairs(t) do str = str .. tostring(v) .. "\n" end
    error(str)
end

local rx,ry,rz,px,py,pz = 0,0,0,0,0,0

local scale = mat.new(4,4,
    1,0,0,0,
    0,1,0,0,
    0,0,1,0,
    0,0,0,1
)

local die = false
local selected_pixels = {}

local function render_pixel(screen,x,y,z,tx,ty,debug)
    local r,g,b,a = tex:getPixel(tx,ty)
    --local r,g,b,a = debug.u,0,debug.v,1

    if a > 0.5 then

        if not screen[y] then screen[y] = {} end
        if screen[y] and screen[y][x] and screen[y][x].w and screen[y][x].w > z then
            if selected_pixels[x] and selected_pixels[x][y] then
                screen[y][x] = {w=0,1,0,0,1}
                print(debug.t3)
            else
                screen[y][x] = {w=z,r,g,b,1}
            end
        elseif not screen[y][x] then
            if selected_pixels[x*res] and selected_pixels[x][y] then
                screen[y][x] = {w=0,1,0,0,1}
            else
                screen[y][x] = {w=z,r,g,b,1}
            end
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
    tris[#tris+1] = {v10,v01,v3,split=true}
end

local function clip_2_vertices(tris,v1,v2,v3)
    local alpha1 = (-v1[3]) / (v3[3]-v1[3])
    local alpha2 = (-v2[3]) / (v3[3]-v2[3])

    local v10 = interpolate_vertex(v1,v3,alpha1)
    local v01 = interpolate_vertex(v2,v3,alpha2)

    tris[#tris+1] = {v10,v01,v3,split=true}
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
    local per = makePerspective(w,h,-0.5,90,50)
    local screen = {}

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
                x=rx,
                y=ry,
                z=rz
            }

            local model = rotated_vertice*mat.new(4,4,
                1,0,0,0,
                0,1,0,0,
                0,0,1,0,
                px,py,pz,1
            )

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
            raster_triangle(
                transform_screen_space(a,w,h),
                transform_screen_space(b,w,h),
                transform_screen_space(c,w,h),
                {tex=tex,w=tex:getWidth(),h=tex:getHeight()},
                function(x,y,z,tx,ty,debug)
                    render_pixel(screen,x,y,z,tx,ty,debug)
                end
            )
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

function love.update(dt)
    if love.keyboard.isDown("up") then
        rx = rx - 100*dt
    end
    if love.keyboard.isDown("down") then
        rx = rx + 100*dt
    end
    if love.keyboard.isDown("left") then
        ry = ry - 100*dt
    end
    if love.keyboard.isDown("right") then
        ry = ry + 100*dt
    end
    if love.keyboard.isDown("a") then
        px = px - 2*dt
    end
    if love.keyboard.isDown("d") then
        px = px + 2*dt
    end
    if love.keyboard.isDown("w") then
        pz = pz - 2*dt
    end
    if love.keyboard.isDown("s") then
        pz = pz + 2*dt
    end
    if love.keyboard.isDown("lctrl") then
        py = py - 2*dt
    end
    if love.keyboard.isDown("lshift") then
        py = py + 2*dt
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
end

function love.resize()
    local res = 10
    worg,horg = love.graphics.getDimensions()
    w = worg / res
    h = horg / res

    per = makePerspective(w,h,-0.5,90,50)
end