-- #section depth_test depth_store interpolate_uvs enable_fs
local function rasterize_triangle(
    fb_front,
    fb_depth,
    fb_width,
    fb_height_m1,
    p0x,
    p0y,
    p0w,
    p0u,
    p0v,
    p1x,
    p1y,
    p1w,
    p1u,
    p1v,
    p2x,
    p2y,
    p2w,
    p2u,
    p2v,
    fixed_colour,
    fragment_shader,
    pipeline_uniforms)
    local math_ceil = math.ceil
    local math_floor = math.floor
    local fb_width_m1 = fb_width - 1

    -- see: https://github.com/exerro/v3d/blob/main/raster_visuals/src/main/kotlin/me/exerro/raster_visuals/rasterize.kt
    -- there's an explanation of the algorithm there
    -- this code has been heavily microoptimised so won't perfectly resemble that
    -- this has also had depth testing and UV interpolation added in, so good
    -- luck understanding anything here :/

    -- #if depth_test depth_store interpolate_uvs
    -- #if interpolate_uvs
    if p0y > p1y then
        p0x, p0y, p0w, p0u, p0v, p1x, p1y, p1w, p1u, p1v = p1x, p1y, p1w, p1u, p1v, p0x, p0y, p0w, p0u, p0v
    end
    if p1y > p2y then
        p1x, p1y, p1w, p1u, p1v, p2x, p2y, p2w, p2u, p2v = p2x, p2y, p2w, p2u, p2v, p1x, p1y, p1w, p1u, p1v
    end
    if p0y > p1y then
        p0x, p0y, p0w, p0u, p0v, p1x, p1y, p1w, p1u, p1v = p1x, p1y, p1w, p1u, p1v, p0x, p0y, p0w, p0u, p0v
    end
    -- #else
    if p0y > p1y then
        p0x, p0y, p0w, p1x, p1y, p1w = p1x, p1y, p1w, p0x, p0y, p0w
    end
    if p1y > p2y then
        p1x, p1y, p1w, p2x, p2y, p2w = p2x, p2y, p2w, p1x, p1y, p1w
    end
    if p0y > p1y then
        p0x, p0y, p0w, p1x, p1y, p1w = p1x, p1y, p1w, p0x, p0y, p0w
    end
    -- #end
    -- #else
    if p0y > p1y then
        p0x, p0y, p1x, p1y = p1x, p1y, p0x, p0y
    end
    if p1y > p2y then
        p1x, p1y, p2x, p2y = p2x, p2y, p1x, p1y
    end
    if p0y > p1y then
        p0x, p0y, p1x, p1y = p1x, p1y, p0x, p0y
    end
    -- #end
    if p0y == p2y then
        return
    end -- skip early if we have a perfectly flat triangle

    local f = (p1y - p0y) / (p2y - p0y)
    local pMx = p0x * (1 - f) + p2x * f
    -- #if depth_test depth_store interpolate_uvs
    local pMw = p0w * (1 - f) + p2w * f
    -- #end
    -- #if interpolate_uvs
    local pMu = (p0u * p0w * (1 - f) + p2u * p2w * f) / pMw
    local pMv = (p0v * p0w * (1 - f) + p2v * p2w * f) / pMw
    -- #end

    if pMx > p1x then
        pMx, p1x = p1x, pMx
        -- #if depth_test depth_store interpolate_uvs
        pMw, p1w = p1w, pMw
        -- #end
        -- #if interpolate_uvs
        pMu, p1u = p1u, pMu
        pMv, p1v = p1v, pMv
    -- #end
    end

    local rowTopMin = math_floor(p0y + 0.5)
    local rowBottomMin = math_floor(p1y + 0.5)
    local rowTopMax = rowBottomMin - 1
    local rowBottomMax = math_ceil(p2y - 0.5)

    if rowTopMin < 0 then
        rowTopMin = 0
    end
    if rowBottomMin < 0 then
        rowBottomMin = 0
    end
    if rowTopMax > fb_height_m1 then
        rowTopMax = fb_height_m1
    end
    if rowBottomMax > fb_height_m1 then
        rowBottomMax = fb_height_m1
    end

    local function rasterise_flat_triangle(
        triLeftGradientX,
        triRightGradientX,
        triLeftGradientW,
        triRightGradientW,
        triLeftGradientUW,
        triRightGradientUW,
        triLeftGradientVW,
        triRightGradientVW,
        triLeftX,
        triRightX,
        triLeftW,
        triRightW,
        triLeftUW,
        triRightUW,
        triLeftVW,
        triRightVW,
        it_a,
        it_b)
        for baseIndex = it_a, it_b, fb_width do
            local columnMinX = math_ceil(triLeftX)
            local columnMaxX = math_ceil(triRightX)
            -- #if depth_test depth_store interpolate_uvs
            local rowTotalDeltaX = triRightX - triLeftX + 1 -- 'cause of awkward optimisations above
            local rowDeltaW = (triRightW - triLeftW) / rowTotalDeltaX
            local rowLeftW = triLeftW + (columnMinX - triLeftX) * rowDeltaW
            -- #end
            -- #if interpolate_uvs
            local rowDeltaU = (triRightUW - triLeftUW) / rowTotalDeltaX
            local rowLeftU = triLeftUW + (columnMinX - triLeftX) * rowDeltaU
            local rowDeltaV = (triRightVW - triLeftVW) / rowTotalDeltaX
            local rowLeftV = triLeftVW + (columnMinX - triLeftX) * rowDeltaV
            -- #end

            if columnMinX < 0 then
                columnMinX = 0
            end
            if columnMaxX > fb_width_m1 then
                columnMaxX = fb_width_m1
            end

            for x = columnMinX, columnMaxX do
                local index = baseIndex + x

                local u, v = 0, 0
                -- #if interpolate_uvs
                u = rowLeftU / rowLeftW
                v = rowLeftV / rowLeftW
                -- #end

                -- #if depth_test
                if rowLeftW > fb_depth[index] then
                    -- #if enable_fs
                    local fs_colour = fragment_shader(pipeline_uniforms, u, v)
                    if fs_colour ~= 0 then
                        fb_front[index] = fs_colour
                        -- #if depth_store
                        fb_depth[index] = rowLeftW
                    -- #end
                    end
                    -- #else
                    fb_front[index] = fixed_colour
                    -- #if depth_store
                    fb_depth[index] = rowLeftW
                -- #end
                -- #end
                end
                -- #else
                -- #if enable_fs
                local fs_colour = fragment_shader(pipeline_uniforms, 0, 0)
                if fs_colour ~= 0 then
                    fb_front[index] = fs_colour
                    -- #if depth_store
                    fb_depth[index] = rowLeftW
                -- #end
                end
                -- #else
                fb_front[index] = fixed_colour
                -- #if depth_store
                fb_depth[index] = rowLeftW
                -- #end
                -- #end
                -- #end

                -- #if depth_test depth_store interpolate_uvs
                rowLeftW = rowLeftW + rowDeltaW
                -- #end
                -- #if interpolate_uvs
                rowLeftU = rowLeftU + rowDeltaU
                rowLeftV = rowLeftV + rowDeltaV
                -- #end
            end

            triLeftX = triLeftX + triLeftGradientX
            triRightX = triRightX + triRightGradientX
            -- #if depth_test depth_store interpolate_uvs
            triLeftW = triLeftW + triLeftGradientW
            triRightW = triRightW + triRightGradientW
            -- #end
            -- #if interpolate_uvs
            triLeftUW = triLeftUW + triLeftGradientUW
            triRightUW = triRightUW + triRightGradientUW
            triLeftVW = triLeftVW + triLeftGradientVW
            triRightVW = triRightVW + triRightGradientVW
            -- #end
        end
    end

    if rowTopMin <= rowTopMax then
        local triDeltaY = p1y - p0y
        local triLeftGradientX = (pMx - p0x) / triDeltaY
        local triRightGradientX = (p1x - p0x) / triDeltaY
        local triLeftGradientW, triRightGradientW
        -- #if depth_test depth_store interpolate_uvs
        triLeftGradientW = (pMw - p0w) / triDeltaY
        triRightGradientW = (p1w - p0w) / triDeltaY
        -- #end
        local triLeftGradientUW, triRightGradientUW
        local triLeftGradientVW, triRightGradientVW
        -- #if interpolate_uvs
        triLeftGradientUW = (pMu * pMw - p0u * p0w) / triDeltaY
        triRightGradientUW = (p1u * p1w - p0u * p0w) / triDeltaY
        triLeftGradientVW = (pMv * pMw - p0v * p0w) / triDeltaY
        triRightGradientVW = (p1v * p1w - p0v * p0w) / triDeltaY
        -- #end

        local triProjection = rowTopMin + 0.5 - p0y
        local triLeftX = p0x + triLeftGradientX * triProjection - 0.5
        local triRightX = p0x + triRightGradientX * triProjection - 1.5
        local triLeftW, triRightW
        -- #if depth_test depth_store interpolate_uvs
        triLeftW = p0w + triLeftGradientW * triProjection
        triRightW = p0w + triRightGradientW * triProjection
        -- #end
        local triLeftUW, triRightUW
        local triLeftVW, triRightVW
        -- #if interpolate_uvs
        triLeftUW = p0u * p0w + triLeftGradientUW * triProjection
        triRightUW = p0u * p0w + triRightGradientUW * triProjection
        triLeftVW = p0v * p0w + triLeftGradientVW * triProjection
        triRightVW = p0v * p0w + triRightGradientVW * triProjection
        -- #end

        local it_a, it_b = rowTopMin * fb_width + 1, rowTopMax * fb_width + 1

        rasterise_flat_triangle(
            triLeftGradientX,
            triRightGradientX,
            triLeftGradientW,
            triRightGradientW,
            triLeftGradientUW,
            triRightGradientUW,
            triLeftGradientVW,
            triRightGradientVW,
            triLeftX,
            triRightX,
            triLeftW,
            triRightW,
            triLeftUW,
            triRightUW,
            triLeftVW,
            triRightVW,
            it_a,
            it_b
        )
    end

    if rowBottomMin <= rowBottomMax then
        local triDeltaY = p2y - p1y
        local triLeftGradientX = (p2x - pMx) / triDeltaY
        local triRightGradientX = (p2x - p1x) / triDeltaY
        local triLeftGradientW, triRightGradientW
        -- #if depth_test depth_store interpolate_uvs
        triLeftGradientW = (p2w - pMw) / triDeltaY
        triRightGradientW = (p2w - p1w) / triDeltaY
        -- #end
        local triLeftGradientUW, triRightGradientUW
        local triLeftGradientVW, triRightGradientVW
        -- #if interpolate_uvs
        triLeftGradientUW = (p2u * p2w - pMu * pMw) / triDeltaY
        triRightGradientUW = (p2u * p2w - p1u * p1w) / triDeltaY
        triLeftGradientVW = (p2v * p2w - pMv * pMw) / triDeltaY
        triRightGradientVW = (p2v * p2w - p1v * p1w) / triDeltaY
        -- #end

        local triProjection = rowBottomMin + 0.5 - p1y
        local triLeftX = pMx + triLeftGradientX * triProjection - 0.5
        local triRightX = p1x + triRightGradientX * triProjection - 1.5
        local triLeftW, triRightW
        -- #if depth_test depth_store interpolate_uvs
        triLeftW = pMw + triLeftGradientW * triProjection
        triRightW = p1w + triRightGradientW * triProjection
        -- #end
        local triLeftUW, triRightUW
        local triLeftVW, triRightVW
        -- #if interpolate_uvs
        triLeftUW = pMu * pMw + triLeftGradientUW * triProjection
        triRightUW = p1u * p1w + triRightGradientUW * triProjection
        triLeftVW = pMv * pMw + triLeftGradientVW * triProjection
        triRightVW = p1v * p1w + triRightGradientVW * triProjection
        -- #end

        local it_a, it_b = rowBottomMin * fb_width + 1, rowBottomMax * fb_width + 1

        rasterise_flat_triangle(
            triLeftGradientX,
            triRightGradientX,
            triLeftGradientW,
            triRightGradientW,
            triLeftGradientUW,
            triRightGradientUW,
            triLeftGradientVW,
            triRightGradientVW,
            triLeftX,
            triRightX,
            triLeftW,
            triRightW,
            triLeftUW,
            triRightUW,
            triLeftVW,
            triRightVW,
            it_a,
            it_b
        )
    end
end
-- #endsection

rasterize_triangle()