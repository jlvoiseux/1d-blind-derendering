@everywhere using LinearAlgebra
@everywhere using Images, TestImages, Colors
@everywhere using ImageView
@everywhere using Distributed
@everywhere using SharedArrays

@everywhere function normalizeVector(v)
    normalizedVector = v./norm(v);
    return normalizedVector
end

@everywhere function intersectSphere(sphere, ray)
    op = sphere[2:4] - ray[1,:]
    eps = 1e-4
    b = dot(op, ray[2,:])
    det = b*b-dot(op, op)+sphere[1]*sphere[1];
    if (det<0)
        intersection = 0
    else()
        det = sqrt(det)
        if b-det .> eps
            intersection = b-det
        elseif b+det .> eps
            intersection = b+det
        else()
            intersection = 0
        end
    end
    return intersection
end

@everywhere function intersectRay(ray, sphereList)
    n = size(sphereList, 1)
    inf = 1e20
    t = 1e20
    id = 0

    for i = n:-1:1
        d = intersectSphere(sphereList[i,:], ray)
        if(d!=0 && d .< t)
            t = d
            id = i
        end
    end

    out = (t .< inf)
    return out, t, id
end

@everywhere function clamp(x)
    if(x .< 0)
        out = 0
    elseif x .> 1
        out = 1
    else()
        out = x
    end
    return out
end


@everywhere function computeRadiance(sphereList, ray, depth, Xi)
    hit, t, id = intersectRay(ray, sphereList)
    out = [0 0 0]
    obj = sphereList[id, :]
    if(hit .== 0)
        out = [0 0 0];
    elseif depth .> 10
        return [0 0 0];
    else()
        x = ray[1,:]+ray[2,:].*t
        n = normalizeVector(x-obj[2:4])
        if dot(n, ray[2,:]) .< 0
            nl = n
        else()
            nl = n*(-1)
        end
        f = obj[8:10]

        p1 = f[1] .> f[2]
        if(f[1] .> f[3])
            p2 = f[1]
        elseif f[2] .> f[3]
            p2 = f[2]
        else()
            p2 = f[3]
        end

        p = (p1 && p2)
        depth = depth + 1
        if(depth .> 5)
            if(rand() .< p)
                f = f.*(1/p)
            else()
                return obj[5:7]
                #return
            end

        # DIFFUSE
        elseif obj[11] .== 1
            r1 = 2*pi*rand()
            r2 = rand()
            r2s = sqrt(r2)
            w = nl
            if(abs(w[1]) .> 0.1)
                u = normalizeVector([0,1,0])
            else()
                u = normalizeVector(cross([1,0,0],w))
            end
            v = cross(w, u)
            d = normalizeVector(u.*cos(r1).*r2s + v.*sin(r1).*r2s + w.*sqrt(1-r2));
            return obj[5:7]+f.*(computeRadiance(sphereList, [x d]', depth, Xi))
            #return

        # MIRROR
        elseif obj[11] .== 2
            return obj[5:7] + f.*(computeRadiance(sphereList, [x ray[2,:]-n.*2 .*dot(n,ray[2,:])]', depth, Xi))
            #return
        # GLASS
        elseif obj[11] .== 3
            reflRay = [x  ray[2,:]-n.*2 .*dot(n,ray[2,:])]'
            into = dot(n, nl) .> 0
            nc = 1
            nt = 1.5
            if into
                nnt = nc/nt
                into1 = 1
            else()
                nnt = nt/nc
                into1 = -1
            end
            ddn = dot(ray[2,:], nl)
            cos2t = 1-nnt.*nnt.*(1-ddn.*ddn)
            if(cos2t .< 0)
                return obj[5:7] + f.*(computeRadiance(sphereList, reflRay, depth, Xi))
                #return
            else()
                tdir = normalizeVector(ray[2,:].*nnt - n.*into1.*(ddn.*nnt + sqrt(cos2t)))
                if into
                    into2 = -ddn
                else()
                    into2 = dot(tdir, n)
                end
                a = nt-nc
                b = nt + nc
                R0 = a.*a / (b.*b);
                c = 1-into2
                Re = R0 + (1-R0).*c.*c.*c.*c.*c
                Tr = 1-Re
                P = 0.25 + 0.5*Re
                RP = Re./P
                TP = Tr./(1-P)
                if(depth .> 2)
                    if(rand() .< P)
                        return obj[5:7] + f.*(computeRadiance(sphereList, reflRay, depth, Xi).*RP)
                        #return
                    else()
                        return obj[5:7] + f.*(computeRadiance(sphereList, [x tdir]', depth, Xi).*TP)
                        #return
                    end
                else()
                    return obj[5:7] + f.*(computeRadiance(sphereList, reflRay, depth, Xi).*Re + computeRadiance(sphereList, [x tdir]', depth, Xi).*Tr)
                end
            end
        end
    end
end

@everywhere function main()
    w = 426 # Image width
    h = 240 # Image height
    samples = 50 # Number of samples
    fov = 0.5135 # Field of view()

    count = SharedArray{Int64}(1)
    renderedImage = SharedArray{Float64}(h, w, 3)
    camera = [50 52 295.6 ; normalizeVector([0 -0.042612 -1])]
    camX = [w.*0.5135./h, 0, 0]
    camY = normalizeVector(cross(camX, vec(camera[2,:]))).*fov

    sphereList = [1e5   1e5+1     40.8      81.6        0   0   0   0.75    0.25    0.25    1 ;
                  1e5   -1e5+99   40.8      81.6        0   0   0   0.25    0.25    0.75    1 ;
                  1e5   50        40.8      1e5         0   0   0   0.75    0.75    0.75    1 ;
                  1e5   50        40.8      -1e5+170    0   0   0   0       0       0       1 ;
                  1e5   50        1e5       81.6        0   0   0   0.75    0.75    0.75    1 ;
                  1e5   50        -1e5+81.6 81.6        0   0   0   0.75    0.75    0.75    1 ;
                  16.5  27        16.5      47          0   0   0   0.999   0.999   0.999   2 ;
                  16.5  73        16.5      78          0   0   0   0.999   0.999   0.999   3 ;
                  600   50        681.6-.27 81.6        12  12  12  0       0       0       1]
    numSpheres = size(sphereList, 1)
    @sync @distributed for y = 1:h
        Xi = [0 0 y.*y.*y]
        imgLine = zeros(w, 3)
        for x = 1:w
            for sy = 1:2
                for sx = 1:2
                    r = zeros(3)
                    for s = 1:samples
                        r1 = 2*rand()
                        r2 = 2*rand()
                        if(r1 .< 1)
                            dx = sqrt(r1)-1
                        else()
                            dx = 1-sqrt(2-r1)
                        end
                        if(r2 .< 1)
                            dy = sqrt(r2)-1
                        else()
                            dy = 1-sqrt(2-r2)
                        end
                        d = camX.*(((sx + 0.5 + dx)./2 + x)./w - 0.5) + camY.*(((sy + 0.5 + dy)./2 + y)./h - 0.5) + camera[2, :]
                        r = r + computeRadiance(sphereList, [camera[1,:]+d.*140  normalizeVector(d)]', 0, Xi).*(1 ./samples)
                    end
                    imgLine[x, 1] = imgLine[x, 1] + 0.25*clamp(r[1])
                    imgLine[x, 2] = imgLine[x, 2] + 0.25*clamp(r[2])
                    imgLine[x, 3] = imgLine[x, 3] + 0.25*clamp(r[3]);
                end
            end
        end
        count[1] = count[1] + 1
        println(count[1])
        renderedImage[y,:, :] = imgLine
    end
    #renderedImage = Images.imrotate(renderedImage, π)
    ImageView.imshow(renderedImage)
end

main()
