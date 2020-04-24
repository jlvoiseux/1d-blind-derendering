function indices = FitImpulses(obs, empty_source, mirror_brdf, num_lin, margin, interf)
    g_mirror = SimpleRendering(obs, empty_source, mirror_brdf, num_lin, margin);
    if interf
        g_mirror = xcorr(g_mirror);
        g_mirror(g_mirror>0.5) = impulses;
    else
        g_mirror(g_mirror>0.5) = impulses;
    end
    indices = g_mirror;    
end