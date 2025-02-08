# Blind Derendering - How to turn a wall into a mirror
Master Thesis document : https://dial.uclouvain.be/memoire/ucl/fr/object/thesis%3A25235

The [rendering equation](https://en.wikipedia.org/wiki/Rendering_equation) sums all the inbound light for a given surface,
applies a function to those light contributions and returns the outbound light for that surface. In essence, [ray-tracing](https://en.wikipedia.org/wiki/Ray_tracing_(graphics))
is a direct implementation of that equation. The function applied depends on the material (for example, a mirror will not
affect the inbound light, such that the outbound light will be similar to the inbound light). The goal of this master thesis
was the **explore the possibilities of recovering inbound light data, given outbound light data**.

This repository is separated in two sections:
- The MATLAB and Julia code used to perform blind deconvolution
- Code in several languages used for rendering purposes, including a pure C++ ray-tracer
