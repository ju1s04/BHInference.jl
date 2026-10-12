
# Adapted from Krang's "Coordinates with inclination" example:
# https://dominic-chang.com/Krang.jl/stable/examples/coordinate-example

using Krang
using CairoMakie

# ---- Black hole and camera ----
metric = Krang.Kerr(0.99)          # spin a = 0.99
θo     = 45 * π / 180              # observer inclination (45 degrees)
sze    = 200                       # image is sze x sze pixels
rmin   = Krang.horizon(metric)     # inner edge: the event horizon
rmax   = 10.0                      # outer edge of the region we trace
ρmax   = 15.0                      # half-width of the screen (units of M)

camera = Krang.SlowLightIntensityCamera(metric, θo, -ρmax, ρmax, -ρmax, ρmax, sze)

# ---- Where does each ray hit a cone of opening angle θs? ----
function coordinate_point(pix::Krang.AbstractPixel, geometry::Krang.ConeGeometry{T,A}) where {T,A}
    n, rmin, rmax = geometry.attributes
    θs = geometry.opening_angle
    coords = ntuple(j -> zero(T), Val(4))

    isindir = false
    for _ in 1:2
        isindir ⊻= true
        ts, rs, ϕs = emission_coordinates(pix, θs, isindir, n)
        if isnan(rs) || rs ≤ rmin || rs ≥ rmax
            continue
        end
        coords = (ts, rs, θs, ϕs)
    end
    return coords
end

# ---- Trace the equatorial plane (θs = π/2), direct image (n = 0) ----
geometry = Krang.ConeGeometry(π / 2, (0, rmin, rmax))
scene    = coordinate_point.(camera.screen.pixels, Ref(geometry))
radii    = [s[2] for s in scene]

# ---- Plot and save ----
fig = Figure(size = (500, 500))
ax  = Axis(fig[1, 1], aspect = 1,
           title = "Equatorial plane radius, a = 0.99, θo = 45°",
           xticksvisible = false, yticksvisible = false,
           xticklabelsvisible = false, yticklabelsvisible = false)
hm = heatmap!(ax, radii, colormap = :afmhot, colorrange = (0, rmax))
Colorbar(fig[1, 2], hm, label = "r (M)")
save("disc_radius.png", fig)
println("Saved disc_radius.png")