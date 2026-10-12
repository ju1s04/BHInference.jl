using Krang

struct ZAMORedshifts{T} <: Krang.AbstractMaterial
    subimgs::Tuple{Int}
    rmin::T
    rmax::T
end

function (m::ZAMORedshifts)(
    pix::Krang.AbstractPixel{T},
    intersection::Krang.Intersection,
) where {T}
    (; rs, θs, νr, νθ) = intersection
    α, β = Krang.screen_coordinate(pix)

    ηtemp = η(metric, α, β, θo)
    λtemp = λ(metric, α, θo)
    curr_p_bl_d = p_bl_d(metric, rs, θs, ηtemp, λtemp, νr, νθ)

    curr_p_bl_u = metric_uu(metric, rs, θs) * curr_p_bl_d
    p_zamo_u = jac_zamo_u_bl_d(metric, rs, θs) * curr_p_bl_u
    redshift = inv(p_zamo_u[1])

    return max(redshift, eps(T))
end

metric = Krang.Kerr(0.99);
θo = 17 * π / 180;
ρmax = 15.0;
rmin = Krang.horizon(metric);
rmax = 10.0;
n = 1;

camera = Krang.IntensityCamera(metric, θo, -ρmax, ρmax, -ρmax, ρmax, 400);
mesh = Krang.Mesh(Krang.ConeGeometry((75 * π / 180)), ZAMORedshifts((n,), rmin, rmax))
scene = Krang.Scene((mesh,))

redshifts = render(camera, scene)

import CairoMakie as CMk

theme = CMk.Theme(
    Axis = (
        xticksvisible = false,
        xticklabelsvisible = false,
        yticksvisible = false,
        yticklabelsvisible = false,
    ),
)

CMk.set_theme!(CMk.merge!(theme, CMk.theme_latexfonts()))

fig = CMk.Figure(resolution = (700, 700));
ax = CMk.Axis(fig[1, 1], title = "Redshifts", titlesize = 20, aspect = 1)
hm = CMk.heatmap!(ax, redshifts, colormap = :afmhot)
CMk.Colorbar(fig[1, 2], hm, label = "Redshifts", labelsize = 20)
CMk.save("redshifts.png", fig)