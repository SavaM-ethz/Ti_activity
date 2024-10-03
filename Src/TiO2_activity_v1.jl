using Plots, MAGEMin_C, Pkg
Pkg.activate(".")

function main()
#    data    = Initialize_MAGEMin("ig", verbose=false);
#    P,T     = 10.0, 1100.0
#   X= [48.43; 15.19; 11.57; 10.13; 6.65; 1.64; 0.59; 1.87; 0.68; 0.0; 3.0];
#    Xoxides = ["SiO2"; "Al2O3"; "CaO"; "MgO"; "FeO"; "Fe2O3"; "K2O"; "Na2O"; "TiO2"; "Cr2O3"; "H2O"];
#    sys_in  = "wt"
#   out     = single_point_minimization(P, T, data, X=X, Xoxides=Xoxides, sys_in=sys_in)
#   Finalize_MAGEMin(data)


    db      = "ig"  # database: ig, igneous (Holland et al., 2018); mp, metapelite (White et al 2014b)
    data    = Initialize_MAGEMin(db, buffer="qfm", verbose=false);
    nT      = 50
    P_const = 2.0
    B_const = -2.0
    P       = [P_const for _ in 1:nT]
    B       = [B_const for _ in 1:nT]
    T       = collect(LinRange(600.0, 1000.0, nT))
    X       = [69.3; 15.6; 1.65; 1.18; 2.50; 1.5; 4.61; 3.49; 0.38; 0.0; 4.0];
    Xoxides = ["SiO2"; "Al2O3"; "CaO"; "MgO"; "FeO"; "O"; "K2O"; "Na2O"; "TiO2"; "Cr2O3"; "H2O"];
    sys_in  = "wt"
    out     = multi_point_minimization(P,T, data, X=X, B=B, Xoxides=Xoxides, sys_in=sys_in);
    Finalize_MAGEMin(data)

    println("Size of out=$(size(out))")

    T_plot = Vector{Float64}(undef,length(out))
    aTi_plot = Vector{Float64}(undef,length(out)) 

    for idx_out in eachindex(out)
        T_plot[idx_out] = out[idx_out].T_C
        aTi_plot[idx_out] = out[idx_out].aTiO2    
    end
    display(plot(T_plot, aTi_plot, seriestype=:scatter))
    return out
end

out=main()

