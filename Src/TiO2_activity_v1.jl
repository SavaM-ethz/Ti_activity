using Plots, MAGEMin_C, Pkg, JLD2
Pkg.activate(".")

# ------------------------- #
#|  Include functionality  |#
# ------------------------- #
include("../Funcs/thermodynamic_functionality.jl")

# ------------------------- #
#|  Define Main function   |#
# ------------------------- #
function main()

    # Switches
    call_MAGEMin        = false      # Calculate titanium activity
    postprocess_MAGEMin = true     # Postprocess data

    # Initialization
    out      = []
    nT       = 50
    T_plot   = Vector{Float64}(undef,nT)
    aTi_plot = Vector{Float64}(undef,nT)

    # Calculate TiO2-activity
    if call_MAGEMin

        # Prepare & run MAGEMin calculation
        db      = "ig"  # database: ig, igneous (Holland et al., 2018); mp, metapelite (White et al 2014b)
        data    = Initialize_MAGEMin(db, buffer="qfm", verbose=false);
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

        # Save results
        JLD2.save("./Data/MAGEMin_out_struct.jld2", "out", out)
    end

    # Postprocess data structure
    if postprocess_MAGEMin
        out = JLD2.load("./Data/MAGEMin_out_struct.jld2", "out")
        println(size(out))
        for idx_out in eachindex(out)
            T_plot[idx_out]   = out[idx_out].T_C
            aTi_plot[idx_out] = out[idx_out].aTiO2
        end

        # Visualize
        display(plot(T_plot, aTi_plot, seriestype = :scatter))
    end

    # Return
    return nothing

end

out=main()

