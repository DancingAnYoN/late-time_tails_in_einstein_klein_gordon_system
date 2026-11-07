# Smoke tests with synthetic fits: can run without the MMIEF package or GLMakie.
include(joinpath(@__DIR__, "..", "scripts", "interpret_qnm_results.jl"))

function write_fixture(dir::String; bad=false)
    open(joinpath(dir,"config.toml"),"w") do io
        println(io,"M0 = 1.0")
    end
    open(joinpath(dir,"diagnostics.csv"),"w") do io
        println(io,"t,T_asymptotic,M_h")
        for t in 0:1:150
            println(io,"$t,$t,1.0")
        end
    end
    open(joinpath(dir,"qnm_comparison.csv"),"w") do io
        println(io,"J_over_M0,Tlo_over_M0,Thi_over_M0,Mf_omegaR,Mf_gamma,omegaR_error_pct,gamma_error_pct,relative_fit_residual")
        for J in (10.0,30.0,50.0), (i,lo) in enumerate((30.0,35.0,40.0,45.0,50.0))
            w = bad ? 0.145 : 0.110455 * (1.0 + 0.001 * (i-3))
            g = bad ? 0.155 : 0.104896 * (1.0 + 0.001 * (3-i))
            resid = bad ? 0.5 : 0.03
            println(io,"$J,$lo,$(lo+65.0),$w,$g,0,0,$resid")
        end
    end
end

mktempdir() do dir
    write_fixture(dir)
    result = QNMInterpretation.interpret_qnm(dir)
    @assert isfile(result.report)
    @assert isfile(result.summary)
    @assert result.observers == 3
    goodreport = read(result.report,String)
    @assert occursin("Plausible candidate (not confirmed)",goodreport)
    @assert occursin("0.110455",goodreport)
    @assert occursin("0.104896",goodreport)
end

mktempdir() do dir
    write_fixture(dir;bad=true)
    result = QNMInterpretation.interpret_qnm(dir)
    badreport = read(result.report,String)
    @assert occursin("Inconclusive / needs inspection",badreport)
    @assert occursin("real frequency differs",badreport)
end

println("Synthetic interpretation self-tests passed")
