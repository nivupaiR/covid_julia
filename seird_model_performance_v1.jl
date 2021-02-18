### A Pluto.jl notebook ###
# v0.12.20

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : missing
        el
    end
end

# ╔═╡ e0c49314-5cda-11eb-2c98-3db586f72f4c
using VegaLite, VegaDatasets, DataFrames, CSV, Dates,Statistics,PlutoUI,Markdown,Gadfly,DifferentialEquations,Random,Distributions,PlutoUI,Query,Plots,Printf

# ╔═╡ 267a42fe-5cfa-11eb-2b77-15f408b05b41
md"	We can model the impact of quarantine and other measures which can delay the pandemic spread. We can also model the delay. We call this SEIQR model. q: range 0-1. 0 means no quarantine. 1=> quarantine. Value in between is a ratio. As soon as people are exposed, a certain fraction q of them go into quarantine and thus avoid further contact with susceptible people. "

# ╔═╡ 79b8e40c-5cfc-11eb-21a8-db8981151d30
md"	The death rate is modeled as μ=μ0+μSlope*t. When μSlope ≠ 0, the death rate changes linearly wirht time. In many situations, the death rate exhibits a time varying behavior. In the short term, it can be modeled as a linear function of time. If the death rate is constant, μSlope=0 and μ0 is the fixed value. "

# ╔═╡ ee41813a-6fe4-11eb-026e-affcce424903
md"	 μSlope=0 for fixed mortality rate. Set this to a nonzero value if the mortality rate is changing linearly with time. μ=μ0+μSlope*t"

# ╔═╡ 92708cfa-5ce5-11eb-101b-0b41cbb434de
@bind μSlope Slider(0:0.001:0.1; default=0, show_value=true)

# ╔═╡ 7a6cf7d4-5ce4-11eb-0941-d16967bac9fb
@bind μ0 Slider(0:0.01:0.1; default=0.02, show_value=true)

# ╔═╡ 9bd62b20-5db1-11eb-3239-1d5f8792ca17
begin
	function newSEIQRD(dW,W,Z,p,t)
		S,E,I,R,D,C,A,Q,J = W
		β,γ,σ,μ0,μSlope,tau,q,m = p

		# @show β,γ,σ,μ0,μSlope,tau,q,m,length(p)
  hDelay = Z(p, t-tau)[3]         # I=U[3]
	qDelay = Z(p, t-tau)[8]	
	mDelay = Z(p, t-tau)[9]	
		
	dW[1] = dS = -β*I*S           # Susceptible
    dW[2] = dE =  β*I*S - (1-q)*σ*E-q*E      # Exposed, q for quarantine factor
    dW[3] = dI = (1-m)*(1-q)*σ*E -γ*hDelay-A*hDelay   # Infected 
	dW[4] = dR = γ*(1-A)*hDelay + γ*(1-A)*mDelay +γ*(1-A)*qDelay    # Recovered
	dW[5] = dD = A*hDelay + A*mDelay +A*qDelay    # Dead
	dW[6] = dC = (1-q)*σ*E         # Confirmed cases
	dW[7] = dA = μSlope 		   # CFR or death rate change
	dW[8] = dQ = q*σ*E -γ*qDelay-A*qDelay            # Exposed->+Q=>+I
	dW[9] = dJ = m*(1-q)*σ*E-γ*mDelay-A*mDelay       # Infectious but mandatory quarantined
		@show A
end
	
	maxTimeSpan1=200
	Z(p, t) = ones(9)  #The size should be >=k where Z(p, t-tau)[k]	
	@show Z(4,5)
	tauQ = 1.2
	lagsQ = [tauQ]
	tspanQ = (0.0,maxTimeSpan1)
	W0 = [1000,0.0,2.0,0.0,0.0,0.0,μ0,0.0,0.0];  # S,E,I,R,D,C,A,Q,J
	function seiqrd_new(x)
	    prob = DDEProblem(newSEIQRD,W0,Z,tspanQ,(x[1],x[2],x[3],x[4],x[5],x[6],x[7],x[8]); constant_lags=lagsQ)
	    sol = solve(prob)
	    return sol
	end
end

# ╔═╡ 22530dc2-5cdd-11eb-2e70-c3de5b46c77f
begin
	function delaySEIQR(dU,U,h,p,t)
		S,E,I,R,D,C,A,Q = U
		β,γ,σ,μ0,μSlope,tau,q = p
		# @show q,tau,length(p)
  hDelay = h(p, t-tau)[3]         # I=U[3]
		
	dU[1] = dS = -β*I*S           # Susceptible
    dU[2] = dE =  β*I*S - σ*E      # Exposed, q for quarantine factor
    dU[3] = dI = (1-q)*σ*E -γ*I    # Infected 
	dU[4] = dR = (1-A)*γ*(hDelay+q*σ*E)   # Recovered
	dU[5] = dD = A*γ*(hDelay+q*σ*E)       # Dead
		# dU[5] = dD = A*γ*I            # Dead
	dU[6] = dC = (1-q)*σ*E        # Confirmed cases
	dU[7] = dA = μSlope+0*hDelay^2 # CFR or death rate change
	dU[8] = dQ = q*σ*E
end
	
	maxTimeSpan=200
	h(p, t) = ones(3)
	tau = 1.2
	lags = [tau]
	ttspan = (0.0,maxTimeSpan)
	U0 = [1000,0.0,2.0,0.0,0.0,0.0,μ0,0];  # S,E,I,R,D,C,A,Q
	function seirqd_manual(x)
	    prob = DDEProblem(delaySEIQR,U0,h,ttspan,(x[1],x[2],x[3],x[4],x[5],x[6],x[7]); constant_lags=lags)
	    sol = solve(prob)
	    return sol
	end

end

# ╔═╡ d23af0a6-5cda-11eb-2059-13640127bf2c

begin
	function SEIR!(dV,V,p,t)
    S,E,I,R,D,C,A = V
    β,γ,σ,μ0,μSlope = p
    dV[1] = dS = -β*I*S         # Susceptible
    dV[2] = dE = β*I*S - σ*E    # Exposed 
    dV[3] = dI = σ*E -γ*I       # Infected 
	dV[4] = dR = (1-A)*γ*I      # Recovered
	dV[5] = dD = A*γ*I          # Dead
	dV[6] = dC = σ*E            # Confirmed cases
	dV[7] = dA = μSlope         # CFR or death rate change
end
	
	tmax = 500
	tspan = (0.0,tmax)
	V0 = [1000,0,2.0,0.0,0.0,0.0,μ0,0];  # S,E,I,R,D,C,A
	function seird_manual(x)
	    prob = ODEProblem(SEIR!,V0,tspan,(x[1],x[2],x[3],x[4],x[5]))
	    sol = solve(prob)
	    return sol
	end
end

# ╔═╡ 2036ae72-6fe5-11eb-06d6-c58059631b16
md"	β is the rate at which a person gets exposed to the virus upon contact."

# ╔═╡ a77d02d4-5cdd-11eb-0d5e-43e087e56f89
@bind β Slider(0:0.0001:0.5; default=0.002, show_value=true)

# ╔═╡ 41b1535e-6fe5-11eb-3a80-c798261766cc


# ╔═╡ c2b6ebb4-5cdd-11eb-0e62-2bfdeb982445
@bind σ Slider(0:0.0001:0.2; default=0.1, show_value=true)

# ╔═╡ c7399308-5cdd-11eb-38cc-a57f64eddb1f
@bind γ Slider(0:0.001:0.5; default=0.01, show_value=true)

# ╔═╡ d90989e2-5cf8-11eb-169a-2740edaa59a5
@bind qFactor Slider(0:0.01:0.99; default=0.5, show_value=true)

# ╔═╡ e530ec56-5cf8-11eb-1b48-5563833340e4
@bind τ Slider(0:1:28; default=7, show_value=true)

# ╔═╡ a840895e-5de9-11eb-2aae-e122c789ecc4
@bind mFactor Slider(0:0.01:0.99; default=0.5, show_value=true)

# ╔═╡ ff6b75f2-6fe8-11eb-0204-4b5fa1862cb7
qV,qF=qFactor,mFactor

# ╔═╡ e25fd174-6fe8-11eb-2341-2f24df1dd367
β,σ,γ,μ0,qF,qV,τ

# ╔═╡ 6a53f356-5ced-11eb-2590-732b3694f58c
begin
	sol_man=seird_manual([β,γ,σ,μ0,μSlope]);
	dF_SEIRD=DataFrame(Time=sol_man.t,S=sol_man[1,:],E=sol_man[2,:],I=sol_man[3,:],R=sol_man[4,:],D=sol_man[5,:],C=sol_man[6,:]);
	l1=Gadfly.layer(dF_SEIRD,x=:Time,y=:S,Geom.line,color=[colorant"black"],Theme(line_width=1pt))
	l2=Gadfly.layer(dF_SEIRD,x=:Time,y=:E,Geom.line,color=[colorant"deepskyblue"],Theme(line_width=1pt))
	l3=Gadfly.layer(dF_SEIRD,x=:Time,y=:I,Geom.line,color=[colorant"tomato"],Theme(line_width=1pt))
	l4=Gadfly.layer(dF_SEIRD,x=:Time,y=:R,Geom.line,color=[colorant"darkgreen"],Theme(line_width=1pt))
	l5=Gadfly.layer(dF_SEIRD,x=:Time,y=:D,Geom.line,color=[colorant"purple"],Theme(line_width=1pt))
	l6=Gadfly.layer(dF_SEIRD,x=:Time,y=:C,Geom.line,color=[colorant"magenta"],Theme(line_width=1pt))
	 sPlot=Gadfly.plot(l1,l2,l3,l4,l5,l6,Theme(point_size=0.5mm),Guide.title("SEIRD"),Guide.ylabel("Number of people"),Coord.cartesian(xmax=100,ymax=maximum(V0[1])),Guide.manual_color_key("Model", ["model:S","model:E","model:I","model:R","model:D","model:C"],["black","deepskyblue","tomato","darkgreen","purple","magenta"]))
end

# ╔═╡ e9ca9140-5d3e-11eb-3752-83e5f2d01b18
begin
	delay_solX=seiqrd_new([β,γ,σ,μ0,μSlope,τ,qFactor,mFactor])
 dF_SEIQRD_new=DataFrame(Time=delay_solX.t,S=delay_solX[1,:],E=delay_solX[2,:],I=delay_solX[3,:],R=delay_solX[4,:],D=delay_solX[5,:],C=delay_solX[6,:],A=delay_solX[7,:],Q=delay_solX[8,:],J=delay_solX[9,:]);
	
	nl1=Gadfly.layer(dF_SEIQRD_new,x=:Time,y=:S,Geom.line,color=[colorant"black"],Theme(line_width=2pt))
	nl2=Gadfly.layer(dF_SEIQRD_new,x=:Time,y=:E,Geom.line,color=[colorant"deepskyblue"],Theme(line_width=2pt))
	nl3=Gadfly.layer(dF_SEIQRD_new,x=:Time,y=:I,Geom.line,color=[colorant"tomato"],Theme(line_width=2pt))
	nl4=Gadfly.layer(dF_SEIQRD_new,x=:Time,y=:R,Geom.line,color=[colorant"darkgreen"],Theme(line_width=2pt))
	nl5=Gadfly.layer(dF_SEIQRD_new,x=:Time,y=:D,Geom.line,color=[colorant"purple"],Theme(line_width=2pt))
	nl6=Gadfly.layer(dF_SEIQRD_new,x=:Time,y=:C,Geom.line,color=[colorant"magenta"],Theme(line_width=2pt))
	nl8=Gadfly.layer(dF_SEIQRD_new,x=:Time,y=:Q,Geom.line,color=[colorant"gold"],Theme(line_width=2pt))
	nl9=Gadfly.layer(dF_SEIQRD_new,x=:Time,y=:J,Geom.line,color=[colorant"lime"],Theme(line_width=2pt))
	
	Gadfly.plot(l1,l2,l3,l4,l5,l6,nl1,nl2,nl3,nl4,nl5,nl6,nl8,nl9,Theme(point_size=2.5mm),Guide.title("SEIRD"),Guide.ylabel("Number of people"),Coord.cartesian(xmax=200,ymax=1000+0*maximum(V0[1])),Guide.manual_color_key("Model", ["model:S","model:E","model:I","model:R","model:D","model:C"],["black","deepskyblue","tomato","darkgreen","purple","magenta"]))
	
end

# ╔═╡ e0c454aa-5cda-11eb-285b-a3152259fc64
begin
	delay_sol=seirqd_manual([β,γ,σ,μ0,μSlope,τ,qFactor]);
	dF_SEIRD_delay=DataFrame(Time=delay_sol.t,S=delay_sol[1,:],E=delay_sol[2,:],I=delay_sol[3,:],R=delay_sol[4,:],D=delay_sol[5,:],C=delay_sol[6,:],A=delay_sol[7,:],Q=delay_sol[8,:]);
	dl1=Gadfly.layer(dF_SEIRD_delay,x=:Time,y=:S,Geom.line,color=[colorant"black"],Theme(line_width=2pt))
	dl2=Gadfly.layer(dF_SEIRD_delay,x=:Time,y=:E,Geom.line,color=[colorant"deepskyblue"],Theme(line_width=2pt))
	dl3=Gadfly.layer(dF_SEIRD_delay,x=:Time,y=:I,Geom.line,color=[colorant"tomato"],Theme(line_width=2pt))
	dl4=Gadfly.layer(dF_SEIRD_delay,x=:Time,y=:R,Geom.line,color=[colorant"darkgreen"],Theme(line_width=2pt))
	dl5=Gadfly.layer(dF_SEIRD_delay,x=:Time,y=:D,Geom.line,color=[colorant"purple"],Theme(line_width=2pt))
	dl6=Gadfly.layer(dF_SEIRD_delay,x=:Time,y=:C,Geom.line,color=[colorant"magenta"],Theme(line_width=2pt))
	dl8=Gadfly.layer(dF_SEIRD_delay,x=:Time,y=:Q,Geom.line,color=[colorant"gold"],Theme(line_width=2pt))
	 Gadfly.plot(l1,l2,l3,l4,l5,l6,dl1,dl2,dl3,dl4,dl5,dl6,dl8,Theme(point_size=0.5mm),Guide.title("SEIRD"),Guide.ylabel("Number of people"),Coord.cartesian(xmax=100,ymax=maximum(2*V0[1])),Guide.manual_color_key("Model", ["model:S","model:E","model:I","model:R","model:D","model:C"],["black","deepskyblue","tomato","darkgreen","purple","magenta"]))
end

# ╔═╡ f5df82a8-6fc0-11eb-1218-2569b589f227
Gadfly.plot(l1,l3,l4,Theme(point_size=0.5mm),Guide.title("SIR"),Guide.ylabel("Number of people"),Coord.cartesian(xmax=500,ymax=maximum(V0[1])),Guide.manual_color_key("Model", ["model:S","model:I","model:R"],["black","tomato","darkgreen"]))

# ╔═╡ Cell order:
# ╠═e0c49314-5cda-11eb-2c98-3db586f72f4c
# ╠═267a42fe-5cfa-11eb-2b77-15f408b05b41
# ╠═9bd62b20-5db1-11eb-3239-1d5f8792ca17
# ╠═22530dc2-5cdd-11eb-2e70-c3de5b46c77f
# ╠═d23af0a6-5cda-11eb-2059-13640127bf2c
# ╠═79b8e40c-5cfc-11eb-21a8-db8981151d30
# ╠═ee41813a-6fe4-11eb-026e-affcce424903
# ╠═92708cfa-5ce5-11eb-101b-0b41cbb434de
# ╠═7a6cf7d4-5ce4-11eb-0941-d16967bac9fb
# ╠═2036ae72-6fe5-11eb-06d6-c58059631b16
# ╠═a77d02d4-5cdd-11eb-0d5e-43e087e56f89
# ╠═41b1535e-6fe5-11eb-3a80-c798261766cc
# ╠═c2b6ebb4-5cdd-11eb-0e62-2bfdeb982445
# ╠═c7399308-5cdd-11eb-38cc-a57f64eddb1f
# ╠═d90989e2-5cf8-11eb-169a-2740edaa59a5
# ╠═e530ec56-5cf8-11eb-1b48-5563833340e4
# ╠═a840895e-5de9-11eb-2aae-e122c789ecc4
# ╠═ff6b75f2-6fe8-11eb-0204-4b5fa1862cb7
# ╠═e25fd174-6fe8-11eb-2341-2f24df1dd367
# ╠═e9ca9140-5d3e-11eb-3752-83e5f2d01b18
# ╠═e0c454aa-5cda-11eb-285b-a3152259fc64
# ╠═6a53f356-5ced-11eb-2590-732b3694f58c
# ╠═f5df82a8-6fc0-11eb-1218-2569b589f227
