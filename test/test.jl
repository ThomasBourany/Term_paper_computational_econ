module modeltest

	using FactCheck, FiscalPolicyModel

	context("Test on accuracy of the algorithm") do

		facts("Are the equilibrium conditions right at the steady state?") do

g = 0.2*ones(20)
tauk=zeros(20)
tauc=zeros(20)
taun=zeros(20)

p=Parameter()
exo=exovariable(g,tauk,tauc,taun,p)
endo=endovariable(exo)
m=Imp_model(exo, endo, p)
model_path(m, 1e-6)

# are the capital and consumption respecting the market clearing conditions at steady state
	@fact m.f(m.endo.k[end], m.p) + (1-m.p["delta"])*m.endo.k[end] -m.endo.k[end] - m.exo.g[end] - m.endo.c[end] --> roughly(0; atol = 1e-5)
	@fact m.deru(m.endo.c[end-1])/m.deru(m.endo.c[end]) - m.p["beta"]*((1+m.exo.tauc[end-1])/(1+m.exo.tauc[end]))*((1-m.exo.tauk[end])(m.derf(k[end]) - m.p["delta"])+1) --> roughly(0; atol = 1e-5)
	@fact m.derf(m.endo.k[end]) - (m.p["rho"]/(1-m.exo.tauk[end])) - m.p["delta"] --> roughly(0; atol = 1e-5)

		end 


	facts("Are the equilibrium conditions right at a random period?") do

n = Int(round(49rand(),0))

g = 0.2*ones(20)
tauk=zeros(20)
tauc=zeros(20)
taun=zeros(20)

p=Parameter()
exo=exovariable(g,tauk,tauc,taun,p)
endo=endovariable(exo)
m=Imp_model(exo, endo, p)
model_path(m, 1e-6)

# are the capital and consumption respecting the market clearing conditions at the period n?
	@fact m.f(m.endo.k[n], m.p) + (1-m.p["delta"])*m.endo.k[n] -m.endo.k[n] - m.exo.g[n] - m.endo.c[n] --> roughly(0; atol = 1e-5)
	@fact m.deru(m.endo.c[n-1])/m.deru(m.endo.c[n]) - m.p["beta"]*((1+m.exo.tauc[n-1])/(1+m.exo.tauc[n]))*((1-m.exo.tauk[n])(m.derf(k[n]) - m.p["delta"])+1) --> roughly(0; atol = 1e-5)
	@fact m.derf(m.endo.k[n]) - (m.p["rho"]/(1-m.exo.tauk[n])) - m.p["delta"] --> roughly(0; atol = 1e-5)


	end

end 