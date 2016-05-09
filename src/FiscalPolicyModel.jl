

module FiscalPolicyModel

export Parameter, exovariable, endovariable, Imp_model, loop_c_k_path, shootingalgorithm, model_path

using Roots, Distributions

function Parameter()
# share of capital in the production function
a = 0.33
#productivity
z = 1.0
# depreciation
d = 0.2
# discount factor
b = 0.95
rho = (1.0/b) - 1.0
# CRRA
sig = 2
# initial gvt exp (calib for US)
g = 0.2
# number of period needed to come back to the Steady State
# Therefore the period 1 is the initial condition, and period S+2 is the steady state
S = 30
## Range of value for endogenous variable
rangelb = 1e-6
rangeub = 100
return Dict("alpha"=>a, "z"=>z, "beta"=>b, "rho"=> rho, "sigma" => sig, "delta"=> d, "S"=>S, "ub"=> rangeub, "lb"=> rangelb)
end

function u(c, p::Dict) ##Utility
### Use the CRRA Utility fct
co=float(c)
sig=float(p["sigma"])
if sig >1.0 
return ((co^(1-sig))-1)/(1-sig)
elseif sig == 1.0
return ln(co)
end
end

function deru(c, p::Dict) ##Utility
### Use the CRRA Utility fct
co=float(c)
sig=float(p["sigma"])
if sig >1.0
return co^(-sig)
elseif sig == 1.0
return 1/co
end
end

function invderu(u, p::Dict)
uo=float(u)
sig=float(p["sigma"])
if sig>1.0
return uo^(-1/sig)
elseif sig == 1.0
return 1/uo
end 
end

function f(k::Real,p::Dict) #Production
### Cobb Douglas
alpha = float(p["alpha"])
ko = float(k)
return ko^alpha
end

function derf(k::Real,p::Dict)
alpha = float(p["alpha"])
ko = float(k)
return alpha*ko^(alpha-1)
end

type Exovar
	## exogenous policy vector
	g::Vector{Real}		#public spending
	tauk::Vector{Real}	#tax on capital
	tauc::Vector{Real}	#tax on consumption
	taun::Vector{Real}	#tax on labor
end

type Endovar
	c::Vector{Real} #conso
	k::Vector{Real}	#capital
	q::Vector{Real}	#pretax price of the good
	eta::Vector{Real}	#pretax profit of HH
	w::Vector{Real}		#pretax wage
	R::Vector{Real}		#gross one-period interest rate
	r::Vector{Real}		#net interest rate
end

function exovariable(var1::Vector,var2::Vector,var3::Vector, var4::Vector, p::Dict)
## Function gathering the fluctuation of the exogenous variable until their last period (var[end])
## Add enough time at the end until the steady state value, (for the endogenous variable to adjust)
var1bar=var1[end]*ones(p["S"])
var2bar=var2[end]*ones(p["S"])
var3bar=var3[end]*ones(p["S"])
var4bar=var4[end]*ones(p["S"])
cvar1 = append!(var1,var1bar)
cvar2 = append!(var2,var2bar)
cvar3 = append!(var3,var3bar)
cvar4 = append!(var4,var4bar)
return Exovar(cvar1, cvar2, cvar3, cvar4)
end



function endovariable(exovar::Exovar)
## generate endogenous variable, empty, to be used latter in the shooting algorithm
n = length(exovar.g)
c = zeros(n)
k = zeros(n)
q = zeros(n)
eta= zeros(n)
w= zeros(n)
R= zeros(n)
r= zeros(n)
return Endovar(c,k,q,eta,w,R,r)
end

type Model
	exo::Exovar
	endo::Endovar
	u::Function
	deru::Function
	invderu::Function
	f::Function
	derf::Function
	p::Dict
	n::Int ## number of period in the model, including the S (adjustment)
end

function Imp_model(exo::Exovar, endo::Endovar, p::Dict)
	n = length(exo.g)
#	u, deru, invderu, f, derf
	return Model(exo, endo, u, deru, invderu, f, derf, p, n)
end

## All the following functions are part of the "final" shooting algorithm. That is why we start 
## from the computed value of the steady state. After that will choose (in a "rule of thumb" manner) the initial value of the control variable
## to approximate this steady state value.  

function steadystate(m::Model)
## We compute the steady state value of the capital stock
taukbar = m.exo.tauk[end]
sol(k) = m.p["delta"] + (m.p["rho"]/(1-taukbar)) - m.derf(k,m.p)
kbar = fzero(k->sol(k),[m.p["lb"],m.p["ub"]])
## The kbar satisfy the equation XX of the paper. 

sol2(k) = m.p["delta"] + m.p["rho"] - m.derf(k,m.p)
kbar2 = fzero(k->sol2(k),[m.p["lb"],m.p["ub"]])
# the kbar2 satisfy the golden rule, without tax. 

## We assign the steady state value to the vector of endogenous variables
m.endo.k[end] = kbar
return kbar, kbar2
end 


function lawofmotion_k(m::Model, c::Vector, k::Vector, i) ##Compute next period capital 
g=m.exo.g
ki= m.f(k[i],m.p) + ((1-m.p["delta"])*k[i]) - g[i] - c[i]
if ki>0
k[i+1]=ki
elseif ki<=0
k[i+1]=0.0
end
return nothing
end


function eulerequation(m::Model, uc::Vector, k::Vector, i) ##Compute next period marg.utility
tauc=m.exo.tauc
tauk=m.exo.tauk
sol(uprime) = (uprime*(m.p["beta"]*((1+tauc[i])/(1+tauc[i+1]))*((1-tauk[i+1])*(m.derf(k[i+1],m.p) - m.p["delta"]) + 1))) - uc[i]

## solving for the next period marg utility
ui=fzeros(u->sol(u), [-Inf, Inf])
ui_bis =[ui;NaN] # in case the euler equation is always infinite and the fzeros gives a 0-element array
uc[i+1]=ui_bis[1]
return nothing
end


function equilibriumquantities(m::Model, c::Vector, k::Vector) ## From the path of consumption and capital, compute the other endogenous variables

## As preliminary assign the newly computed c and k quantities (in fact it does not matter when you do it.)
m.endo.c = c
m.endo.k = k

tc = m.exo.tauc # rename conso-tax for simplicity
tk = m.exo.tauk

#Loop
for i in 1:m.n
## First consumption (even though we already computed it in the loop)
#if i >1
#m.endo.c[i-1] = f(k[i-1],m.p) + (1-m.p["delta"])*(k[i-1]) - k[i] - m.exo.g[i-1]
#end 
#if i=m.n
#m.endo.c[m.n] = f(k[m.n],m.p) + (1-m.p["delta"])*(k[m.n]) - m.exo.g[m.n]
#end

## Second compute the sequence of asset prices
m.endo.q[i]=(m.p["beta"]^i) *deru(c[i], m.p) /(1+tc[i])

## Third compute the rental rate of HH from capital
m.endo.eta[i]= m.derf(m.endo.k[i],m.p)

## Fourth, the wage 
m.endo.w[i]= m.f(k[i],m.p)- k[i] * m.derf(k[i],m.p)

## Fifth the gross interest rate btw t and t+1
if i >1
m.endo.R[i-1]= ((1-tk[i])*(m.derf(k[i],m.p) - m.p["delta"]) +1)
end 

if i==m.n
m.endo.R[m.n]=m.endo.R[m.n-1] 
end

## Sixth the net interest rate
if i >1
m.endo.r[i-1]= m.endo.R[i-1]-1
end 
if i==m.n
m.endo.r[m.n]=m.endo.r[m.n-1] 
end

end 
## setting the last interest rate as equal to the one before the last (at steady-state, we suppose it does not fluctuate much). 

return nothing
end 



function loop_c_k(m::Model, c0::Real, kbar::Real)

k_temp = zeros(m.n) # update the k-vector, to zero for each new trial of the shooting 
k_temp[1]=kbar
c_temp=  zeros(m.n)
c_temp[1] = c0 
uc_temp = zeros(m.n) 
uc_temp[1]=m.deru(c0,m.p)
for j in 1:m.n-1
lawofmotion_k(m, c_temp, k_temp, j)  	##Compute next period capital 
eulerequation(m, uc_temp, k_temp, j)	##Compute next period marg.utility
c_temp[j+1] = invderu(uc_temp[j+1],m.p)		## inverse to get the nextperiod consumption
#uc[i] = m.deru(c[i])
end
ks = k_temp[end]
cs = c_temp[end]

## only return the last period consumption and capital
return ks, cs

end



function loop_c_k_path(m::Model, c0::Real, kbar::Real)
### the same loop, but return the path of the capital and consumption

k_temp = zeros(m.n) # update the k-vector, to zero for each new trial of the shooting 
k_temp[1]=kbar
c_temp=  zeros(m.n)
c_temp[1] = c0 
uc_temp = zeros(m.n) 
uc_temp[1]=m.deru(c0,m.p)
for j in 1:m.n-1
lawofmotion_k(m, c_temp, k_temp, j)  	##Compute next period capital 
eulerequation(m, uc_temp, k_temp, j)	##Compute next period marg.utility
c_temp[j+1] = invderu(uc_temp[j+1],m.p)		## inverse to get the nextperiod consumption
#uc[i] = m.deru(c[i])
end
k_path = k_temp
c_path = c_temp

return k_path, c_path
end



	function shootingalgorithm(m::Model)

s = steadystate(m)
kbarnt = s[2]  ## the steady state without tax as starting point (initial condition for k)
kbar = s[1] ## the steady state with tax for the target of the algorithm 


#################################################
############################ 1st scale

## Draw a grid at the 1e-3 level
c0bis=[1e-3n for n in 1:1000]
kbar_bis = zeros(1000)
cbar_bis = zeros(1000)

	for index in 1:1000
c0=c0bis[index]
# compute the 1000-paths of consumption and capital drawn from the loop of the model. 
result = loop_c_k(m, c0, kbarnt)
kbar_bis[index] = result[1]
cbar_bis[index] = result[2]
	end

### find the index (among the 1000) which gives a final capital close to the steady state value
kdiff = kbar_bis - kbar*ones(1000)
ind_temp=findmin(abs(kdiff))[2] ## find its index
if kdiff[ind_temp] > 0
ind = ind_temp
elseif kdiff[ind_temp] <0
ind = ind_temp-1
end 
##OR
#ind = findmax(cbar_bis)[2] ## find the index where the consumption is not null (slightly faster)
c1 = 1e-3*ind
c1ter=round(c1, 5)
ks = kbar_bis[ind]


#################################################
############################ 2nd scale
#### Draw a grid at the 1e-6 level

c1bis = [c1ter+1e-6n for n in 1:1000]
kbar_ter = zeros(1000)
cbar_ter = zeros(1000)

	for index2 in 1:1000
c1=c1bis[index2]
# again, compute the 1000-paths of consumption and capital drawn from the loop of the model
# (with difference in initial condition much smaller than in the loop above). 
result= loop_c_k(m, c1, kbarnt)
kbar_ter[index2] = result[1]
cbar_ter[index2] = result[2]
	end 

## Find the "shoot" giving the target-capital the closest to the steady-state value. 
kdiff = kbar_ter - kbar*ones(1000)
ind_temp=findmin(abs(kdiff))[2] ## find its index

# you want to take the "shoot", for the value yielding a Ks superior to Kbar
if kdiff[ind_temp] > 0
ind2 = ind_temp
elseif kdiff[ind_temp] <0
ind2 = ind_temp-1
end 
c2ter = c1ter+1e-6*ind2 ## compute the closest initial consumption level
ks= kbar_ter[ind2]

#################################################
############################ 3rd scale
######### Draw a grid at the 1e-9 level

c2bis = [c2ter+1e-9n for n in 1:1000]
kbar_quar = zeros(1000)
cbar_quar = zeros(1000)

	for index3 in 1:1000
c2=c2bis[index3]
# again, compute the 1000-paths of consumption and capital drawn from the loop of the model
# (with difference in initial condition again much smaller than in the loop above). 
result= loop_c_k(m, c2, kbarnt)
kbar_quar[index3] = result[1]
cbar_quar[index3] = result[2]
	end 

## Find the "shoot" giving the target-capital the closest to the steady-state value. 
kdiff = kbar_quar - kbar*ones(1000)
ind_temp=findmin(abs(kdiff))[2] ## find its index

# you want to take the "shoot", for the value yielding a Ks superior to Kbar
if kdiff[ind_temp] > 0
ind3 = ind_temp
elseif kdiff[ind_temp] <0
ind3 = ind_temp-1
end 
c3ter = c2ter+1e-9*ind3  ## compute the closest initial consumption level
ks= kbar_quar[ind3]


#################################################
############################ 4th scale
########## Draw a grid at the 1e-12 level

c3bis = [c3ter+1e-12n for n in 1:1000]
kbar_quin = zeros(1000)
cbar_quin = zeros(1000)

	for index4 in 1:1000
c3=c3bis[index4]
# again, compute the 1000-paths of consumption and capital drawn from the loop of the model
# (with difference in initial condition again much smaller than in the loop above). 
result= loop_c_k(m, c3, kbarnt)
kbar_quin[index4] = result[1]
cbar_quin[index4] = result[2]
	end 

## Find the "shoot" giving the target-capital the closest to the steady-state value. 
kdiff = kbar_quin - kbar*ones(1000)
ind_temp=findmin(abs(kdiff))[2] ## find its index

# you want to take the "shoot", for the value yielding a Ks superior to Kbar
if kdiff[ind_temp] > 0
ind4 = ind_temp
elseif kdiff[ind_temp] <0
ind4 = ind_temp-1
end 
c4ter = c3ter+1e-12*ind4  ## compute the closest initial consumption level
ks= kbar_quin[ind4]

## Results and display a small message after the result

c0 = c4ter
ks = kbar_quin[ind4]

signif=0.0
for i in 1:9
if abs(kbar-ks) < 10.0^(-i)
signif = 10.0^(-i)
end 
end

println("The shooting algorithm needs 4 loops (at 4 different scales (1e-3, 1e-6, 1e-9, 1e-12) to get the right c0 = $c0")
println("The obtained final capital is Ks = $ks, which is close to the true kbar = $kbar at the $signif level")
return c0, ks, kbar, signif
end

### Compute the rest. 
function model_path(m::Model, tol::Real=1e-6)

s = steadystate(m)
kbarnt = s[2]  ## the steady state without tax as starting point (initial condition for k)

shooting = shootingalgorithm(m)

c0 = shooting[1]
ks = shooting[2]
kbar = shooting[3]


	if abs(ks - kbar) < tol 
	# it should be if the loop-above is correct (test it?)
	## assign the new quantities & compute the equilibrium
k_temp = zeros(m.n)
c_temp = zeros(m.n)
path= loop_c_k_path(m, c0, kbarnt)
k_temp = path[1]	
c_temp = path[2]
equilibriumquantities(m, c_temp, k_temp)

println("The path of the endogenous variables has been updated, following the fiscal policy shocks")

	else
	println("The path of the endogenous variables has not been updated, because the shooting algorithm did not find an accurate optimal path")

	end 

return nothing 
#return m 

end



	function runAll()
	#	println("running tests:")
	#	include("test/runtests.jl")
	#	println("")
		println("Result of the shooting algorithm:")
		#shootingalgorithm(m)
		println("")
		ok = input("enter y to close this session.")
		if ok == "y"
			quit()
		end
	end





				end 

