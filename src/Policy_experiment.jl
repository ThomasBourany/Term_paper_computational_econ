

module Policy_experiment

include("src/FiscalPolicyModel.jl")
using FiscalPolicyModel
using PyPlot


############################################ Preliminary
############################################
### Compute the steady state value, without policy experiment. 
function modelref()

g = 0.2*ones(20)
tauk=zeros(20)
tauc=zeros(20)
taun=zeros(20)

p=Parameter()
exo=exovariable(g,tauk,tauc,taun,p)
endo=endovariable(exo)
m=Imp_model(exo, endo, p)
#Modify all the variables of the model, following the exogenous variations. 
model_path(m, 1e-6)

mref=Imp_model(m.exo, m.endo, m.p)
return mref

end

############################################
### Plot of convergence of the shooting algorithm, starting with initial condition below the steady state

function graph1(m)
shooting = shootingalgorithm(m)

c0 = shooting[1]
ks = shooting[2]
kbar = shooting[3]

c0_below=[0.1, 0.5, 0.9c0, 0.99c0, (1-1e-3)c0, (1-1e-4)c0, (1-1e-5)c0, (1-5*1e-6)c0, (1-4*1e-6)c0, (1-3*1e-6)c0, (1-2*1e-6)c0,  (1-1e-6)c0, (1-5*1e-7)c0, (1-2*1e-7)c0, (1-1e-7)c0, (1-5*1e-8)c0,(1-4*1e-8)c0, (1-3*1e-8)c0, (1-2*1e-8)c0, (1-1e-8)c0]

# reminder truec0 = 0.64264525131
x=[n for n in 1:m.n] ## periods of time
#plot(x, m.endo.k)

kbarvector = kbar*ones(m.n)
cbarvector = 0.6426452513109608*ones(m.n)

fig, axes = subplots(2,2)
#    ax[:set_ylim](0, 10)
#    ax[:set_xlim](, maximum(grid))
#    lb = "Initial condition"
    jet = ColorMap("jet")[:__call__]
    fig[:suptitle]("Shooting algorithm, convergence from below the steady state level")
    ax = axes[1,1]
    ax[:plot](x, kbarvector, color=jet(0))
    ax[:set_title]("Convergence, capital")
    ax = axes[2,1]
    ax[:plot](x, kbarvector, color=jet(0))
	ax[:set_title]("Capital path, zoom")

	ax = axes[1,2]
    ax[:plot](x, cbarvector, color=jet(0))
    ax[:set_title]("Convergence, consumption")
    ax = axes[2,2]
    ax[:plot](x, cbarvector, color=jet(0))
    ax[:set_title]("Consumption path, zoom")
# lw = 2 lign width
# alpha = transparence
# jet = ColorMap("jet")[:__call__] (color of the jet() 0 if 

	for index in 1:length(c0_below)
c0=c0_below[index]
path= loop_c_k_path(m, c0, kbar)
k_temp = path[1]	
c_temp = path[2]
ax = axes[1,1]
ax[:plot](x, k_temp, color=jet(index/length(c0_below)))

if index > 9 
ax = axes[2,1]
ax[:plot](x, k_temp, color=jet(index/length(c0_below)))
end 

ax = axes[1,2]
ax[:plot](x, c_temp, color=jet(index/length(c0_below)))

if index > 9 
ax = axes[2,2]
ax[:plot](x, c_temp, color=jet(index/length(c0_below)))
end 
	end
#ax[:legend](loc="upper left")
end 


############################################
### ### Plot of convergence, starting with initial condition above the steady state


function graph2(m)

shooting = shootingalgorithm(m)

c0 = shooting[1]
ks = shooting[2]
kbar = shooting[3]

c0_above=[0.95, 0.8, c0*1.1, c0*1.01, (1+1e-3)c0, (1+1e-4)c0, (1+1e-5)c0, (1+5*1e-6)c0, (1+4*1e-6)c0, (1+3*1e-6)c0, (1+2*1e-6)c0,  (1+1e-6)c0, (1+5*1e-7)c0, (1+2*1e-7)c0, (1+1e-7)c0, (1+5*1e-8)c0,(1+4*1e-8)c0, (1+3*1e-8)c0, (1+2*1e-8)c0, (1+1e-8)c0]

# reminder true_c0 = 0.64264525131
x=[n for n in 1:m.n] ## periods of time
#plot(x, m.endo.k)

kbarvector = kbar*ones(m.n)
cbarvector = 0.6426452513109608*ones(m.n)

fig, axes = subplots(2,2)
#    ax[:set_ylim](0, 10)
#    ax[:set_xlim](, maximum(grid))
#    lb = "Initial condition"
    jet = ColorMap("jet")[:__call__]
     fig[:suptitle]("Shooting algorithm, convergence from above the steady state level")
    ax = axes[1,1]
    ax[:plot](x, kbarvector, color=jet(0))
    ax[:set_title]("Convergence, capital")
    ax = axes[2,1]
    ax[:plot](x, kbarvector, color=jet(0))
	ax[:set_title]("Capital path, zoom")

	ax = axes[1,2]
    ax[:plot](x, cbarvector, color=jet(0))
    ax[:set_title]("Convergence, consumption")
    ax = axes[2,2]
    ax[:plot](x, cbarvector, color=jet(0))
    ax[:set_title]("Consumption path, zoom")
# lw = 2 lign width
# alpha = transparence
# jet = ColorMap("jet")[:__call__] (color of the jet() 0 if 

	for index in 1:length(c0_above)
c0=c0_above[index]
path= loop_c_k_path(m, c0, kbar)
k_temp = path[1]	
c_temp = path[2]
ax = axes[1,1]
ax[:plot](x, k_temp, color=jet(index/length(c0_above)))

if index > 9 
ax = axes[2,1]
ax[:plot](x, k_temp, color=jet(index/length(c0_above)))
end 

ax = axes[1,2]
ax[:plot](x, c_temp, color=jet(index/length(c0_above)))

if index > 9 
ax = axes[2,2]
ax[:plot](x, c_temp, color=jet(index/length(c0_above)))
end

	end

end 

#### Draw the graph
graph1(mref)
graph2(mref)


############################################
#### Plot the IRF, for the usual variables 

## Need to feed the function with a model with exogenous shocks and a reference model, without policy shocks. 

function graph3(m, mref, tau =g)

fprod=zeros(m.n)
for j in 1:m.n
fprod[j]=m.f(m.endo.k[j], m.p)
end

fprodref=zeros(m.n)
for j in 1:m.n
fprodref[j]=m.f(mref.endo.k[j], m.p)
end


invest = zeros(m.n)
for j in 1:m.n-1
invest[j]=m.endo.k[j+1] - (1-m.p["delta"])*m.endo.k[j]
end
invest[m.n]=invest[m.n-1]

investref = zeros(m.n)
for j in 1:m.n-1
investref[j]=mref.endo.k[j+1] - (1-mref.p["delta"])*mref.endo.k[j]
end
investref[m.n]=investref[m.n-1]


x=[n for n in 1:m.n]

fig,axes = subplots(3,3)

ax = axes[1,1]
ax[:set_title]("Capital")
ax[:plot](x, m.endo.k)
ax[:plot](x, mref.endo.k, "k--")

ax = axes[1,2]
ax[:set_title]("Consumption")
ax[:plot](x, m.endo.c)
ax[:plot](x, mref.endo.c, "k--")

ax = axes[1,3]
ax[:set_ylim](-0.1, 0.6)
if tau == g
ax[:set_title]("Government spending")
ax[:plot](x, m.exo.g)
ax[:plot](x, mref.exo.g, "k--")
elseif tau == tauc
ax[:set_title]("Consumption tax")
ax[:plot](x, m.exo.tauc)
ax[:plot](x, mref.exo.tauc, "k--")
elseif tau == tauk
ax[:set_title]("Capital tax")
ax[:plot](x, m.exo.tauk)
ax[:plot](x, mref.exo.tauk, "k--")
end 

ax = axes[2,1]
ax[:set_title]("Pretax-price")
ax[:plot](x, m.endo.q)
ax[:plot](x, mref.endo.q, "k--")

ax = axes[2,2]
ax[:set_title]("Gross-interest rate")
ax[:plot](x, m.endo.R)
ax[:plot](x, mref.endo.R, "k--")

ax = axes[2,3]
ax[:set_title]("Rental price of capital")
ax[:plot](x, m.endo.eta)
ax[:plot](x, mref.endo.eta, "k--")


ax = axes[3,1]
ax[:set_title]("Production y")
ax[:plot](x, fprod)
ax[:plot](x, fprodref, "k--")

ax = axes[3,2]
ax[:set_title]("Investment")
ax[:plot](x, invest)
ax[:plot](x, investref, "k--")

ax = axes[3,3]
ax[:set_title]("Wage")
ax[:plot](x, m.endo.w)
ax[:plot](x, mref.endo.w, "k--")



fig[:canvas][:draw]()

end 

############################################ Second part, 
############################################ First policy shock
## Second, increase once-for-all in g at period 10. 

function policy1()

tauk=zeros(20)
tauc=zeros(20)
taun=zeros(20)
g = [0.2*ones(10);0.4*ones(10)]

p=Parameter()
exo=exovariable(g,tauk,tauc,taun,p)
endo=endovariable(exo)
m1=Imp_model(exo, endo, p)
model_path(m1, 1e-6)

graph3(m1, mref, g)

graph1(m1)

end

############################################ Third part, 
############################################ second policy shock
## Third, increase once-for-all the consumption tax at period 10. 

function policy2()

tauk=zeros(20)
tauc=[zeros(10);0.2*ones(10)]
taun=zeros(20)
g = 0.2*ones(20)

p=Parameter()
exo=exovariable(g,tauk,tauc,taun,p)
endo=endovariable(exo)
m2=Imp_model(exo, endo, p)

model_path(m2, 1e-6)

graph3(m2, mref, tauc)

graph1(m2)

end 


############################################ Fourth part, 
############################################ Third policy shock
## Fourth, increase once-for-all the capital tax at period 10. 

function policy3()

tauc=zeros(20)
tauk=[zeros(10);0.2*ones(10)]
taun=zeros(20)
g = 0.2*ones(20)


p=Parameter()
exo=exovariable(g,tauk,tauc,taun,p)
endo=endovariable(exo)
m3=Imp_model(exo, endo, p)
model_path(m3, 1e-6)

graph3(m3, mref, tauk)


end 


############################################ Fifth part, 
############################################ Fourth policy shock
## Fifth, increase, one time-pulse of government spending at period 10

function policy4()

tauk=zeros(20)
tauc=zeros(20)
taun=zeros(20)
g = 0.2*ones(20)
g[10] = 0.4

p=Parameter()
exo=exovariable(g,tauk,tauc,taun,p)
endo=endovariable(exo)
m4=Imp_model(exo, endo, p)
model_path(m4, 1e-6)

graph3(m4, mref, g)

end 

############################################ Sixth part, 
############################################ Fifth policy shock
## White-noisy government spending

function policy5()

g=[0.7rand(19);0.2]

tauk=zeros(20)
tauc=zeros(20)
taun=zeros(20)


p=Parameter()
exo=exovariable(g,tauk,tauc,taun,p)
endo=endovariable(exo)
m5=Imp_model(exo, endo, p)
model_path(m5, 1e-6)

graph3(m5, mref, g)

end 

############################################ Seventh part, 
############################################ Sixth policy shock
## White-noisy consumption tax

function policy6()

tauc=[0.7rand(19);0.2]
g = 0.2*ones(20)
tauk=zeros(20)
taun=zeros(20)


p=Parameter()
exo=exovariable(g,tauk,tauc,taun,p)
endo=endovariable(exo)
m6=Imp_model(exo, endo, p)
model_path(m6, 1e-6)

graph3(m6, mref, tauc)

end


############################################ Eighth part, 
############################################ Seventh & Eighth policy shock
## Concave vs. convexe increasing spending

function policy7()

h1(x)= (0.00366313)*exp(x/5)
h2(x)= 0.2(1-exp(-0.2x) + 0.000915782*x)
g1 = zeros(20)
g2 = zeros(20)

for i in 1:20
g1[i] = h1(i)
g2[i] = h2(i)
end 
g1[end] = round(g1[end], 5)
g2[end] = round(g2[end], 5)


tauk=zeros(20)
tauc=zeros(20)
taun=zeros(20)

# 1# convex increase 

p=Parameter()
exo1=exovariable(g1,tauk,tauc,taun,p)
endo1=endovariable(exo1)
m7=Imp_model(exo1, endo1, p)
model_path(m7, 1e-6)
g=g1
graph3(m7, mref, g)


# 2# concave increase

p=Parameter()
exo2=exovariable(g2,tauk,tauc,taun,p)
endo2=endovariable(exo2)
m8=Imp_model(exo2, endo2, p)
model_path(m8, 1e-6)
g = g2
graph3(m8, mref, g)

end 






    function runAll()
    
modelref()
policy1()
policy2()
policy3()
policy4()
policy5()
policy6()
policy7()

        ok = input("enter y to close this session.")
        if ok == "y"
            quit()
        end
    end

end 




