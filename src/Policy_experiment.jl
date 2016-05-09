

module Policy_experiment

############################################ Preliminary
############################################
### Try the normal steady state value, without policy experiment. 
include("src/FiscalPolicyModel.jl")
using FiscalPolicyModel

g = 0.2*ones(20)
tauk=zeros(20)
tauc=zeros(20)
taun=zeros(20)

p=Parameter()
exo=exovariable(g,tauk,tauc,taun,p)
endo=endovariable(exo)
m=Imp_model(exo, endo, p)
#Modify all the variables of the model, following the exogenous variations. 
model_path(m, 1e-4)

mref=Imp_model(m.exo, m.endo, m.p)


############################################
### Plot of convergence of the shooting algorithm
shooting = shootingalgorithm(m)

c0 = shooting[1]
ks = shooting[2]
kbar = shooting[3]

c0_below=[0.1, 0.5, 0.6, 0.64, (1-1e-3)c0, (1-1e-4)c0, (1-1e-5)c0, (1-5*1e-6)c0, (1-4*1e-6)c0, (1-3*1e-6)c0, (1-2*1e-6)c0,  (1-1e-6)c0, (1-5*1e-7)c0, (1-2*1e-7)c0, (1-1e-7)c0, (1-5*1e-8)c0,(1-4*1e-8)c0, (1-3*1e-8)c0, (1-2*1e-8)c0, (1-1e-8)c0]

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

############################################
### From above
c0_above=[0.95, 0.8, 0.7, 0.65, (1+1e-3)c0, (1+1e-4)c0, (1+1e-5)c0, (1+5*1e-6)c0, (1+4*1e-6)c0, (1+3*1e-6)c0, (1+2*1e-6)c0,  (1+1e-6)c0, (1+5*1e-7)c0, (1+2*1e-7)c0, (1+1e-7)c0, (1+5*1e-8)c0,(1+4*1e-8)c0, (1+3*1e-8)c0, (1+2*1e-8)c0, (1+1e-8)c0]

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


############################################
#### Plot the IRF, for the usual variables 


#Modify all the variables of the model, following the exogenous variations. 


fprod=zeros(m.n)
for j in 1:m.n
fprod[j]=mref.f(m.endo.k[j], m.p)
end

x=[n for n in 1:m.n]
plot(x, mref.endo.k)
plot(x, mref.endo.c)
plot(x, mref.exo.g)
plot(x, fprod)


fig,axes = subplots(2,3,figsize=(10,5))

ax = axes[1,1]
ax[:plot](x, m.endo.k)
ax[:plot](x, zeros(m.n))
ax = axes[1,2]
ax[:plot](x, m.endo.c)
ax[:plot](x, zeros(m.n))
ax = axes[1,3]
ax[:plot](x, m.exo.g)
ax[:plot](x, zeros(m.n))
ax = axes[2,1]
ax[:plot](x, fprod)
ax[:plot](x, zeros(m.n))
ax = axes[2,2]
ax[:plot](x, m.endo.R)
ax[:plot](x, zeros(m.n))
ax = axes[2,3]
ax[:plot](x, m.endo.eta)
ax[:plot](x, zeros(m.n))

fig[:canvas][:draw]()



############################################ Second part, 
############################################ First policy shock
## Second, increase once-for-all in g at period 10. 

tauk=zeros(20)
tauc=zeros(20)
taun=zeros(20)
g = [0.2*ones(10);0.4*ones(10)]

p=Parameter()
exo=exovariable(g,tauk,tauc,taun,p)
endo=endovariable(exo)
m=Imp_model(exo, endo, u, deru, invderu, f, derf, p)
model_path(m, 1e-4)



x=[n for n in 1:m.n]
plot(x, m.endo.k)
plot(x, m.endo.c)
plot(x, m.exo.g)



fig,axes = subplots(2,3,figsize=(10,5))

ax = axes[1,1]
ax[:plot](x, m.endo.k)
#ax[:plot](x, zeros(m.n))
ax = axes[1,2]
ax[:plot](x, m.endo.c)
#ax[:plot](x, zeros(m.n))
ax = axes[1,3]
ax[:plot](x, m.exo.g)
#ax[:plot](x, zeros(m.n))
ax = axes[2,1]
ax[:plot](x, fprod)
#ax[:plot](x, zeros(m.n))
ax = axes[2,2]
ax[:plot](x, m.endo.R)
#ax[:plot](x, zeros(m.n))
ax = axes[2,3]
ax[:plot](x, m.endo.eta)
#ax[:plot](x, zeros(m.n))

fig[:canvas][:draw]()





############################################ Third part, 
############################################ second policy shock
## Third, increase once-for-all the consumption tax at period 10. 



tauk=zeros(20)
tauc=[zeros(10);0.2*ones(10)]
taun=zeros(20)
g = 0.2*ones(20)


p=Parameter()
exo=exovariable(g,tauk,tauc,taun,p)
endo=endovariable(exo)
m=Imp_model(exo, endo, u, deru, invderu, f, derf, p)

model_path(m, 1e-4)


fig,axes = subplots(2,3,figsize=(10,5))

ax = axes[1,1]
ax[:plot](x, m.endo.k)
#ax[:plot](x, zeros(m.n))
ax = axes[1,2]
ax[:plot](x, m.endo.c)
#ax[:plot](x, zeros(m.n))
ax = axes[1,3]
ax[:plot](x, m.exo.tauc)
#ax[:plot](x, zeros(m.n))
ax = axes[2,1]
ax[:plot](x, fprod)
#ax[:plot](x, zeros(m.n))
ax = axes[2,2]
ax[:plot](x, m.endo.R)
#ax[:plot](x, zeros(m.n))
ax = axes[2,3]
ax[:plot](x, m.endo.eta)
#ax[:plot](x, zeros(m.n))

fig[:canvas][:draw]()






############################################ Fourth part, 
############################################ Third policy shock
## Fourth, increase once-for-all the capital tax at period 10. 

tauc=zeros(20)
tauk=[zeros(10);0.2*ones(10)]
taun=zeros(20)
g = 0.2*ones(20)


p=Parameter()
exo=exovariable(g,tauk,tauc,taun,p)
endo=endovariable(exo)
m=Imp_model(exo, endo, u, deru, invderu, f, derf, p)

model_path(m, 1e-4)


fig,axes = subplots(2,3,figsize=(10,5))

ax = axes[1,1]
ax[:plot](x, m.endo.k)
#ax[:plot](x, zeros(m.n))
ax = axes[1,2]
ax[:plot](x, m.endo.c)
#ax[:plot](x, zeros(m.n))
ax = axes[1,3]
ax[:plot](x, m.exo.tauc)
#ax[:plot](x, zeros(m.n))
ax = axes[2,1]
ax[:plot](x, fprod)
#ax[:plot](x, zeros(m.n))
ax = axes[2,2]
ax[:plot](x, m.endo.R)
#ax[:plot](x, zeros(m.n))
ax = axes[2,3]
ax[:plot](x, m.endo.eta)
#ax[:plot](x, zeros(m.n))

fig[:canvas][:draw]()






############################################ Fifth part, 
############################################ Fourth policy shock
## Fifth, increase, one time-pulse of government spending at period 10


tauk=zeros(20)
tauc=zeros(20)
taun=zeros(20)
g = 0.2*ones(20)
g[10] = 0.4

p=Parameter()
exo=exovariable(g,tauk,tauc,taun,p)
endo=endovariable(exo)
m=Imp_model(exo, endo, u, deru, invderu, f, derf, p)
model_path(m, 1e-4)



x=[n for n in 1:m.n]
plot(x, m.endo.k)
plot(x, m.endo.c)
plot(x, m.exo.g)



fig,axes = subplots(2,3,figsize=(10,5))

ax = axes[1,1]
ax[:plot](x, m.endo.k)
#ax[:plot](x, zeros(m.n))
ax = axes[1,2]
ax[:plot](x, m.endo.c)
#ax[:plot](x, zeros(m.n))
ax = axes[1,3]
ax[:plot](x, m.exo.g)
#ax[:plot](x, zeros(m.n))
ax = axes[2,1]
ax[:plot](x, fprod)
#ax[:plot](x, zeros(m.n))
ax = axes[2,2]
ax[:plot](x, m.endo.R)
#ax[:plot](x, zeros(m.n))
ax = axes[2,3]
ax[:plot](x, m.endo.eta)
#ax[:plot](x, zeros(m.n))

fig[:canvas][:draw]()















#### other drafts, to be cleaned

y=[0.001n for n in 1:1000]



# look at the path of capital


x=[n for n in 1:m.n]

fk=zeros(m.n)
for j in 1:m.n
fk[j]=m.derf(k_temp[j], m.p)
end



c_temp3 = m.f(kbar, m.p) + (1-m.p["delta"])*kbar -kbar - g[end]

c_temp2 = zeros(m.n)
for j in 1:m.n-1
c_temp2[j]=m.f(k_temp[j], m.p) + (1-m.p["delta"])*k_temp[j] -k_temp[j+1] -g[j]
end

invest = zeros(m.n)
for j in 1:m.n-1
invest[j]=k_temp[j+1] - (1-m.p["delta"])*k_temp[j]
end



#fig,axes = subplots(1,4,figsize=(10,5))
#ax = axes[1,1]
plot(x, k_temp)
plot(x, c_temp)
#plot(x, c_temp2)
plot(x, uc_temp)
plot(x, fk)


plot(x, fkap)
plot(x, g)
plot(x, invest)
plot(x, c_temp)


ax = axes[1,2]
ax[:plot](x, c_temp)
ax = axes[1,3]
axes[:plot](x, uc_temp)

ax = axes[1,4]
axes[:plot](x, fk)


end 
