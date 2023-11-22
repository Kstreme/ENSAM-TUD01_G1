using Plots
using DifferentialEquations

E(x)=1                 #young module equation
J(x)= 1.1-x                      #inertia equation
M(x)=-(1-x);                  #momentum equation


f(x)=-M(x)/(E(x)*J(x))

function poisson!(du,u,p,x)
    du[1]=u[2]
    du[2]=f(x)
end

function bc!(r,u,p,t)
    r[1] = u[1][1]
    r[2] = u[1][2]
end

xspan=(0.0 , 1.0)


bvp= BVProblem(poisson!,bc!,[0,0],xspan)
u_sol= solve(bvp,BS3(),saveat = 0.01)

plot(u_sol.t,-u_sol[1,:],label="displacement")
plot!(u_sol.t,0*u_sol.t,label="beam")
plot!(u_sol.t,-0.05*M.(u_sol.t),label="momentum")