using DifferentialEquations
using Plots 
l=1
F=1
EJ(x)=1+0*x

function M(x)
    -F/(EJ(x))*(l-x) #momentum
end 

function beam!(du,u,p,x)   
    du[1]= u[2]
    du[2]= M(x)
end 

function bc!(residual,u,p,t)
    residual[1]=u[1]u[1]
    residual[2]=u[2]u[1]
end

xspan=(0,1)
bvp = BVProblem(beam!,bc!,[0,0],xspan)
u_sol = solve(bvp,BS3(),saveat=0.01)

plot(u_sol[1,:],label="deflaction")
plot!(xspan,M(xspan),label="momentum")

