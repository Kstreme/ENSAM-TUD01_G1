using Plots
using DifferentialEquations

q(x)=1
EJ=1

xspan=(0.0 , 1.0)

f(x)=q(x)/EJ



function beam!(du,u,p,x)
    du[1]=u[2]
    du[2]=u[3]
    du[3]=u[4]
    du[4]=f(x)
end



function bc!(r,u,p,x)
    r[1]= u[1][1]
    r[2]= u[3][1]
    r[3]= u[1][end]
    r[4]= u[3][end]
end



bvp= BVProblem(beam!,bc!,[0,0,0,0],xspan)
u_sol= solve(bvp,BS3(),saveat = 0.01)



plot(u_sol.t,-u_sol[1,:],label="displacement")
plot!(u_sol.t,0*u_sol.t,label="beam")
plot!(u_sol.t,0.01*f.(u_sol.t),label="distributed load")