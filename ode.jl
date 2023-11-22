using DifferentialEquations
using Plots 

F=1
function M(x)
    1
end #momentum
EJ=1

function beam!(du,u,p,t)
   
    du[1]= u[2]
    du[2]= -M(x)
end 

function bc!(residual,u,p,x)
    residual[1]=u[1]u[1]
    residual[2]=u[2]u[1]
end

h=0.01
xspan(0,1)
bvp1 = BVProblem(beam!,bc!,[0,0],xspan)
sol1 = solve(bvp1,GeneralMIRK4(),dt=0.005)

plot(M(x))

