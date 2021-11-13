using NLPModelsJuMP
using JuMP

function Jason1()
 P1=JuMP.Model()
 JuMP.@variable(P1,x[1:2],start=1.0)
 JuMP.@constraint(P1,x[1]^2+x[2]^2-1<=0)
 JuMP.@constraint(P1,x[1]>=0)
 JuMP.@constraint(P1,x[2]>=0)
 JuMP.@NLobjective(P1,Min,0.5*(x[1]-2/3)^2)
 Player1=MathProgNLPModel(P1)

 P2=JuMP.Model()
 JuMP.@variable(P2,x[1:2],start=1.0)
 JuMP.@constraint(P2,x[2]+x[1]-1<=0)
 JuMP.@constraint(P2,x[1]>=0)
 JuMP.@constraint(P2,x[2]>=0)
 JuMP.@NLobjective(P2,Min,0.5*(x[2]-1/3)^2)
 Player2=MathProgNLPModel(P2)
 sol=x->(x[1]==2/3 && x[2]==1/3)

 return GNEP(2,[1,1],[Player1, Player2], sol=sol)
end

function Jason2()
	P1=JuMP.Model()
	JuMP.@variable(P1,x[1:3],start=1.0)
	JuMP.@constraint(P1,x[1]+x[2]+x[3]-3<=0)
	JuMP.@constraint(P1,x[1]>=0)
	JuMP.@constraint(P1,x[2]>=0)
	JuMP.@constraint(P1,x[3]>=0)
	JuMP.@NLobjective(P1,Min,0.5*(x[1]-1)^2+0.5*(x[2]-1)^2)
	Player1=MathProgNLPModel(P1)

	P2=JuMP.Model()
	JuMP.@variable(P2,x[1:3],start=1.0)
	JuMP.@constraint(P2,-x[1]-x[2]-x[3]+3<=0)
	JuMP.@constraint(P2,x[1]^2+x[2]^2+x[3]^2-9<=0)
	JuMP.@constraint(P2,x[1]>=0)
	JuMP.@constraint(P2,x[2]>=0)
	JuMP.@constraint(P2,x[3]>=0)
	JuMP.@NLobjective(P2,Min,0.5*x[1]*x[2]*x[3]^2)
	Player2=MathProgNLPModel(P2)
	sol = x->(t=1-x[1]; (x[2]==1-t) && (x[3]==1+2*t))

 return GNEP(2,[2,1],[Player1, Player2], sol=sol)
end

function Tangi1()
	ul(i)=[0.1;0.1;0.1][i]
	ux(i)=[Inf;Inf;Inf][i]

	P1=JuMP.Model()
	JuMP.@variable(P1,x[i=1:3],upperbound=ux(i),lowerbound=ul(i), start=1.0)
	JuMP.@constraint(P1,x[1]+x[2]+x[3]<=3)
	JuMP.@NLobjective(P1,Min,0.5*(x[1]-1)^2+0.5*(x[2]-1)^2)
	Player1=MathProgNLPModel(P1)

	P2=JuMP.Model()
	JuMP.@variable(P2,x[i=1:3],upperbound=ux(i),lowerbound=ul(i), start=1.0)
	JuMP.@NLconstraint(P2,-x[1]^2*x[3]^2<=-0.5)
	JuMP.@NLconstraint(P2,-x[2]^2*x[3]^2<=-0.5)
	JuMP.@NLobjective(P2,Min,0.5*x[1]*x[2]*x[3]^2)
	Player2=MathProgNLPModel(P2)

	sol= x->(t=x[1]; (x[1]==5 && x[2]==9) || (x[2]==15-t && t <= 10 && t >= 9))
 return GNEP(2,[2,1],[Player1, Player2], sol=sol)
end

function Jason2cvx()
        ul(i)=[log(0.1);log(0.1);log(0.1)][i]
        ux(i)=[Inf;Inf;Inf][i]

	P1=JuMP.Model()
	JuMP.@variable(P1,y[i=1:3],upperbound=ux(i),lowerbound=ul(i),start=1.0)
	#JuMP.@NLconstraint(P1,exp(y[1]+log(1/3))+exp(y[2]+log(1/3))+exp(y[3]+log(1/3))<=1)
	#JuMP.@NLconstraint(P1,log(exp(y[1]+log(1/3))+exp(y[2]+log(1/3))+exp(y[3]+log(1/3)))<=0)
	JuMP.@NLconstraint(P1,log(exp(y[1])+exp(y[2])+exp(y[3]))<=log(3))
	#JuMP.@constraint(P1,y[1]>=log(0.1))
	#JuMP.@constraint(P1,y[2]>=log(0.1))
	#JuMP.@constraint(P1,y[3]>=log(0.1))
	JuMP.@NLobjective(P1,Min,0.5*(y[1])^2+0.5*(y[2])^2)
	Player1=MathProgNLPModel(P1)

	P2=JuMP.Model()
	JuMP.@variable(P2,y[i=1:3],upperbound=ux(i),lowerbound=ul(i),start=1.0)
	JuMP.@constraint(P2,2*y[1]+2*y[3]>=log(0.5))
	JuMP.@constraint(P2,2*y[2]+2*y[3]>=log(0.5))
	#JuMP.@constraint(P2,y[1]>=log(0.1))
	#JuMP.@constraint(P2,y[2]>=log(0.1))
	#JuMP.@constraint(P2,y[3]>=log(0.1))
	JuMP.@NLobjective(P2,Min,y[1]+y[2]+2*y[3])
	Player2=MathProgNLPModel(P2)
	sol = x->((x[1]==x[2] && x[1]==(0.5 + 1/sqrt(2))) || x == [1.,1.,sqrt(2)/2])
	#norm(xk-log.([1/2+1/sqrt(2),1/2+1/sqrt(2),1/(sqrt(2)*(1/2+1/sqrt(2)))]),Inf)
	#norm(xk-log.([1.,1.,sqrt(2)/2]),Inf)
 return GNEP(2,[2,1],[Player1, Player2], sol=sol)
end

function Jason3()
        ul(i)=[0;0][i]
        ux(i)=[Inf;Inf][i]

 P1=JuMP.Model()
 JuMP.@variable(P1,x[i=1:2],upperbound=ux(i),lowerbound=ul(i),start=1.0)
 JuMP.@constraint(P1,x[1]^2+x[2]^2-1<=0)
 #JuMP.@constraint(P1,x[1]>=0)
 #JuMP.@constraint(P1,x[2]>=0)
 #JuMP.@NLobjective(P1,Min,x[1]^2-x[1]*x[2]-x[2])
 JuMP.@NLobjective(P1,Min,x[1]^2-x[1]*x[2]-x[1]+x[2]^2)
 Player1=MathProgNLPModel(P1)

 P2=JuMP.Model()
 JuMP.@variable(P2,x[i=1:2],upperbound=ux(i),lowerbound=ul(i),start=1.0)
 JuMP.@constraint(P2,x[1]^2+x[2]^2-1<=0)
 #JuMP.@constraint(P2,x[1]>=0)
 #JuMP.@constraint(P2,x[2]>=0)
 #JuMP.@NLobjective(P2,Min,x[2]^2-x[1]*x[2]/2-2*x[2])
 JuMP.@NLobjective(P2,Min,x[2]^2-x[1]*x[2]/2-2*x[2]+x[1]^2)
 Player2=MathProgNLPModel(P2)
 sol=x->(t=x[1]; x[2]==sqrt(1-t^2) && t<=4/5 && t>=0)

 return GNEP(2,[1,1],[Player1, Player2], sol=sol)
end

function TangiCEx()
	P1=JuMP.Model()
	JuMP.@variable(P1,x[1:4],start=1.0)
	JuMP.@constraint(P1,x[1]+x[3]>=0)
	JuMP.@constraint(P1,x[2]+x[4]>=0)
	JuMP.@NLobjective(P1,Min,x[1]+2*x[2])
	Player1=MathProgNLPModel(P1)

	P2=JuMP.Model()
	JuMP.@variable(P2,x[1:4],start=1.0)
	JuMP.@constraint(P2,x[1]+x[3]>=0)
	JuMP.@constraint(P2,x[2]+x[4]>=0)
	JuMP.@NLobjective(P2,Min,2*x[3]+4*x[4])
	Player2=MathProgNLPModel(P2)
	sol = x->true
	#0 is a solution / may be unbounded!!

 return GNEP(2,[2,2],[Player1, Player2], sol=sol)
end

function GNEP1()
        ul(i)=[0;0][i]
        ux(i)=[Inf;Inf][i]

	P1=JuMP.Model()
	JuMP.@variable(P1,x[i=1:2],upperbound=ux(i),lowerbound=ul(i),start=1.0)
	JuMP.@constraint(P1,x[1]^2+x[2]^2-1<=0)
	JuMP.@constraint(P1,x[1]+x[2]-1<=0)
	JuMP.@NLobjective(P1,Min,0.5*(x[1]-1)^2)
	Player1=MathProgNLPModel(P1)

	P2=JuMP.Model()
	JuMP.@variable(P2,x[i=1:2],upperbound=ux(i),lowerbound=ul(i),start=1.0)
	JuMP.@constraint(P2,x[1]+x[2]-1<=0)
	JuMP.@NLobjective(P2,Min,0.5*(x[2]-0.5)^2)
	Player2=MathProgNLPModel(P2)
	sol = x->(t=x[1]; x[2]==1-t && 0.5<=t && t<=1)
	#VE: (3/4,1/4)

 return GNEP(2,[1,1],[Player1, Player2], sol=sol)
end

function AbsVal()
	P1=JuMP.Model()
	JuMP.@variable(P1,x[1:2],start=1.0)
	JuMP.@constraint(P1,x[1]>=-1)
	JuMP.@constraint(P1,x[1]<=1)
	JuMP.@NLobjective(P1,Min,0.)
	Player1=MathProgNLPModel(P1)

	P2=JuMP.Model()
	JuMP.@variable(P2,x[1:2],start=1.0)
	JuMP.@constraint(P2,x[1]>=-1)
	JuMP.@constraint(P2,x[1]<=1)
	JuMP.@constraint(P2,x[2]>=-x[1])
	JuMP.@constraint(P2,x[2]<=x[1])
	JuMP.@NLobjective(P2,Min,-x[2])
	Player2=MathProgNLPModel(P2)
	sol = x->(t=x[1]; x[2]==abs(t) && -1<=t && t<=1)
	#VE: (3/4,1/4)

 return GNEP(2,[1,1],[Player1, Player2], sol=sol)
end

function A8()
        ul(i)=[0;0;0][i]
        ux(i)=[1;1;2][i]

	P1=JuMP.Model()
	JuMP.@variable(P1,x[i=1:3],upperbound=ux(i),lowerbound=ul(i), start=1.0)
	JuMP.@constraint(P1,x[1]+x[2]-1<=0)
	JuMP.@constraint(P1,x[1]+x[2]-x[3]>=0)
	#JuMP.@constraint(P1,x[1]>=0)
	#JuMP.@constraint(P1,x[2]>=0)
	#JuMP.@constraint(P1,x[3]>=0)
	#JuMP.@constraint(P1,x[3]<=2)
	JuMP.@NLobjective(P1,Min,-x[1])
	Player1=MathProgNLPModel(P1)

	P2=JuMP.Model()
	JuMP.@variable(P2,x[i=1:3],upperbound=ux(i),lowerbound=ul(i), start=1.0)
	JuMP.@constraint(P2,x[1]+x[2]-1<=0)
	JuMP.@constraint(P2,x[1]+x[2]-x[3]>=0)
	#JuMP.@constraint(P1,x[1]>=0)
	#JuMP.@constraint(P1,x[2]>=0)
	#JuMP.@constraint(P1,x[3]>=0)
	#JuMP.@constraint(P1,x[3]<=2)
	JuMP.@NLobjective(P2,Min,0.5*(x[2]-0.5)^2)
	Player2=MathProgNLPModel(P2)

	P3=JuMP.Model()
	JuMP.@variable(P3,x[i=1:3],upperbound=ux(i),lowerbound=ul(i), start=1.0)
	#JuMP.@constraint(P1,x[1]>=0)
	#JuMP.@constraint(P1,x[2]>=0)
	#JuMP.@constraint(P1,x[3]>=0)
	#JuMP.@constraint(P1,x[3]<=2)
	JuMP.@NLobjective(P3,Min,0.5*(x[3]-1.5*x[1])^2)
	Player3=MathProgNLPModel(P3)
	sol = x->(t=x[1]; x[2]==1-t && x[3]==1.5*t && 0.5<=t && t<=2/3)

 return GNEP(3,[1,1,1],[Player1, Player2, Player3], sol=sol)
end

function A8slack()
        ul(i)=[0;0;0;0;0;0;0][i]
        ux(i)=[1;1;2;Inf;Inf;Inf;Inf][i]

	P1=JuMP.Model()
	JuMP.@variable(P1,x[i=1:7],upperbound=ux(i),lowerbound=ul(i), start=1.0)
	JuMP.@constraint(P1,x[1]+x[4]-1==-x[2])
	JuMP.@constraint(P1,x[1]+x[4]-x[7]==x[3])
	JuMP.@NLobjective(P1,Min,-x[1])
	Player1=MathProgNLPModel(P1)

	P2=JuMP.Model()
	JuMP.@variable(P2,x[i=1:7],upperbound=ux(i),lowerbound=ul(i), start=1.0)
	JuMP.@constraint(P2,x[1]+x[4]-1==-x[5])
	JuMP.@constraint(P2,x[1]+x[4]-x[7]==x[6])
	JuMP.@NLobjective(P2,Min,0.5*(x[2]-0.5)^2)
	Player2=MathProgNLPModel(P2)

	P3=JuMP.Model()
	JuMP.@variable(P3,x[i=1:7],upperbound=ux(i),lowerbound=ul(i), start=1.0)
	JuMP.@NLobjective(P3,Min,0.5*(x[7]-1.5*x[1])^2)
	Player3=MathProgNLPModel(P3)
	sol = x->(t=x[1]; x[2]==1-t && x[3]==1.5*t && 0.5<=t && t<=2/3)

 return GNEP(3,[3,3,1],[Player1, Player2, Player3], sol=sol)
end

function Riverbasin()
	a1=0.01; a2=0.05; a3=0.01; b=0.01; e1=2.9; e2=2.88; e3=2.85
        ul(i)=[0;0;0][i]
        ux(i)=[Inf;Inf;Inf][i]

	P1=JuMP.Model()
	JuMP.@variable(P1,x[i=1:3],upperbound=ux(i),lowerbound=ul(i), start=1.0)
	JuMP.@constraint(P1,3.25*x[1]+1.25*x[2]+4.125*x[3]<=100)
	JuMP.@constraint(P1,2.29115*x[1]+1.5625*x[2]+2.8125*x[3]<=100)
	JuMP.@NLobjective(P1,Min,(a1*x[1]+b*(x[1]+x[2]+x[3])-e1)*x[1])
	Player1=MathProgNLPModel(P1)

	P2=JuMP.Model()
	JuMP.@variable(P2,x[i=1:3],upperbound=ux(i),lowerbound=ul(i), start=1.0)
	JuMP.@constraint(P2,3.25*x[1]+1.25*x[2]+4.125*x[3]<=100)
	JuMP.@constraint(P2,2.29115*x[1]+1.5625*x[2]+2.8125*x[3]<=100)
	JuMP.@NLobjective(P2,Min,(a2*x[2]+b*(x[1]+x[2]+x[3])-e2)*x[2])
	Player2=MathProgNLPModel(P2)

	P3=JuMP.Model()
	JuMP.@variable(P3,x[i=1:3],upperbound=ux(i),lowerbound=ul(i), start=1.0)
	JuMP.@constraint(P3,3.25*x[1]+1.25*x[2]+4.125*x[3]<=100)
	JuMP.@constraint(P3,2.29115*x[1]+1.5625*x[2]+2.8125*x[3]<=100)
	JuMP.@NLobjective(P3,Min,(a3*x[3]+b*(x[1]+x[2]+x[3])-e3)*x[3])
	Player3=MathProgNLPModel(P3)
	sol = x->(t=x[1]; x[2]==1-t && x[3]==1.5*t && 0.5<=t && t<=2/3)

 return GNEP(3,[1,1,1],[Player1, Player2, Player3], sol=sol)
end

function Harker()
	ul(i)=[0;0][i]
	ux(i)=[10;10][i]

	P1=JuMP.Model()
	JuMP.@variable(P1,x[i=1:2],upperbound=ux(i),lowerbound=ul(i), start=1.0)
	JuMP.@constraint(P1,x[1]+x[2]<=15)
	JuMP.@NLobjective(P1,Min,x[1]^2+8/3*x[1]*x[2]-34*x[1])
	Player1=MathProgNLPModel(P1)

	P2=JuMP.Model()
	JuMP.@variable(P2,x[i=1:2],upperbound=ux(i),lowerbound=ul(i), start=1.0)
	JuMP.@constraint(P2,x[1]+x[2]<=15)
	JuMP.@NLobjective(P2,Min,x[2]^2+5/4*x[1]*x[2]-24.25*x[2])
	Player2=MathProgNLPModel(P2)

	sol= x->(t=x[1]; (x[1]==5 && x[2]==9) || (x[2]==15-t && t <= 10 && t >= 9))
 return GNEP(2,[1,1],[Player1, Player2], sol=sol)
end

function ContrHarker()
	ul(i)=[0;0][i]
	ux(i)=[10;10][i]

        theta = 0.9

	P1=JuMP.Model()
	JuMP.@variable(P1,x[i=1:2],upperbound=ux(i),lowerbound=ul(i), start=1.0)
	JuMP.@constraint(P1,x[1]+0.9*x[2]<=14.4)
	JuMP.@NLobjective(P1,Min,x[1]^2+8/3*x[1]*x[2]-34*x[1])
	Player1=MathProgNLPModel(P1)

	P2=JuMP.Model()
	JuMP.@variable(P2,x[i=1:2],upperbound=ux(i),lowerbound=ul(i), start=1.0)
	JuMP.@constraint(P2,0.9*x[1]+x[2]<=14.1)
	JuMP.@NLobjective(P2,Min,x[2]^2+5/4*x[1]*x[2]-24.25*x[2])
	Player2=MathProgNLPModel(P2)

	sol= x->(t=x[1]; (x[1]==5 && x[2]==9) || (x[2]==15-t && t <= 10 && t >= 9))
 return GNEP(2,[1,1],[Player1, Player2], sol=sol)
end

function Nabutani1()
	ul(i)=[0;0;0][i]
	ux(i)=[Inf;Inf;Inf][i]

	P1=JuMP.Model()
	JuMP.@variable(P1,x[i=1:3],upperbound=ux(i),lowerbound=ul(i), start=1.0)
	JuMP.@constraint(P1,x[1]+2*x[2]-x[3]<=14)
	JuMP.@constraint(P1,3*x[1]+2*x[2]+x[3]<=30)
	JuMP.@NLobjective(P1,Min,x[1]^2+x[1]*x[2]+x[2]^2+(x[1]+x[2])*x[3]-25*x[1]-38*x[2])
	Player1=MathProgNLPModel(P1)

	P2=JuMP.Model()
	JuMP.@variable(P2,x[i=1:3],upperbound=ux(i),lowerbound=ul(i), start=1.0)
	JuMP.@constraint(P2,x[1]+2*x[2]-x[3]<=14)
	JuMP.@constraint(P2,3*x[1]+2*x[2]+x[3]<=30)
	JuMP.@NLobjective(P2,Min,x[3]^2+(x[1]+x[2])*x[3]-25*x[3])
	Player2=MathProgNLPModel(P2)

	sol= x->(t=x[1]; x[2] == 11-t && x[3]==8-t && t >= 0 && t<= 2)
 return GNEP(2,[2,1],[Player1, Player2], sol=sol)
end

function Ex1()
	ul(i)=[0;0][i]
	ux(i)=[Inf;Inf][i]

	P1=JuMP.Model()
	JuMP.@variable(P1,x[i=1:2],upperbound=ux(i),lowerbound=ul(i), start=1.0)
	JuMP.@constraint(P1,x[1]+x[2]<=1)
	JuMP.@NLobjective(P1,Min,x[1]^2-x[1]*x[2]-x[1])
	Player1=MathProgNLPModel(P1)

	P2=JuMP.Model()
	JuMP.@variable(P2,x[i=1:2],upperbound=ux(i),lowerbound=ul(i), start=1.0)
	JuMP.@constraint(P2,x[1]+x[2]<=1)
	JuMP.@NLobjective(P2,Min,x[2]^2-1/2*x[1]*x[2]-2*x[2])
	Player2=MathProgNLPModel(P2)

	sol= x->(t=x[1]; x[2]==1-t && t <= 2/3 && t>=0)
 return GNEP(2,[1,1],[Player1, Player2], sol=sol)
end

function Ex1bis()
	ul(i)=[0;0][i]
	ux(i)=[Inf;Inf][i]

	P1=JuMP.Model()
	JuMP.@variable(P1,x[i=1:2],upperbound=ux(i),lowerbound=ul(i), start=1.0)
	JuMP.@constraint(P1,x[1]^2+x[2]^2<=1)
	JuMP.@NLobjective(P1,Min,x[1]^2-x[1]*x[2]-x[1])
	Player1=MathProgNLPModel(P1)

	P2=JuMP.Model()
	JuMP.@variable(P2,x[i=1:2],upperbound=ux(i),lowerbound=ul(i), start=1.0)
	JuMP.@constraint(P2,x[1]^2+x[2]^2<=1)
	JuMP.@NLobjective(P2,Min,x[2]^2-1/2*x[1]*x[2]-2*x[2])
	Player2=MathProgNLPModel(P2)

	sol= x->(t=x[1]; x[2]==sqrt(1-t^2) && t <= 4/5 && t>=0)
 return GNEP(2,[1,1],[Player1, Player2], sol=sol)
end

function ContreExDS()
	ul(i)=[-Inf;-Inf][i]
	ux(i)=[Inf;Inf][i]
	L=2

	P1=JuMP.Model()
	JuMP.@variable(P1,x[i=1:2],upperbound=ux(i),lowerbound=ul(i), start=1.0)
	JuMP.@constraint(P1,x[1]<=1)
	JuMP.@constraint(P1,-x[1]+L*x[2]<=0)
	JuMP.@NLobjective(P1,Min,x[1])
	Player1=MathProgNLPModel(P1)

	P2=JuMP.Model()
	JuMP.@variable(P2,x[i=1:2],upperbound=ux(i),lowerbound=ul(i), start=1.0)
	JuMP.@constraint(P2,x[1]-x[2]==0)
	JuMP.@NLobjective(P2,Min,(x[1]-x[2])^2)
	Player2=MathProgNLPModel(P2)

	sol= x->(t=x[1]; x[2]==t && t==L*x[2]) #no solutions if L>1, (1,1) if L=1, any solution if L<1
 return GNEP(2,[1,1],[Player1, Player2], sol=sol)
end

function Ex4()
	ul(i)=[0.;-Inf;0.;-Inf][i]
	ux(i)=[Inf;Inf;Inf;Inf][i]

	P1=JuMP.Model()
	JuMP.@variable(P1,x[i=1:4],upperbound=ux(i),lowerbound=ul(i), start=1.0)
	JuMP.@constraint(P1,x[1]+x[3]<=1)
	JuMP.@constraint(P1,x[1]-2*x[3]-x[2]<=0)
	JuMP.@constraint(P1,-x[1]-x[3]-x[2]<=0)
	JuMP.@NLobjective(P1,Min,x[2])
	Player1=MathProgNLPModel(P1)

	P2=JuMP.Model()
	JuMP.@variable(P2,x[i=1:4],upperbound=ux(i),lowerbound=ul(i), start=1.0)
	JuMP.@constraint(P2,x[1]+x[3]<=1)
	JuMP.@constraint(P2,-x[1]+x[3]-x[4]<=0)
	JuMP.@constraint(P2,x[1]-x[3]-x[4]<=0)
	JuMP.@NLobjective(P2,Min,x[4])
	Player2=MathProgNLPModel(P2)

	sol= x->(true)#((x[3]=1-x[1]) || ((x[3]==x[1]) && (x[1]<=0.5)))
 return GNEP(2,[2,2],[Player1, Player2], sol=sol)
end

function Ex5()
        ul(i)=[-Inf;-Inf;0][i]
        ux(i)=[Inf;Inf;Inf][i]

	P1=JuMP.Model()
	JuMP.@variable(P1,x[i=1:3],upperbound=ux(i),lowerbound=ul(i), start=1.0)
	JuMP.@constraint(P1,x[1]+x[2]+x[3]<=1)
	JuMP.@NLobjective(P1,Min,0.5*(x[1]-1)^2-x[1]*x[2])
	Player1=MathProgNLPModel(P1)

	P2=JuMP.Model()
	JuMP.@variable(P2,x[i=1:3],upperbound=ux(i),lowerbound=ul(i), start=1.0)
	JuMP.@constraint(P2,x[1]+x[2]+x[3]<=1)
	JuMP.@NLobjective(P2,Min,0.5*(x[2]-1)^2+x[1]*x[2])
	Player2=MathProgNLPModel(P2)

	P3=JuMP.Model()
	JuMP.@variable(P3,x[i=1:3],upperbound=ux(i),lowerbound=ul(i), start=1.0)
	JuMP.@constraint(P3,-x[1]-x[2]+x[3]<=0)
	JuMP.@NLobjective(P3,Min,0.5*(x[3]-1)^2)
	Player3=MathProgNLPModel(P3)

	sol = x->(true)

 return GNEP(3,[1,1,1],[Player1, Player2, Player3], sol=sol)
end

function Ex6()
	ul(i)=[-Inf;-Inf][i]
	ux(i)=[Inf;Inf][i]

	P1=JuMP.Model()
	JuMP.@variable(P1,x[i=1:2],upperbound=ux(i),lowerbound=ul(i), start=1.0)
	JuMP.@constraint(P1,1/6*x[1]^2+x[2]-5/2<=0)
	JuMP.@NLobjective(P1,Min,0.5*x[1]^2+32/5*x[1])
	Player1=MathProgNLPModel(P1)

	P2=JuMP.Model()
	JuMP.@variable(P2,x[i=1:2],upperbound=ux(i),lowerbound=ul(i), start=1.0)
	JuMP.@NLobjective(P2,Min,0.5*x[2]^2+x[1]*x[2]-4/5*x[2])
	Player2=MathProgNLPModel(P2)

	sol= x->(true)
 return GNEP(2,[1,1],[Player1, Player2], sol=sol)
end

function Ex7()
	ul(i)=[-Inf;-Inf][i]
	ux(i)=[Inf;Inf][i]

	P1=JuMP.Model()
	JuMP.@variable(P1,x[i=1:2],upperbound=ux(i),lowerbound=ul(i), start=1.0)
	JuMP.@constraint(P1,x[1]^2+x[2]<=1)
	JuMP.@NLobjective(P1,Min,x[1])
	Player1=MathProgNLPModel(P1)

	P2=JuMP.Model()
	JuMP.@variable(P2,x[i=1:2],upperbound=ux(i),lowerbound=ul(i), start=1.0)
	JuMP.@NLobjective(P2,Min,0.5*x[2]^2)
	Player2=MathProgNLPModel(P2)

	sol= x->((x[1]==-1.) && (x[2] == 0.))
 return GNEP(2,[1,1],[Player1, Player2], sol=sol)
end

function Ex8()
	ul(i)=[-Inf;-Inf][i]
	ux(i)=[Inf;Inf][i]

	P1=JuMP.Model()
	JuMP.@variable(P1,x[i=1:2],upperbound=ux(i),lowerbound=ul(i), start=1.0)
	JuMP.@constraint(P1,x[1]^2+x[2]<=0)
	JuMP.@NLobjective(P1,Min,0.5*x[1]^2-2*x[1])
	Player1=MathProgNLPModel(P1)

	P2=JuMP.Model()
	JuMP.@variable(P2,x[i=1:2],upperbound=ux(i),lowerbound=ul(i), start=1.0)
	JuMP.@NLobjective(P2,Min,0.5*x[2]^2+(2-x[1]^2)*x[2])
	Player2=MathProgNLPModel(P2)

	sol= x->((x[1]==-1.) && (x[2] == 0.))
 return GNEP(2,[1,1],[Player1, Player2], sol=sol)
end

"""
Arrow and Debreu model of a competitive economy proposed in:
Kenneth J. Arrow and Gerard Debreu. Existence of an equilibrium for a competitive
economy. Econometrica: Journal of the Econometric Society, pages 265–290, 1954.

Also known as problem A10 in:
Francisco Facchinei and Christian Kanzow. Penalty methods for the solution of general-
ized Nash equilibrium problems (with complete test problems). Sapienza University of
Rome, 2009.

This problem is specialized in several variants: a), b), c), d) and e)
"""
function A10a()

	include("A10a.jl")
	arr_player = Array{NLPModelsJuMP.MathProgNLPModel,1}(undef,F+C+1)
        ul(i)=zeros((F+C+1)*P)[i]
	ux(i)=vcat(sqrt(21)*ones(P*F),6*ones(P*C),1.1*ones(P))[i]

 #Firms problems
 for j=1:F
	firm_j = JuMP.Model()
	JuMP.@variable(firm_j,x[i=1:(F+C+1)*P],upperbound=ux(i),lowerbound=ul(i), start=1.0)
	#JuMP.@NLobjective(firm_j,Min,x[(F+C)*P+1:n]'*x[P*(j-1)+1:P*(j-1)+P])
	JuMP.@NLobjective(firm_j,Min,-x[22]*x[P*(j-1)+1]-x[23]*x[P*(j-1)+2]-x[24]*x[P*(j-1)+3])
	JuMP.@constraint(firm_j,sum(x[P*(j-1)+1:P*(j-1)+P].^2)<=10*j)
	arr_player[j] = MathProgNLPModel(firm_j)
 end
 #Consumers problems
 for i=1:C
	cons_i = JuMP.Model()
	JuMP.@variable(cons_i,x[i=1:(F+C+1)*P],upperbound=ux(i),lowerbound=ul(i), start=1.0)
  #error max not min
	#JuMP.@NLobjective(cons_i,Min, -0.5*(x[1]*Q[i,1,1]*x[1]+x[1]*Q[i,1,2]*x[2]+x[1]*Q[i,1,1]*x[3] + x[2]*Q[i,2,1]*x[1]+x[2]*Q[i,2,2]*x[2]+x[2]*Q[i,2,1]*x[3] + x[3]*Q[i,3,1]*x[1]+x[3]*Q[i,3,2]*x[2]+x[3]*Q[i,3,1]*x[3]) + x[P*(F+i-1)+1]*b[i,1,1] + x[P*(F+i-1)+2]*b[i,2,1] +x[P*(F+i-1)+3]*b[i,3,1])
  JuMP.@NLobjective(cons_i,Min, 0.5*(x[1]*Q[i,1,1]*x[1]+x[1]*Q[i,1,2]*x[2]+x[1]*Q[i,1,1]*x[3] + x[2]*Q[i,2,1]*x[1]+x[2]*Q[i,2,2]*x[2]+x[2]*Q[i,2,1]*x[3] + x[3]*Q[i,3,1]*x[1]+x[3]*Q[i,3,2]*x[2]+x[3]*Q[i,3,1]*x[3]) - x[P*(F+i-1)+1]*b[i,1,1] - x[P*(F+i-1)+2]*b[i,2,1] - x[P*(F+i-1)+3]*b[i,3,1])
	JuMP.@constraint(cons_i,x[(F+C)*P+1:(F+C+1)*P]'*x[P*(F+i-1)+1:P*(F+i-1)+P]-x[(F+C)*P+1:(F+C+1)*P]'*xi[i,:,1]<=0)
	arr_player[F+i] = MathProgNLPModel(cons_i)
 end
 #Market player problem
	mark = JuMP.Model()
	JuMP.@variable(mark,x[i=1:(F+C+1)*P],upperbound=ux(i),lowerbound=ul(i), start=1.0)
	JuMP.@NLobjective(mark,Min,-x[22]*(x[7]+x[10]+x[13]+x[16]+x[19]-x[1]-x[4]-xi[1,1,1]-xi[2,1,1]-xi[3,1,1]-xi[4,1,1]-xi[5,1,1])-x[23]*(x[8]+x[11]+x[14]+x[17]+x[20]-x[2]-x[5]-xi[1,2,1]-xi[2,2,1]-xi[3,2,1]-xi[4,2,1]-xi[5,2,1])-x[24]*(x[9]+x[12]+x[15]+x[18]+x[21]-x[3]-x[5]-xi[1,3,1]-xi[2,3,1]-xi[3,3,1]-xi[4,3,1]-xi[5,3,1]))
	JuMP.@constraint(mark,sum(x[(F+C)*P+1:(F+C+1)*P])==1)
	arr_player[F+C+1] = MathProgNLPModel(mark)

 sol = x->(true)

 #initial point in the paper:
 x0 = zeros(24);x0[22:24].=1/3;
 return GNEPmod.GNEP(F+C+1,P*ones(Int64,F+C+1),arr_player, sol=sol, x0=x0)
end

"""
Inspired by an internet switching model introduced by

Alex Kesselman, Stefano Leonardi, and Vincenzo Bonifaci. Game-theoretic analysis of
internet switching with selfish users. In Xiaotie Deng and Yinyu Ye, editors, Internet and
Network Economics, pages 236–245, Berlin, Heidelberg, 2005. Springer Berlin Heidelberg.

and also analyzed as problem A14 in
Francisco Facchinei and Christian Kanzow. Penalty methods for the solution of general-
ized Nash equilibrium problems (with complete test problems). Sapienza University of
Rome, 2009.

This is my* specialization.
"""
function TwoGroupsSwitch(B1,B2)

	N1=2;N2=3;
	arr_player = Array{NLPModelsJuMP.MathProgNLPModel,1}(undef,N1+N2)
	ul(i)=0.01*ones(Int64,N1+N2)[i]
	ux(i)=Inf*ones(Int64,N1+N2)[i]

	for k=1:N1
	P1=JuMP.Model()
	JuMP.@variable(P1,x[i=1:N1+N2],upperbound=ux(i),lowerbound=ul(i), start=1.0)
	JuMP.@constraint(P1,sum(x[1:N1+N2])<=B1)
	JuMP.@NLobjective(P1,Min,-x[k]/(x[1]+x[2]+x[3]+x[4]+x[5]))
	arr_player[k]=MathProgNLPModel(P1)
	end

	for j=N1+1:N1+N2
	P2=JuMP.Model()
	JuMP.@variable(P2,x[i=1:N1+N2],upperbound=ux(i),lowerbound=ul(i), start=1.0)
	JuMP.@constraint(P2,sum(x[1:N1+N2].^2)<=B2)
	JuMP.@NLobjective(P2,Min,-x[j]/(x[1]+x[2]+x[3]+x[4]+x[5]))
	arr_player[j]=MathProgNLPModel(P2)
	end

	sol= x->((x[1]==-1.) && (x[2] == 0.))
 return GNEP(N1+N2, ones(Int64,N1+N2), arr_player, sol=sol)
end

using LinearAlgebra
"""
Environmental protocol control proposed in
Breton, M., Zaccour, G., & Zahaf, M. (2006).
A game-theoretic formulation of joint implementation of environmental projects.
European Journal of Operational Research, 168(1), 221-239.
"""
function EPC()

    include("EPC.jl")

    arr_player = Array{NLPModelsJuMP.MathProgNLPModel,1}(undef,N)

    ub = 50*ones(Float64,N*(N+1))
    ub[union(1,5,9)] = p

	ul(i)=zeros(Float64,N*(N+1))[i]
	ux(i)=ub[i]

	for k=1:N
	P1=JuMP.Model()
        nu = 1+(k-1)*(N+1) #index of e_k
        j  = mod(k,N)+1
        l  = mod(k+1,N)+1
	JuMP.@variable(P1,x[i=1:N*(N+1)],upperbound=ux(i),lowerbound=ul(i), start=1.0)
	JuMP.@constraint(P1,dot(x[nu:(N+nu)],vcat(1.0,-gam))<=E[k])
	JuMP.@constraint(P1,-dot(x[union(1,2,6,10)],vcat(1.0,-gam[1]*ones(N)))<=0.0)
	JuMP.@constraint(P1,-dot(x[union(5,3,7,11)],vcat(1.0,-gam[2]*ones(N)))<=0.0)
	JuMP.@constraint(P1,-dot(x[union(9,4,8,12)],vcat(1.0,-gam[3]*ones(N)))<=0.0)
    JuMP.@NLobjective(P1,Min,-x[nu]*(p[k]-0.5*x[nu])
                             +0.5*(x[nu+k]^2+x[nu+j]^2+x[nu+l]^2
                                   +2*x[nu+j]*x[1+(j-1)*(N+1)+j]
                                   +2*x[nu+l]*x[1+(l-1)*(N+1)+l])
                             +lambda[1]*(x[1]-gam[1]*(x[2]+x[6]+x[10]))
                             +lambda[2]*(x[5]-gam[2]*(x[3]+x[7]+x[11]))
                             +lambda[3]*(x[9]-gam[3]*(x[4]+x[8]+x[12])))
    arr_player[k]=MathProgNLPModel(P1)
	end

	sol= x->((x[1]==-1.) && (x[2] == 0.))
 return GNEP(N, (N+1)*ones(Int64,N), arr_player, sol=sol)
end

"""
Environmental protocol control proposed in
Breton, M., Zaccour, G., & Zahaf, M. (2006).
A game-theoretic formulation of joint implementation of environmental projects.
European Journal of Operational Research, 168(1), 221-239.

Variant of EPC() studied in
Cojocaru, M. G., Wild, E., & Small, A. (2018).
On describing the solution sets of generalized Nash games with shared constraints.
Optimization and Engineering, 19(4), 845-870.
"""
function EPCshared()

    include("EPC.jl")

    ub = 50*ones(Float64,N*(N+1))
    ub[union(1,5,9)] = p

    arr_player = Array{NLPModelsJuMP.MathProgNLPModel,1}(undef,N)

	ul(i)=zeros(Float64,N*(N+1))[i]
	ux(i)=ub[i]

	for k=1:N
	P1=JuMP.Model()
        nu = 1+(k-1)*(N+1) #index of e_k
        j  = mod(k,N)+1
        l  = mod(k+1,N)+1
	JuMP.@variable(P1,x[i=1:N*(N+1)],upperbound=ux(i),lowerbound=ul(i), start=1.0)
	JuMP.@constraint(P1,dot(x[nu:(N+nu)],vcat(1.0,-gam))<=E[k])
	JuMP.@constraint(P1,-dot(x[union(1,2,6,10)],vcat(1.0,-gam[1]*ones(N)))
                        -dot(x[union(5,3,7,11)],vcat(1.0,-gam[2]*ones(N)))
                        -dot(x[union(9,4,8,12)],vcat(1.0,-gam[3]*ones(N)))<=0.0)
    JuMP.@NLobjective(P1,Min,-x[nu]*(p[k]-0.5*x[nu])
                             +0.5*(x[nu+k]^2+x[nu+j]^2+x[nu+l]^2
                                   +2*x[nu+j]*x[1+(j-1)*(N+1)+j]
                                   +2*x[nu+l]*x[1+(l-1)*(N+1)+l])
                             +lambda[1]*(x[1]-gam[1]*(x[2]+x[6]+x[10]))
                             +lambda[2]*(x[5]-gam[2]*(x[3]+x[7]+x[11]))
                             +lambda[3]*(x[9]-gam[3]*(x[4]+x[8]+x[12])))
     arr_player[k]=MathProgNLPModel(P1)
	end

	sol= x->((x[1]==-1.) && (x[2] == 0.))
 return GNEP(N, (N+1)*ones(Int64,N), arr_player, sol=sol)
end

"""
Environmental protocol control proposed in
Breton, M., Zaccour, G., & Zahaf, M. (2006).
A game-theoretic formulation of joint implementation of environmental projects.
European Journal of Operational Research, 168(1), 221-239.

Variant of EPC() to find the emission of each country
"""
function regularEPC() #it is a Nash game

    include("EPC.jl")
    ub = 100*ones(Float64,6)

    arr_player = Array{NLPModelsJuMP.MathProgNLPModel,1}(undef,N)

	ul(i)=zeros(Int64,2*N)[i]
	ux(i)=ub[i]

	for k=1:N
	P1=JuMP.Model()
    nu = 1+(k-1)*2 #index of e_k
	JuMP.@variable(P1,x[i=1:2*N],upperbound=ux(i),lowerbound=ul(i), start=1.0)
	JuMP.@constraint(P1,x[nu]-gam[k]*x[nu+1]<=E[k])
    JuMP.@NLobjective(P1,Min,-x[nu]*(p[k]-0.5*x[nu])+0.5*(x[nu+1]^2)
	                         +lambda[1]*(x[1]-gam[1]*(x[2]))
							 +lambda[2]*(x[3]-gam[2]*(x[4]))
							 +lambda[3]*(x[5]-gam[3]*(x[6])))
	arr_player[k]=MathProgNLPModel(P1)
	end

	sol= x->((x[1]==-1.) && (x[2] == 0.))
 return GNEP(N, 2*ones(Int64,N), arr_player, sol=sol)
end

"""
Environmental protocol control proposed in
Breton, M., Zaccour, G., & Zahaf, M. (2006).
A game-theoretic formulation of joint implementation of environmental projects.
European Journal of Operational Research, 168(1), 221-239.

Variant of EPC() studied in *** new ***.
"""
function EPCnonshared()

    include("EPC.jl")

    ub = 50*ones(Float64,N*(N+1))
    ub[union(1,5,9)] = p

    arr_player = Array{NLPModelsJuMP.MathProgNLPModel,1}(undef,N)

	ul(i)=zeros(Float64,N*(N+1))[i]
	ux(i)=ub[i]

	for k=1:N
	P1=JuMP.Model()
        nu = 1+(k-1)*(N+1) #index of e_k
        j  = mod(k,N)+1
        l  = mod(k+1,N)+1
	JuMP.@variable(P1,x[i=1:N*(N+1)],upperbound=ux(i),lowerbound=ul(i), start=1.0)
	JuMP.@constraint(P1,dot(x[nu:(N+nu)],vcat(1.0,-gam))<=E[k])
	JuMP.@constraint(P1,-dot(x[union(1,2,6,10)],vcat(1.0,-gam[1]*ones(N)))
                        -dot(x[union(5,3,7,11)],vcat(1.0,-gam[2]*ones(N)))
                        -dot(x[union(9,4,8,12)],vcat(1.0,-gam[3]*ones(N)))<=0.0)
    if (k==1) || (k==2)
        JuMP.@NLconstraint(P1,-(x[2]+x[6])/(x[2]+x[3]+x[4]+x[6])
    						  -(x[3]+x[7])/(x[6]+x[7]+x[8]+x[3])<=-gap)
    end
    JuMP.@NLobjective(P1,Min,-x[nu]*(p[k]-0.5*x[nu])
                             +0.5*(x[nu+k]^2+x[nu+j]^2+x[nu+l]^2
                                   +2*x[nu+j]*x[1+(j-1)*(N+1)+j]
                                   +2*x[nu+l]*x[1+(l-1)*(N+1)+l])
                             +lambda[1]*(x[1]-gam[1]*(x[2]+x[6]+x[10]))
                             +lambda[2]*(x[5]-gam[2]*(x[3]+x[7]+x[11]))
                             +lambda[3]*(x[9]-gam[3]*(x[4]+x[8]+x[12])))
     arr_player[k]=MathProgNLPModel(P1)
	end

	sol= x->((x[1]==-1.) && (x[2] == 0.))
 return GNEP(N, (N+1)*ones(Int64,N), arr_player, sol=sol)
end

"""
Environmental protocol control proposed in
Breton, M., Zaccour, G., & Zahaf, M. (2006).
A game-theoretic formulation of joint implementation of environmental projects.
European Journal of Operational Research, 168(1), 221-239.

Variant of EPC() studied in
Cojocaru, M. G., Wild, E., & Small, A. (2018).
On describing the solution sets of generalized Nash games with shared constraints.
Optimization and Engineering, 19(4), 845-870.

N=5
"""
function EPC5shared()

    include("EPC5.jl")

    ub = 50*ones(Float64,N*(N+1))
    ub[union(1,7,13,19,25)] = p

    arr_player = Array{NLPModelsJuMP.MathProgNLPModel,1}(undef,N)

	ul(i)=zeros(Float64,N*(N+1))[i]
	ux(i)=ub[i]

	for k=1:N
	P1=JuMP.Model()
        nu = 1+(k-1)*(N+1) #index of e_k
        j  = mod(k,N)+1
        l  = mod(k+1,N)+1
        o  = mod(k+2,N)+1
        m  = mod(k+3,N)+1
	JuMP.@variable(P1,x[i=1:N*(N+1)],upperbound=ux(i),lowerbound=ul(i), start=1.0)
	JuMP.@constraint(P1,dot(x[nu:(N+nu)],vcat(1.0,-gam))<=E[k])
	JuMP.@constraint(P1,-dot(x[union(1, 2, 8,14,20,26)],vcat(1.0,-gam[1]*ones(N)))
                        -dot(x[union(7, 3, 9,15,21,27)],vcat(1.0,-gam[2]*ones(N)))
                        -dot(x[union(13,4,10,16,22,28)],vcat(1.0,-gam[3]*ones(N)))
                        -dot(x[union(19,5,11,17,23,29)],vcat(1.0,-gam[4]*ones(N)))
                        -dot(x[union(25,6,12,18,24,30)],vcat(1.0,-gam[5]*ones(N)))<=0.0)
    JuMP.@NLobjective(P1,Min,-x[nu]*(p[k]-0.5*x[nu])
                             +0.5*(x[nu+k]^2+x[nu+j]^2+x[nu+l]^2+x[nu+o]^2+x[nu+m]^2
                                   +2*x[nu+j]*x[1+(j-1)*(N+1)+j]
                                   +2*x[nu+l]*x[1+(l-1)*(N+1)+l]
                                   +2*x[nu+o]*x[1+(o-1)*(N+1)+o]
                                   +2*x[nu+m]*x[1+(m-1)*(N+1)+m])
                             +lambda[1]*(x[1]-gam[1]*(x[2]+x[8]+x[14]+x[20]+x[26]))
                             +lambda[2]*(x[7]-gam[2]*(x[3]+x[9]+x[15]+x[21]+x[27]))
                             +lambda[3]*(x[13]-gam[3]*(x[4]+x[10]+x[16]+x[22]+x[28]))
                             +lambda[4]*(x[19]-gam[4]*(x[5]+x[11]+x[17]+x[23]+x[29]))
                             +lambda[5]*(x[25]-gam[5]*(x[6]+x[12]+x[18]+x[24]+x[30])))
     arr_player[k]=MathProgNLPModel(P1)
	end

	sol= x->((x[1]==-1.) && (x[2] == 0.))
 return GNEP(N, (N+1)*ones(Int64,N), arr_player, sol=sol)
end

"""
Environmental protocol control proposed in
Breton, M., Zaccour, G., & Zahaf, M. (2006).
A game-theoretic formulation of joint implementation of environmental projects.
European Journal of Operational Research, 168(1), 221-239.

Variant of EPC() to find the emission of each country
"""
function regularEPC5() #it is a Nash game

    include("EPC5.jl")
    ub = 100*ones(Float64,2*N)

    arr_player = Array{NLPModelsJuMP.MathProgNLPModel,1}(undef,N)

	ul(i)=zeros(Int64,2*N)[i]
	ux(i)=ub[i]

	for k=1:N
	P1=JuMP.Model()
    nu = 1+(k-1)*2 #index of e_k
	JuMP.@variable(P1,x[i=1:2*N],upperbound=ux(i),lowerbound=ul(i), start=1.0)
	JuMP.@constraint(P1,x[nu]-gam[k]*x[nu+1]<=E[k])
    JuMP.@NLobjective(P1,Min,-x[nu]*(p[k]-0.5*x[nu])+0.5*(x[nu+1]^2)
	                         +lambda[1]*(x[1]-gam[1]*(x[2]))
							 +lambda[2]*(x[3]-gam[2]*(x[4]))
							 +lambda[3]*(x[5]-gam[3]*(x[6]))
                             +lambda[4]*(x[7]-gam[4]*(x[8]))
							 +lambda[5]*(x[9]-gam[5]*(x[10])))
	arr_player[k]=MathProgNLPModel(P1)
	end

	sol= x->((x[1]==-1.) && (x[2] == 0.))
 return GNEP(N, 2*ones(Int64,N), arr_player, sol=sol)
end

"""
Environmental protocol control proposed in
Breton, M., Zaccour, G., & Zahaf, M. (2006).
A game-theoretic formulation of joint implementation of environmental projects.
European Journal of Operational Research, 168(1), 221-239.

Variant of EPC() studied in *** new ***.
"""
function EPC5nonshared()

    include("EPC5.jl")

    ub = 50*ones(Float64,N*(N+1))
    ub[union(1,7,13,19,25)] = p

    arr_player = Array{NLPModelsJuMP.MathProgNLPModel,1}(undef,N)

	ul(i)=zeros(Float64,N*(N+1))[i]
	ux(i)=ub[i]

	for k=1:N
	P1=JuMP.Model()
        nu = 1+(k-1)*(N+1) #index of e_k
        j  = mod(k,N)+1
        l  = mod(k+1,N)+1
        o  = mod(k+2,N)+1
        m  = mod(k+3,N)+1
	JuMP.@variable(P1,x[i=1:N*(N+1)],upperbound=ux(i),lowerbound=ul(i), start=1.0)
	JuMP.@constraint(P1, dot(x[nu:(N+nu)], vcat(1.0,-gam))<=E[k])
	JuMP.@constraint(P1,-dot(x[union(1, 2, 8,14,20,26)],vcat(1.0,-gam[1]*ones(N)))
                        -dot(x[union(7, 3, 9,15,21,27)],vcat(1.0,-gam[2]*ones(N)))
                        -dot(x[union(13,4,10,16,22,28)],vcat(1.0,-gam[3]*ones(N)))
                        -dot(x[union(19,5,11,17,23,29)],vcat(1.0,-gam[4]*ones(N)))
                        -dot(x[union(25,6,12,18,24,30)],vcat(1.0,-gam[5]*ones(N)))<=0.0)
    if (k==3) || (k==4)
        JuMP.@NLconstraint(P1,-(x[16]*gam[4])/(x[14]+x[15]+x[16]+x[17]+x[18])
    						  -(x[23]*gam[3])/(x[20]+x[21]+x[22]+x[23]+x[24])<=-gap)
    end
    JuMP.@NLobjective(P1,Min,-x[nu]*(p[k]-0.5*x[nu])
                             +0.5*(x[nu+k]^2+x[nu+j]^2+x[nu+l]^2+x[nu+o]^2+x[nu+m]^2
                                   +2*x[nu+j]*x[1+(j-1)*(N+1)+j]
                                   +2*x[nu+l]*x[1+(l-1)*(N+1)+l]
                                   +2*x[nu+o]*x[1+(o-1)*(N+1)+o]
                                   +2*x[nu+m]*x[1+(m-1)*(N+1)+m])
                             +lambda[1]*(x[1]-gam[1]*(x[2]+x[8]+x[14]+x[20]+x[26]))
                             +lambda[2]*(x[7]-gam[2]*(x[3]+x[9]+x[15]+x[21]+x[27]))
                             +lambda[3]*(x[13]-gam[3]*(x[4]+x[10]+x[16]+x[22]+x[28]))
                             +lambda[4]*(x[19]-gam[4]*(x[5]+x[11]+x[17]+x[23]+x[29]))
                             +lambda[5]*(x[25]-gam[5]*(x[6]+x[12]+x[18]+x[24]+x[30])))
     arr_player[k]=MathProgNLPModel(P1)
	end

	sol= x->((x[1]==-1.) && (x[2] == 0.))
 return GNEP(N, (N+1)*ones(Int64,N), arr_player, sol=sol)
end

"""
Environmental protocol control proposed in
Breton, M., Zaccour, G., & Zahaf, M. (2006).
A game-theoretic formulation of joint implementation of environmental projects.
European Journal of Operational Research, 168(1), 221-239.

Variant of EPC() studied in *** new ***.
"""
function EPC5nonsharedgreen()

    include("EPC5.jl")

    ub = 50*ones(Float64,N*(N+1))
    ub[union(1,7,13,19,25)] = p

    arr_player = Array{NLPModelsJuMP.MathProgNLPModel,1}(undef,N)

	ul(i)=zeros(Float64,N*(N+1))[i]
	ux(i)=ub[i]

	for k=1:N
	P1=JuMP.Model()
        nu = 1+(k-1)*(N+1) #index of e_k
        j  = mod(k,N)+1
        l  = mod(k+1,N)+1
        o  = mod(k+2,N)+1
        m  = mod(k+3,N)+1
	JuMP.@variable(P1,x[i=1:N*(N+1)],upperbound=ux(i),lowerbound=ul(i), start=1.0)
	JuMP.@constraint(P1, dot(x[nu:(N+nu)], vcat(1.0,-gam))<=E[k])
	JuMP.@constraint(P1,-dot(x[union(1, 2, 8,14,20,26)],vcat(1.0,-gam[1]*ones(N)))
                        -dot(x[union(7, 3, 9,15,21,27)],vcat(1.0,-gam[2]*ones(N)))
                        -dot(x[union(13,4,10,16,22,28)],vcat(1.0,-gam[3]*ones(N)))
                        -dot(x[union(19,5,11,17,23,29)],vcat(1.0,-gam[4]*ones(N)))
                        -dot(x[union(25,6,12,18,24,30)],vcat(1.0,-gam[5]*ones(N)))<=0.0)
    if (k==1)
        #JuMP.@NLconstraint(P1,x[3]-1/(1e-3+x[7]-97.5)^2<=0.0)
        #JuMP.@NLconstraint(P1,x[4]-1/(1e-3+x[13]-97.5)^2<=0.0)
        #JuMP.@NLconstraint(P1,x[5]-1/(1e-3+x[19]-97.5)^2<=0.0)
        #JuMP.@NLconstraint(P1,x[6]-1/(1e-3+x[25]-97.5)^2<=0.0)
        JuMP.@NLconstraint(P1,x[3]*(x[7]-97.5)^2-1<=0.0)
        JuMP.@NLconstraint(P1,x[4]*(x[13]-97.5)^2-1<=0.0)
        JuMP.@NLconstraint(P1,x[5]*(x[19]-97.5)^2-1<=0.0)
        JuMP.@NLconstraint(P1,x[6]*(x[25]-97.5)^2-1<=0.0)
    end
    JuMP.@NLobjective(P1,Min,-x[nu]*(p[k]-0.5*x[nu])
                             +0.5*(x[nu+k]^2+x[nu+j]^2+x[nu+l]^2+x[nu+o]^2+x[nu+m]^2
                                   +2*x[nu+j]*x[1+(j-1)*(N+1)+j]
                                   +2*x[nu+l]*x[1+(l-1)*(N+1)+l]
                                   +2*x[nu+o]*x[1+(o-1)*(N+1)+o]
                                   +2*x[nu+m]*x[1+(m-1)*(N+1)+m])
                             +lambda[1]*(x[1]-gam[1]*(x[2]+x[8]+x[14]+x[20]+x[26]))
                             +lambda[2]*(x[7]-gam[2]*(x[3]+x[9]+x[15]+x[21]+x[27]))
                             +lambda[3]*(x[13]-gam[3]*(x[4]+x[10]+x[16]+x[22]+x[28]))
                             +lambda[4]*(x[19]-gam[4]*(x[5]+x[11]+x[17]+x[23]+x[29]))
                             +lambda[5]*(x[25]-gam[5]*(x[6]+x[12]+x[18]+x[24]+x[30])))
     arr_player[k]=MathProgNLPModel(P1)
	end

	sol= x->((x[1]==-1.) && (x[2] == 0.))
 return GNEP(N, (N+1)*ones(Int64,N), arr_player, sol=sol)
end

"""
A scalable game
"""
function NeighborsGame(N :: Int64)

    if N <= 0
        throw("Error: N must be positive.")
    end

    lb = zeros(Float64,N)
    ub = 10*ones(Float64,N)

    arr_player = Array{NLPModelsJuMP.MathProgNLPModel,1}(undef,N)

	ul(i)=lb[i]
	ux(i)=ub[i]

    C = 0.5 * rand(N,4)

	for k=1:N

        km2 = mod(k-3,N)+1
        km1 = mod(k-2,N)+1
        kp1 = mod(k,N)+1
        kp2 = mod(k+1,N)+1

    	P1=JuMP.Model()
    	JuMP.@variable(P1,x[i=1:N],upperbound=ux(i),lowerbound=ul(i), start=1.0)
    	JuMP.@constraint(P1,x[k]^2+x[km2]^2+x[km1]^2+x[kp1]^2+x[kp2]^2
                            -C[k,1]*x[k]*x[km2]-C[k,2]*x[k]*x[km1]-C[k,3]*x[k]*x[kp1]-C[k,4]*x[k]*x[kp2]<=0)
        JuMP.@NLobjective(P1,Min,x[k]^2+x[km2]^2+x[km1]^2+x[kp1]^2+x[kp2]^2-x[k]*x[kp1]*x[kp2]-x[k]*x[km1]*x[km2])
    	arr_player[k]=MathProgNLPModel(P1)
	end

	sol= x->((x[1]==-1.) && (x[2] == 0.))
 return GNEP(N, ones(Int64,N), arr_player, sol=sol)
end
