

function test1(solveMIP::Function, solver::MathProgBase.AbstractMathProgSolver = JuMP.UnsetSolver())
    @testset "Custo fixo" begin
                
        #Custo Unitário
        c = [2; 3; 2; 5; 6]


        #Custo Fixo
        f = [1; 3; 1; 5; 10]

        m= Model(solver = solver)

        @variable(m, x[i=1:5]>=0)

        @variable(m,y[i=1:5], Bin)

        @constraint(m, sum(x[j] for j=1:5) >=10)

        @constraint(m,x[1]<=5*y[1])

        @constraint(m,x[2]<=4*y[2])

        @constraint(m,x[3]<=3*y[3])

        @constraint(m,x[4]<=2*y[4])

        @constraint(m,x[5]<=1*y[5])

        @objective(m, Min, sum(f[j]*y[j]+c[j]*x[j] for j=1:5))

        sol = solveMIP(m)
        @test getobjectivevalue(m) == 27
        @test getvalue(x) == [5 2 3 0 0]                
end                
        
#-------------------------------------

function test1(solveMIP::Function, solver::MathProgBase.AbstractMathProgSolver = JuMP.UnsetSolver())
    @testset "Cobertura de pontos" begin
                     
        pontosextra= 0
        
        #um ponto não é coberto por nenhum subconjunto        

        S1=[1 0 0 1 ones(1,pontosextra)]
        S2=[1 1 0 0 ones(1,pontosextra)]
        S3=[0 1 0 1 ones(1,pontosextra)]

        c=[4 3 2]

        A=[S1' S2' S3']

        m = Model(solver = solver)
        @variable(m, x[i=1:3], Bin)
        @constraints(m, begin
          constrain[i=1:4+pontosextra], sum(A[i,j]*x[j] for j=1:3)>= 1
          end)
        @objective(m, Min, sum(c[j]*x[j] for j=1:3))
                        
        m.ext[:status] == :Infeasible
end
#-------------------------------------

function test1(solveMIP::Function, solver::MathProgBase.AbstractMathProgSolver = JuMP.UnsetSolver())
    @testset "Cobertura de pontos" begin
              
        pontosextra= 50

        S1=[1 0 1 1 ones(1,pontosextra)]
        S2=[1 1 0 0 ones(1,pontosextra)]
        S3=[0 1 0 1 ones(1,pontosextra)]

        c=[4 3 2]

        A=[S1' S2' S3']

        m = Model(solver = solver)
        @variable(m, x[i=1:3], Bin)
        @constraints(m, begin
          constrain[i=1:4+pontosextra], sum(A[i,j]*x[j] for j=1:3)>= 1
          end)
        @objective(m, Min, sum(c[j]*x[j] for j=1:3))
        
        sol = solveMIP(m)
        @test getobjectivevalue(m) == 6
        @test getvalue(x) == [1 0 1] 
                                
end
                        
                        
#---------------------
                        

#Produção com custo fixo Inviavel
function test1(solveMIP::Function, solver::MathProgBase.AbstractMathProgSolver = JuMP.UnsetSolver())
    @testset "Custo fixo" begin


        vari=50


        m= Model(solver = solver)

        @variable(m, x[i=1:vari]>=0)

        @variable(m,y[i=1:vari], Bin)

        @constraint(m, sum(x[j] for j=1:vari) >= 2*vari)

        @constraints(m, begin
          constrain[i=1:vari], x[i] <= 1*y[i]
          end)

        @objective(m, Min, sum(j*y[j]+(vari-j)*x[j] for j=1:vari))
                                        
         
        sol = solveMIP(m)
        m.ext[:status] == :Infeasible  
end
#------------------------------------------
                                
                                
function test1(solveMIP::Function, solver::MathProgBase.AbstractMathProgSolver = JuMP.UnsetSolver())
    @testset "Custo fixo" begin
 
        vari=500

        m= Model(solver = solver)

        @variable(m, x[i=1:vari]>=0)

        @variable(m,y[i=1:vari], Bin)

        @constraint(m, sum(x[j] for j=1:vari) >= vari)

        @constraints(m, begin
          constrain[i=1:vari], x[i] <= 10*y[i]
          end)

        @objective(m, Min, sum(j*y[j]+(vari-j)*x[j] for j=1:vari))

        sol = solveMIP(m)
                @test getobjectivevalue(m) == 36025

                #Os últimos Y's tem que estar "ligados"
end
#-----------------------

#Expansão Unbounded
                                
function test1(solveMIP::Function, solver::MathProgBase.AbstractMathProgSolver = JuMP.UnsetSolver())
    @testset "Expansão Unbouded" begin
 

        Cinv = 13.16
           M = 200

               m = Model(solver = solver)
               @variable(m, x[i=1:2])
               @variable(m, u, Bin)
               @objective(m, Min, 4*x[1] + 3*x[2] - u*Cinv)

               @constraint(m, 2*x[1] + 1*x[2] <= 4 +u*M)
               @constraint(m, 1*x[1] + 2*x[2] <= 4 +u*M)

               @constraint(m, 1*x[1] + 0.1*x[2] <= 4 +(1-u)*M)
               @constraint(m, 0.4*x[1] + 1*x[2] <= 4 +(1-u)*M)
        sol = solveMIP(m)                                                
        m.ext[:status] == :Infeasible  #Unbounded
end                           
                        
#--------------------------
                        
function test1(solveMIP::Function, solver::MathProgBase.AbstractMathProgSolver = JuMP.UnsetSolver())
    @testset "PL da minha cabeça" begin
                                
                                
                @variable(m, x[i=1:4]>=0)

                @constraint(m, x[1]+x[2]+x[3]<=3)
                @constraint(m, x[4]+2*x[1]+6*x[3]<=10)
                @constraint(m, 4*x[3]+x[1]+3*x[2]<=5)
                                
                @objective(m, Max, 4*x[1]+5*x[2]+2*x[3]-3*x[4])

                sol = solveMIP(m)
                    @test getobjectivevalue(m) == 13
                    @test getvalue(x) == [2 1 0 0]
                                
end
                            
                            
                                                    
#--------------------------
                        
function test1(solveMIP::Function, solver::MathProgBase.AbstractMathProgSolver = JuMP.UnsetSolver())
    @testset "PL da minha cabeça" begin
                                
                                
                @variable(m, x[i=1:4]>=0)

                @constraint(m, x[1]+x[2]+x[3]<=3)
                @constraint(m, x[4]+2*x[1]+6*x[3]<=10)
                @constraint(m, 4*x[3]+x[1]+3*x[2]<=5)
                @constraint(m, x[3]==2)                
                @objective(m, Max, 4*x[1]+5*x[2]+2*x[3]-3*x[4])

                sol = solveMIP(m)                                                
                    m.ext[:status] == :Infeasible
                                
end
