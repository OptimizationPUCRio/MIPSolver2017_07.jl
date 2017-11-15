

function test1(solveMIP::Function, solver::MathProgBase.AbstractMathProgSolver = JuMP.UnsetSolver())
    @testset "Bagulhão da P1 (não me pergunte pq)" begin
    #não                


    QtdComp = 3

    Pcomp = [20 8 7]
    Qcomp = [5 10 10]
    Demanda = 25
    #---------#

    #Definição da sua produção#
    Pmax =100
    Pmin = 0

    custo = 0

    Qmax = 6
    Qmin = 0

    Qdiscre = 3

    Prodmáx = 10

    Bm=100

    #----------#


    deltap = (Pmax - Pmin)/(2^(Qdiscre-1))

    deltaof = (Qmax - Qmin)/(2^(Qdiscre-1))

    QtdVariaveis = 2*QtdComp + 2 + 1 +4*QtdComp

    deltasOferta = [deltaof 2*deltaof 4*deltaof]

    deltasPreço = [deltap 2*deltap 4*deltap]

    #Restrição referente a discretização do Bid#
    A = [zeros(1,QtdVariaveis-Qdiscre) deltasOferta]

     #Restrição referente a discretização do Preço#
    A = [A;zeros(1,QtdVariaveis-2*Qdiscre) deltasPreço zeros(1,Qdiscre)]

    #Total produzido = Demanda#
    A = [A;ones(1,QtdComp+1) zeros(1,QtdVariaveis-(QtdComp+1));-ones(1,QtdComp+1) zeros(1,QtdVariaveis-(QtdComp+1))]

    #Produzido <= Ofertado#
    A=[A;1 zeros(1,QtdVariaveis-1-Qdiscre) -deltasOferta]

    #Produção dos Competidores <= Oferta deles#
    A=[A;zeros(QtdComp,1) eye(QtdComp) zeros(QtdComp,QtdVariaveis-QtdComp-1)]

    #Primeira restrição dual: Spot + Dual(produção vc) - Preço de oferta <= 0#
    A=[A;zeros(1,QtdComp+1) 1 -1 zeros(1,QtdComp+2*Qdiscre) -deltasPreço zeros(1,Qdiscre)]

    #Dual da produção deles#
    A=[A;zeros(QtdComp,QtdComp+1) ones(QtdComp,1) zeros(QtdComp,1) -eye(QtdComp) zeros(QtdComp,4*Qdiscre)]

    #Primal = Dual#
    DP=[Pmin Pcomp -Demanda Qmin Qcomp deltasPreço deltasOferta zeros(1,2*Qdiscre)]
    DP=[DP;-DP]
    A=[A;DP]


    #Ordem: Gvc G1 G2 G3 Spot PIvc PI1 PI2 PI3 z1 z2 z3 w1 w2 w3 x1 x2 x3 y1 y2 y3 #




    #Restrições referentes a linearização binária#
    #Se Xk = 1, Zk = produção#
    A=[A;ones(Qdiscre,1) zeros(Qdiscre,1+2*QtdComp+1) -eye(Qdiscre) zeros(Qdiscre,Qdiscre) Bm*eye(Qdiscre) zeros(Qdiscre,Qdiscre)]

    #Se Xk = 0, Zk = 0#
    A=[A;zeros(Qdiscre,1) zeros(Qdiscre,1+2*QtdComp+1) eye(Qdiscre) zeros(Qdiscre,Qdiscre) -Bm*eye(Qdiscre) zeros(Qdiscre,Qdiscre)]

    #Se yk = 1, wk = dual da produção#
    A=[A;zeros(Qdiscre,1) zeros(Qdiscre,1+QtdComp) ones(Qdiscre,1) zeros(Qdiscre,QtdComp) zeros(Qdiscre,Qdiscre) -eye(Qdiscre) zeros(Qdiscre,Qdiscre) Bm*eye(Qdiscre)]

    #Se yk = 0, wk = 0#
    A=[A;zeros(Qdiscre,1) zeros(Qdiscre,1+2*QtdComp+1) zeros(Qdiscre,Qdiscre) -eye(Qdiscre) zeros(Qdiscre,Qdiscre) -Bm*eye(Qdiscre)]

    #Bin menor que 1#
    A=[A;zeros(2*QtdComp,QtdVariaveis-2*Qdiscre) eye(2*QtdComp)]

    #Ordem: Gvc G1 G2 G3 Spot PIvc PI1 PI2 PI3 z1 z2 z3 w1 w2 w3 x1 x2 x3 y1 y2 y3 #

    #Primal#
    b = [Qmax-Qmin;Pmax-Pmin;Demanda;-Demanda;Qmin;Qcomp']

    #Dual#
    b=[b;Pmin;Pcomp';0;0]

    #Binários#
    b=[b;Bm*ones(Qdiscre,1);zeros(Qdiscre,1);Bm*ones(Qdiscre,1);zeros(Qdiscre,1)]


    #Bin <1#
    b=[b;ones(2*QtdComp,1)]



    #Ordem: Gvc G1 G2 G3 Spot PIvc PI1 PI2 PI3 z1 z2 z3 w1 w2 w3 x1 x2 x3 y1 y2 y3 #
    c = [Pmin-custo; zeros(QtdComp+1,1); Qmin; zeros(QtdComp,1); deltasPreço' ;deltasOferta' ;zeros(2*Qdiscre,1)]

    @variable(m, y[i=QtdVariaveis-2*Qdiscre+1:QtdVariaveis]>=0, Bin)

    @variable(m, dual[i=QtdComp+3:2*QtdComp+3]>=0)

    p, n = size(A)
    Astd = zeros(p,p+n)
    cstd = zeros(p+n)
    Astd = [A eye(p)]
    cstd = [c ; zeros(p)]

    #---------------------------


    p,k = size(Astd)
    @variable(m, x[i=1:QtdComp+2]>=0)
    @variable(m, xc[i=2*QtdComp+4:QtdVariaveis-2*Qdiscre]>=0)
    @variable(m, z[i=QtdVariaveis+1:k]>=0)
    @constraints(m, begin
      constrain[i=1:p], sum(Astd[i,j]*x[j] for j=1:QtdComp+2)+sum(Astd[i,k]*y[k] for k=QtdVariaveis-2*Qdiscre+1:QtdVariaveis)+sum(Astd[i,l]*z[l] for l=QtdVariaveis+1:k)+sum(Astd[i,h]*dual[h] for h=QtdComp+3:2*QtdComp+3)+sum(Astd[i,u]*xc[u] for u=2*QtdComp+4:QtdVariaveis-2*Qdiscre)<= b[i]
    end)
    @objective(m, Max, sum(cstd[j]*x[j] for j=1:QtdComp+2)+sum(cstd[k]*y[k] for k=QtdVariaveis-2*Qdiscre+1:QtdVariaveis)+sum(cstd[h]*dual[h] for h=QtdComp+3:2*QtdComp+3)+sum(cstd[u]*xc[u] for u=2*QtdComp+4:QtdVariaveis-2*Qdiscre))

        
   sol = solveMIP(m)
   @test getobjectivevalue(m) == 90  atol = exp10(-5)    
   # vc tem que produzir 4.5, confiram
        
   end     
end       
        
        
        
        
        
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
