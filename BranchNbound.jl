using JuMP
using Cbc

start=time()

Cinv = 13.16
  M = 200

m = Model(solver = CbcSolver())
        @variable(m, x[i=1:2]>=0)
        @variable(m, u, Bin)
        @objective(m, Max, 4*x[1] + 3*x[2] - u*Cinv)

        @constraint(m, 2*x[1] + 1*x[2] <= 4 +u*M)
        @constraint(m, 1*x[1] + 2*x[2] <= 4 +u*M)

        @constraint(m, 1*x[1] + 0.1*x[2] <= 4 +(1-u)*M)
        @constraint(m, 0.4*x[1] + 1*x[2] <= 4 +(1-u)*M)




# Estruturas usadas no Branck N Bound

  #Definição da direção do problema
  if getobjectivesense(model) == :Max

    direcao = 1

    ZglobalINT = -Inf

  else

    direcao = 2

    ZglobalINT = Inf

  end
  #---------------------------

  #Estruta do nó
    mutable struct node

        model::JuMP.Model

    end
  #---------------------------

  #Estruta da lista

  vectorIndex=Vector{Int}(0)

  for i = 1:length(m.colCat) #PROCURA UMA VARIÁVEL BINÁRIA

    if m.colCat[i] == :Bin

    push!(vectorIndex,i)

    end

  end

  lista = Vector{node}(0)

  push!(lista,node(model))

  model=deepcopy(m)

  #---------------------------

  #Variáveis utilizadas como Flags ou saídas do problema

  k=0

  flag=0

  Solu=0

  u

  #-------------------------------


  #BRANCH N Bound

  while length(lista) >0

    flag=0

    model=lista[1].model

    status= solve(model, relaxation = true)


    #PODA POR VIABILIDADE

    if status!=:Optimal

      shift!(lista)

      flag=1

    end

    if flag==0


      if sum(model.colVal[vectorIndex]-round(model.colVal[vectorIndex])) == 0  #PODA POR OTIMALIDADE

        Solu=1

        if direcao == 1

          if getobjectivevalue(model) > ZglobalINT

            ZglobalINT = getobjectivevalue(model)

            u=model.colVal[:]

          end

        else

          if getobjectivevalue(model) < ZglobalINT

            ZglobalINT = getobjectivevalue(model)

            u=model.colVal[:]

          end

        end

        shift!(lista)

        flag=3

      else #PODA POR LIMITE

        if direcao == 1

          if getobjectivevalue(model) < ZglobalINT

            shift!(lista)

            flag =2

          end

        else

          if getobjectivevalue(model) > ZglobalINT

            shift!(lista)

            flag =2

          end

        end

      end

      if flag==0

        j=0

        for i in vectorIndex #Ve se a variáevl é N BIn

          if m.colVal[i]-round(m.colVal[i]) !=0

                j=i

          end

        end

        shift!(lista)

        #PRIMEIRO BRANCH

        modelLeft=deepcopy(model)

        modelLeft.colLower[j]=ceil(modelLeft.colVal[j])     #ceil(modelLeft.colVal[j])

        modelLeft.colUpper[j]=ceil(modelLeft.colVal[j])      #ceil(modelLeft.colVal[j])

        NodeL=node(modelLeft)

        push!(lista,NodeL)

        #-------------------------


        #SEGUNDO BRANCH

        modelRight=deepcopy(model)

        modelRight.colLower[j]=floor(modelRight.colVal[j])     #floor(modelRight.colVal[j])

        modelRight.colUpper[j]=floor(modelRight.colVal[j])     #floor(modelRight.colVal[j])

        NodeR=node(modelRight)

        push!(lista,NodeR)

        #-------------------------


      end

    end

  end

  m.colVal=u

  fim =time()

  if Solu == 1

    println("Z convergence = ", ZglobalINT)

    println("tempo = ", fim-start)

  else

    println("Inveasible")

    println("tempo = ", fim-start)

  end
