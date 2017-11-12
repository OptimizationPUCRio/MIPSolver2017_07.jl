sing JuMP
using Cbc

start=time()

# Estruturas usadas no Branck N Bound

  #Definição da direção do problema
  if getobjectivesense(m) == :Max

    direcao = 1

    ZglobalINT = -Inf

    Zbound = Inf

  else

    direcao = 2

    ZglobalINT = Inf

    Zbound = -Inf

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

  model=deepcopy(m)

  push!(lista,node(model))

  #---------------------------

  #Variáveis utilizadas como Flags ou saídas do problema

  k=0

  flag=0

  Solu=0

  bestx=0

  maxit=1

  iter = 0
  #-------------------------------


  #BRANCH N Bound

  while abs(Zbound-ZglobalINT) > exp10(-3) && iter <= maxit

    flag=0

    model=lista[1].model

    status= solve(model, relaxation = true)

    if direcao == 1

        if getobjectivevalue(model) < Zbound

          Zbound = getobjectivevalue(model)

        end
    else

        if getobjectivevalue(model) > Zbound

          Zbound = getobjectivevalue(model)

        end

    end

    println("LB = ", getobjectivevalue(model))

    #PODA POR VIABILIDADE

    if status!=:Optimal

      shift!(lista)

      flag=1

    end

    if flag==0


      if sum(model.colVal[vectorIndex])-sum(round.(model.colVal[vectorIndex])) == 0  #PODA POR OTIMALIDADE

        Solu=1

        if direcao == 1

          if getobjectivevalue(model) > ZglobalINT

            ZglobalINT = getobjectivevalue(model)

            bestx=model.colVal[:]

          end

        else

          if getobjectivevalue(model) < ZglobalINT

            ZglobalINT = getobjectivevalue(model)

            bestx=model.colVal[:]

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

          if model.colVal[i]-round(model.colVal[i]) !=0

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
    println("Z convergence = ", ZglobalINT)
    iter=iter+1
  end

  fim =time()

  if Solu == 1

    m.colVal=bestx

    println("Z convergence = ", ZglobalINT)

    println("tempo = ", fim-start)

  else

    println("Inveasible")

    println("tempo = ", fim-start)

  end
