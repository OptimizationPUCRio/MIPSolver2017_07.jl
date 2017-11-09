using JuMP
using Cbc

start=time()

model = Model(solver= CbcSolver())
@variable(model, x[i=1:3], Bin)
        @constraint(model, 6*x[1] + 5*x[2] + 5*x[3] <= 5)
        @constraint(model, x[1] == 1)
        @objective(model, Max, 6*x[1] + 4*x[2] + 3*x[3])


function SolveMIP(model :: Jump.Model)

  mutable struct node

      model::JuMP.Model

  end

  lista = Vector{node}(0)

  if getobjectivesense(model) == :Max

    direcao = 1

    Zglobal = -Inf

  else

    direcao = 2

    Zglobal = Inf

  end

  push!(lista,node(model))

  flag=0
  
  Solu =0
  
  while length(lista) >0

    flag=0

    model=lista[1].model

    status = solve(model, relaxation = true)

    #PODA POR VIABILIDADE

    if status!=:Optimal

      shift!(lista)

      flag=1

    end

    if flag==0

      if sum(model.colVal[:]-round(model.colVal[:])) == 0  #PODA POR OTIMALIDADE
        
        Solu=1
        
        if direcao == 1

          if getobjectivevalue(model) > Zglobal

            Zglobal = getobjectivevalue(model)

          end

        else

          if getobjectivevalue(model) < Zglobal

            Zglobal = getobjectivevalue(model)
            x=model.colVal[:]
          end

        end

        shift!(lista)

        flag=3

      else #PODA POR LIMITE

        if direcao == 1

          if getobjectivevalue(model) < Zglobal

            shift!(lista)

            flag =2

          end

        else

          if getobjectivevalue(model) > Zglobal

            shift!(lista)

            flag =2

          end

        end

      end

      if flag==0

        j=0

          for i = 1:length(x) #PROCURA UMA VARIÁVEL NÃO BINÁRIA (sempre pega a última)
                #Fazer para MILP (if model.colvat[i] == Bin)
            if model.colVal[i] !=1
              if model.colVal[i] !=0

                j=i

              end
            end
          end

        shift!(lista)


        #PRIMEIRO BRANCH

        modelLeft=deepcopy(model)

        modelLeft.colLower[j]=1      #ceil(modelLeft.colVal[j])

        modelLeft.colUpper[j]=1      #ceil(modelLeft.colVal[j])

        NodeL=node(modelLeft)

        push!(lista,NodeL)

        #-------------------------


        #SEGUNDO BRANCH

        modelRight=deepcopy(model)

        modelRight.colLower[j]=0      #floor(modelRight.colVal[j])

        modelRight.colUpper[j]=0      #floor(modelRight.colVal[j])

        NodeR=node(modelRight)

        push!(lista,NodeR)

        #-------------------------


      end

    end

  end
  
end
  fim =time()

#m.colval=x joga os valores no problema inicial

  if Solu == 1

    println("Z convergence = ", Zglobal)

    println("tempo = ", fim-start)

  else

    println("Inveasible")

    println("tempo = ", fim-start)

  end

