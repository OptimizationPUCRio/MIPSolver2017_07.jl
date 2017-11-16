using JuMP
  using Gurobi

mutable struct node

    model::JuMP.Model

    Level::Integer

end

function SolveMIP(m::JuMP.Model)
  start=time()

    # Estruturas usadas no Branck N Bound

      #Estruta da lista
      if getobjectivesense(m) == :Max

        direcao = 1

        ZglobalINT = -Inf

        Zbound = Inf

        Zbound2 = Inf
      else

        direcao = 2

        ZglobalINT = Inf

        Zbound = -Inf

        Zbound2 = -Inf
      end
      #---------------------------

      #Estruta do nó

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

      push!(lista,node(model,0))

      #---------------------------

      #Variáveis utilizadas como Flags ou saídas do problema

      Level=0

      k=0

      flag=0

      Solu=0

      bestx=0

      maxit=200

      iter = 0
  
      jg_time0 = time_ns()

      #-------------------------------


      #BRANCH N Bound

      while abs(Zbound2-ZglobalINT) > exp10(-3) && iter <= maxit && length(lista)>0 && (time_ns()-jg_time0)/1e9 < 180

        flag=0

        model=lista[1].model

        Status= solve(model, relaxation = true)

        if iter==0

          Zbound = getobjectivevalue(model)

        end

        if direcao == 1

            if getobjectivevalue(model) < Zbound && lista[1].Level == Level

              Zbound = getobjectivevalue(model)

            elseif  lista[1].Level != Level

              Zbound2 = Zbound

              Zbound =Inf

              Level= Level +1

            end
        else

            if getobjectivevalue(model) > Zbound && lista[1].Level == Level

              Zbound = getobjectivevalue(model)

            elseif  lista[1].Level != Level

              Zbound2 = Zbound

              Zbound =-Inf

              Level= Level +1

            end

        end

        #PODA POR VIABILIDADE

       if Status ==:Infeasible

          shift!(lista)

          flag=1

          if iter == 0

            m.ext[:status] = :Infeasible

            return

          end
          
        end

        if Status ==:Unbounded

          shift!(lista)

          flag=1

          if iter == 0

            m.ext[:status] = :Unbounded

            return

          end
          
        end

        if flag==0

          if abs(sum(model.colVal[vectorIndex])-sum(round.(model.colVal[vectorIndex]))) <= 0.01  #PODA POR OTIMALIDADE
          
            Solu=Solu+1
            
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


            #PRIMEIRO BRANCH

            modelRight=deepcopy(model)

            modelRight.colLower[j]=floor(modelRight.colVal[j])     #floor(modelRight.colVal[j])

            modelRight.colUpper[j]=floor(modelRight.colVal[j])     #floor(modelRight.colVal[j])

            NodeR=node(modelRight,lista[1].Level+1)

            push!(lista,NodeR)

            #-------------------------


            #SEGUNDO BRANCH

            modelLeft=deepcopy(model)

            modelLeft.colLower[j]=ceil(modelLeft.colVal[j])     #ceil(modelLeft.colVal[j])

            modelLeft.colUpper[j]=ceil(modelLeft.colVal[j])      #ceil(modelLeft.colVal[j])

            NodeL=node(modelLeft,lista[1].Level+1)

            push!(lista,NodeL)

            #-------------------------

            shift!(lista)

          end

        end

        iter=iter+1
      end

      fim =time()

      m.ext[:time] = fim-start
      m.ext[:intsols] = Solu
      m.ext[:nodes] = iter
      
      if Solu != 0

        m.ext[:status] = :Optimal
        
        m.colVal=bestx

        m.objVal = ZglobalINT
        
        println("Iterações= ", iter)

        println("Z convergence = ", ZglobalINT)

        println("tempo = ", fim-start)

      else

        m.ext[:status] = :Infeasible
        
        println("Iterações= ", iter)

        println("Infeasible")

        println("tempo = ", fim-start)

      end
end
