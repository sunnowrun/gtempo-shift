open("logs/convergence.log","w") do io
    redirect_stdio(stdout=io,stderr=io) do
        Main.ShiftStudy.convergence()
    end
end
