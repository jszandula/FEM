import FEM_grid as fem
import Jacobian as jacobian
import Simulation as simulation
import numpy as np
import csv

def main():
    grid = fem.Grid() 
    grid.generateNodes()         
    grid.generateElements()
    #grid.printNodes()
    grid.printElements()
    gauss = simulation.GaussElimination()
    jacob = jacobian.Elem4()
    jacob.generateInitialShapeFunction()    # Calculate shape functions Neta, Nksi, N
    steps = int(grid.global_data.simulation_time/grid.global_data.dtau)     # Calculate number of steps: time/step_time
    f = open("mes_arg.csv", "w", newline="")
    writer = csv.writer(f)
    for step in range(steps):   
        jacob.clearAll()
        grid.clearAll()
        for j in range(0,grid.global_data.nE):
            # print("Element", j+1)
            for i in range(0,jacob.pc):
                jacob.calculateJacobi(j, grid, i)   # Calculates also determination
                jacob.calculateRevJacobi()
                jacob.calculateNxNy(i)
                jacob.calculateH(i, grid, j)
                jacob.calculateC(i, grid, j)
                for a in range(0,4):
                    for b in range(0,4):
                        grid.Elem[j].localH[a][b] += jacob.H[a][b]
                        grid.Elem[j].localC[a][b] += jacob.C[a][b]
            jacob.calculateHbcP(j,grid)

            grid.Elem[j].localHbc = jacob.Hbc
            grid.Elem[j].localP = jacob.P_vector

            for row in range(4):
                for column in range(4):
                    grid.Elem[j].localH[row][column] += jacob.Hbc[row][column]  # Adding boundary conditions to local H matrix
            
            jacob.calculateGlobalH(j, grid)
            jacob.calculateGlobalC(j, grid)
            jacob.calculateGlobalP(j,grid)

        jacob.calculateNewHC(grid)
        T0 = gauss.gaussElimination(jacob.globalH, jacob.globalP)
        #print(T0)
        for i in range(grid.global_data.nN):
            grid.ND[i].t0 = T0[i]

        writer.writerow([T0[2], T0[18], T0[66], T0[82], T0[130], T0[146], T0[186]])    
        print("Time: ", round((step+1)*grid.global_data.dtau,1) ,"\tT0 =  ", T0[2],"\t", T0[18],"\t", T0[66],"\t", T0[82],"\t", T0[130],"\t", T0[146],"\t", T0[186])
    
        # jacob.displayGlobalH(grid)
        # jacob.displayGlobalC(grid)
        # jacob.displayGlobalP(grid)               


if __name__ == "__main__":
    main()

