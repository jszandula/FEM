import FEM_grid as fem
import math
import numpy as np

class Elem4:
    def __init__(self):
        if fem.Grid.global_data.pc==2:
            self.ksi = [-1/math.sqrt(3), 1/math.sqrt(3), 1/math.sqrt(3), -1/math.sqrt(3)]
            self.eta = [-1/math.sqrt(3), -1/math.sqrt(3), 1/math.sqrt(3), 1/math.sqrt(3)]
            self.ksiWeight =[1,1,1,1]
            self.etaWeight = [1,1,1,1]
            self.Nksi = np.zeros((4,4))
            self.Neta = np.zeros((4,4))
            self.N = np.zeros((4,4))
            self.scheme = 2
            self.pc=4
            self.integral_points = [[self.ksi[0], 1.0], [self.ksi[1], 1.0]]
            self.int_points_weight = [1.0, 1.0]
        elif fem.Grid.global_data.pc==3:
            self.ksi = [-1*math.sqrt(3.0/5.0),0, math.sqrt(3.0/5.0), -1*math.sqrt(3.0/5.0), 0, math.sqrt(3.0/5.0), -1*math.sqrt(3.0/5.0),0, math.sqrt(3.0/5.0)]
            self.eta = [-1*math.sqrt(3.0/5.0), -1*math.sqrt(3.0/5.0), -1*math.sqrt(3.0/5.0), 0, 0, 0,math.sqrt(3.0/5.0),math.sqrt(3.0/5.0),math.sqrt(3.0/5.0)]
            self.ksiWeight = [5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0, 5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0 , 5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0]
            self.etaWeight = [5.0 / 9.0, 5.0 / 9.0, 5.0 / 9.0, 8.0 / 9.0, 8.0 / 9.0, 8.0 / 9.0 , 5.0 / 9.0, 5.0 / 9.0, 5.0 / 9.0]
            self.Nksi = np.zeros((9,4))
            self.Neta = np.zeros((9,4))
            self.N = np.zeros((9,4))
            self.scheme = 3
            self.pc=9
            self.integral_points = [[self.ksi[0], 1.0], [self.ksi[1], 1.0], [self.ksi[2], 1.0]]
            self.int_points_weight = [5.0/9.0, 8.0/9.0, 5.0/9.0]
        elif fem.Grid.global_data.pc==4:
            self.ksi = [-1 * math.sqrt((3.0/7.0) + (2.0/7.0) * math.sqrt(6.0/5.0)), -1 * math.sqrt((3.0 / 7.0) - (2.0 / 7.0) * math.sqrt(6.0 / 5.0)), math.sqrt((3.0 / 7.0) - (2.0 / 7.0) * math.sqrt(6.0 / 5.0)), math.sqrt((3.0 / 7.0) + (2.0 / 7.0) * math.sqrt(6.0 / 5.0)),
                 -1 * math.sqrt((3.0/7.0) + (2.0/7.0) * math.sqrt(6.0/5.0)), -1 * math.sqrt((3.0 / 7.0) - (2.0 / 7.0) * math.sqrt(6.0 / 5.0)), math.sqrt((3.0 / 7.0) - (2.0 / 7.0) * math.sqrt(6.0 / 5.0)), math.sqrt((3.0 / 7.0) + (2.0 / 7.0) * math.sqrt(6.0 / 5.0)),
                 -1 * math.sqrt((3.0/7.0) + (2.0/7.0) * math.sqrt(6.0/5.0)), -1 * math.sqrt((3.0 / 7.0) - (2.0 / 7.0) * math.sqrt(6.0 / 5.0)), math.sqrt((3.0 / 7.0) - (2.0 / 7.0) * math.sqrt(6.0 / 5.0)), math.sqrt((3.0 / 7.0) + (2.0 / 7.0) * math.sqrt(6.0 / 5.0)),
                 -1 * math.sqrt((3.0/7.0) + (2.0/7.0) * math.sqrt(6.0/5.0)), -1 * math.sqrt((3.0 / 7.0) - (2.0 / 7.0) * math.sqrt(6.0 / 5.0)), math.sqrt((3.0 / 7.0) - (2.0 / 7.0) * math.sqrt(6.0 / 5.0)), math.sqrt((3.0 / 7.0) + (2.0 / 7.0) * math.sqrt(6.0 / 5.0))]

            self.eta = [-1 * math.sqrt((3.0/7.0) + (2.0/7.0) * math.sqrt(6.0/5.0)), -1 * math.sqrt((3.0/7.0) + (2.0/7.0) * math.sqrt(6.0/5.0)), -1 * math.sqrt((3.0/7.0) + (2.0/7.0) * math.sqrt(6.0/5.0)), -1 * math.sqrt((3.0/7.0) + (2.0/7.0) * math.sqrt(6.0/5.0)),
                -1 * math.sqrt((3.0 / 7.0) - (2.0 / 7.0) * math.sqrt(6.0 / 5.0)), -1 * math.sqrt((3.0 / 7.0) - (2.0 / 7.0) * math.sqrt(6.0 / 5.0)), -1 * math.sqrt((3.0 / 7.0) - (2.0 / 7.0) * math.sqrt(6.0 / 5.0)), -1 * math.sqrt((3.0 / 7.0) - (2.0 / 7.0) * math.sqrt(6.0 / 5.0)),
                math.sqrt((3.0 / 7.0) - (2.0 / 7.0) * math.sqrt(6.0 / 5.0)), math.sqrt((3.0 / 7.0) - (2.0 / 7.0) * math.sqrt(6.0 / 5.0)), math.sqrt((3.0 / 7.0) - (2.0 / 7.0) * math.sqrt(6.0 / 5.0)), math.sqrt((3.0 / 7.0) - (2.0 / 7.0) * math.sqrt(6.0 / 5.0)),
                math.sqrt((3.0 / 7.0) + (2.0 / 7.0) * math.sqrt(6.0 / 5.0)), math.sqrt((3.0 / 7.0) + (2.0 / 7.0) * math.sqrt(6.0 / 5.0)), math.sqrt((3.0 / 7.0) + (2.0 / 7.0) * math.sqrt(6.0 / 5.0)), math.sqrt((3.0 / 7.0) + (2.0 / 7.0) * math.sqrt(6.0 / 5.0))]

            self.ksiWeight = [(18.0 - math.sqrt(30.0))/36, (18.0 + math.sqrt(30.0)) / 36, (18.0 + math.sqrt(30.0)) / 36, (18.0 - math.sqrt(30.0)) / 36,
                    (18.0 - math.sqrt(30.0))/36, (18.0 + math.sqrt(30.0)) / 36, (18.0 + math.sqrt(30.0)) / 36, (18.0 - math.sqrt(30.0)) / 36,
                    (18.0 - math.sqrt(30.0))/36, (18.0 + math.sqrt(30.0)) / 36, (18.0 + math.sqrt(30.0)) / 36, (18.0 - math.sqrt(30.0)) / 36,
                    (18.0 - math.sqrt(30.0))/36, (18.0 + math.sqrt(30.0)) / 36, (18.0 + math.sqrt(30.0)) / 36, (18.0 - math.sqrt(30.0)) / 36]

            self.etaWeight = [(18.0 - math.sqrt(30.0)) / 36, (18.0 - math.sqrt(30.0)) / 36, (18.0 - math.sqrt(30.0)) / 36, (18.0 - math.sqrt(30.0)) / 36,
                (18.0 + math.sqrt(30.0)) / 36, (18.0 + math.sqrt(30.0)) / 36, (18.0 + math.sqrt(30.0)) / 36, (18.0 + math.sqrt(30.0)) / 36,
                (18.0 + math.sqrt(30.0)) / 36, (18.0 + math.sqrt(30.0)) / 36, (18.0 + math.sqrt(30.0)) / 36, (18.0 + math.sqrt(30.0)) / 36,
                (18.0 - math.sqrt(30.0)) / 36, (18.0 - math.sqrt(30.0)) / 36, (18.0 - math.sqrt(30.0)) / 36, (18.0 - math.sqrt(30.0)) / 36]
            self.Nksi = np.zeros((16,4))
            self.Neta = np.zeros((16,4))
            self.N = np.zeros((16,4))
            self.scheme = 4
            self.pc=16
            self.integral_points = [[self.ksi[0], 1.0], [self.ksi[1], 1.0], [self.ksi[2], 1.0], [self.ksi[3], 1.0]]
            self.int_points_weight = [self.ksiWeight[0], self.ksiWeight[1], self.ksiWeight[2], self.ksiWeight[3]]
        self.J = [[0,0],[0,0]]
        self.detJ = 0
        self.revJ = [[0,0],[0,0]]
        self.Nx = [0,0,0,0]
        self.Ny = [0,0,0,0]
        self.finalNx = np.zeros((4,4))
        self.finalNy = np.zeros((4,4))
        self.H = np.zeros((4,4))
        self.Hbc = np.zeros((4,4))
        self.C = np.zeros((4,4))
        self.P_vector = np.zeros(4)
        self.x = []
        self.y = []
        self.globalH = np.zeros((fem.Grid.global_data.nN, fem.Grid.global_data.nN))
        self.globalC = np.zeros((fem.Grid.global_data.nN, fem.Grid.global_data.nN))
        self.globalP = np.zeros((fem.Grid.global_data.nN))
        self.alpha = fem.Grid.global_data.alfa
        self.t_ambient = fem.Grid.global_data.t_ambient
        self.dtau = fem.Grid.global_data.dtau
        
    def generateInitialShapeFunction(self):
        for i in range(0, self.pc):
            for j in range(0, 4):
                if j==0:
                    self.Nksi[i][j] = -0.25 * (1 - self.eta[i])
                if j==1:
                    self.Nksi[i][j] = 0.25 * (1 - self.eta[i])
                if j==2:
                    self.Nksi[i][j] = 0.25 * (1 + self.eta[i])
                if j==3:
                    self.Nksi[i][j] = -0.25 * (1 + self.eta[i])

        for i in range(0, self.pc):
            for j in range(0, 4):
                if j==0:
                    self.Neta[i][j] = -0.25 * (1 - self.ksi[i])
                if j==1:
                    self.Neta[i][j] = -0.25 * (1 + self.ksi[i])
                if j==2:
                    self.Neta[i][j] = 0.25 * (1 + self.ksi[i])
                if j==3:
                    self.Neta[i][j] = 0.25 * (1 - self.ksi[i])

        for i in range(0,self.pc):
            for j in range(0, 4):
                if j==0:
                    self.N[i][j] = 0.25 * (1 - self.eta[i]) * (1-self.ksi[i])
                if j==1:
                    self.N[i][j] = 0.25 * (1 - self.eta[i])* (1+self.ksi[i])
                if j==2:
                    self.N[i][j] = 0.25 * (1 + self.eta[i])* (1+self.ksi[i])
                if j==3:
                    self.N[i][j] = 0.25 * (1 + self.eta[i])* (1-self.ksi[i])

    def generateShapeFunction(self, ksi, eta):
        shape_function = [0,0,0,0]
        shape_function[0] = 0.25 * (1 - eta) * (1-ksi)
        shape_function[1] = 0.25 * (1 - eta)* (1+ksi)
        shape_function[2] = 0.25 * (1 + eta)* (1+ksi)
        shape_function[3] = 0.25 * (1 + eta)* (1-ksi)

        return shape_function

    def calculateJacobi(self, elemID, grid, pc):
        
        n = grid.Elem[elemID].ID
        self.x = (grid.ND[n[0]-1].x, grid.ND[n[1]-1].x, grid.ND[n[2]-1].x, grid.ND[n[3]-1].x)
        self.y = (grid.ND[n[0]-1].y, grid.ND[n[1]-1].y, grid.ND[n[2]-1].y, grid.ND[n[3]-1].y)
        #print(self.x) 
        #print(self.y)
        #print("\nPunkt calkowania {}".format(pc+1))
        
        for i in range(0,2):
            for j in range(0,2):
                if i==0 and j==0: 
                    self.J[i][j] = self.Nksi[pc][0]*self.x[0] + self.Nksi[pc][1]*self.x[1] + self.Nksi[pc][2]*self.x[2] + self.Nksi[pc][3]*self.x[3]
                if i==0 and j==1:
                    self.J[i][j] = self.Nksi[pc][0]*self.y[0] + self.Nksi[pc][1]*self.y[1] + self.Nksi[pc][2]*self.y[2] + self.Nksi[pc][3]*self.y[3]
                if i==1 and j==0:
                    self.J[i][j] = self.Neta[pc][0]*self.x[0] + self.Neta[pc][1]*self.x[1] + self.Neta[pc][2]*self.x[2] + self.Neta[pc][3]*self.x[3]
                if i==1 and j==1:
                    self.J[i][j] = self.Neta[pc][0]*self.y[0] + self.Neta[pc][1]*self.y[1] + self.Neta[pc][2]*self.y[2] + self.Neta[pc][3]*self.y[3]

        # print("Jakobian:")
        # for i in range(0,2):
        #     for j in range(0,2):
        #         print(self.J[i][j], end=" ")
        #     print("\n")

        self.detJ = self.J[0][0]*self.J[1][1] - self.J[0][1]*self.J[1][0]
        # print("Wyznacznik Jakobianu:\n{}\n".format(self.detJ))
                
        
    def calculateRevJacobi(self):
        self.revJ[0][0] = 1/self.detJ * self.J[1][1]
        if(self.J[0][1] == 0):
            self.revJ[0][1] = 1/self.detJ * self.J[0][1]
        else:
            self.revJ[0][1] = 1/self.detJ * -self.J[0][1]
        if(self.J[1][0] == 0):
            self.revJ[1][0] = 1/self.detJ * self.J[1][0]
        else:
            self.revJ[1][0] = 1/self.detJ * -self.J[1][0]
        self.revJ[1][1] = 1/self.detJ * self.J[0][0]

        # print("Odwrocony Jakobian:")
        # for i in range(0,2):
        #     for j in range(0,2):
        #         print(self.revJ[i][j], end=" ")
        #     print("\n")

    def calculateNxNy(self, pc):
        for i in range(0,4):
            self.Nx[i] = (self.revJ[0][0]*self.Nksi[pc][i]) + (self.revJ[0][1]*self.Neta[pc][i])

        for i in range(0,4):
            self.Ny[i] = (self.revJ[1][0]*self.Nksi[pc][i]) + (self.revJ[1][1]*self.Neta[pc][i])

        # print("Nx")
        # for i in range(0,4):
        #     print(self.Nx[i])
        
        # print("Ny")
        # for i in range(0,4):
        #     print(self.Ny[i])

        for i in range(0,4):
            for j in range(0,4):
                self.finalNx[i][j] = self.Nx[i]*self.Nx[j]

        for i in range(0,4):
            for j in range(0,4):
                self.finalNy[i][j] = self.Ny[i]*self.Ny[j]
            
        # print("d{{N}}/dx * d{{N}}T/dx")
        # for i in range(0,4):
        #     for j in range(0,4):
        #         print(self.finalNx[i][j], end=' ')
        #     print() 

        # print("d{{N}}/dy * d{{N}}T/dy")
        # for i in range(0,4):
        #     for j in range(0,4):
        #         print(self.finalNy[i][j], end=' ')
        #     print()
    
    def calculateH(self,pc, grid, id):
        for i in range(0,4):
            for j in range(0,4):
                self.H[i][j] = (self.finalNx[i][j] + self.finalNy[i][j]) * grid.Elem[id].k * self.detJ * self.ksiWeight[pc] * self.etaWeight[pc]
     
        # print("\nHmatrix")
        # for i in range(0,4):
        #     for j in range(0,4):
        #         print(self.H[i][j], end=' ')
        #     print()
        
    def calculateGlobalH(self, elemID, grid):
        for i in range(0,4):
            for j in range(0,4):
                self.globalH[grid.Elem[elemID].ID[i]-1][grid.Elem[elemID].ID[j]-1] += grid.Elem[elemID].localH[i][j]

    def calculateC(self, pc, grid, id):
        for i in range(0,4):
            for j in range(0,4):
                self.C[i][j] = (self.N[pc][i] * self.N[pc][j]) * self.detJ * grid.Elem[id].cp * grid.Elem[id].ro * self.etaWeight[pc] * self.ksiWeight[pc]

        # print("\nCmatrix")
        # for i in range(0,4):
        #     for j in range(0,4):
        #         print(self.C[i][j], end=' ')
        #     print()    
        
    def calculateGlobalC(self, elemID, grid):
        for i in range(0,4):
            for j in range(0,4):
                self.globalC[grid.Elem[elemID].ID[i]-1][grid.Elem[elemID].ID[j]-1] += grid.Elem[elemID].localC[i][j]

    def calculateBCEdges(self, elemID, grid):

        edgesBC = [None for i in range(4)]

        for i in range(4):
            if  grid.ND[grid.Elem[elemID].ID[i]-1].BC and grid.ND[grid.Elem[elemID].ID[(i+1)%4]-1].BC:
                edgesBC[i] = [grid.ND[grid.Elem[elemID].ID[i]-1], grid.ND[grid.Elem[elemID].ID[(i+1)%4]-1]]

        return edgesBC

    def calculateDetJ(self, node1, node2): 
        return math.sqrt(math.pow(node2.x - node1.x,2) + math.pow(node2.y - node1.y,2))/2

    def calculateHbcP(self, elemID, grid):
        self.Hbc = np.zeros((4,4))
        self.P_vector = np.zeros(4)
        fields_matrices = []
        p_vectors = []
        edgesBC = self.calculateBCEdges(elemID,grid)
        detJ = np.zeros((4))

        for i in range(4):
            if edgesBC[i] != None:
                detJ[i] = self.calculateDetJ(edgesBC[i][0], edgesBC[i][1])

        if edgesBC[0] != None:
            #print("Powierzchnia 1")
            field1_N = np.zeros((4,4))
            p_1 = np.zeros(4)
            for integral_point in range(self.scheme):
                ksi = self.integral_points[integral_point][0]
                eta = self.integral_points[integral_point][1]*(-1)
                weight = self.int_points_weight[integral_point]

                N = self.generateShapeFunction(ksi, eta)

                N1 = N[0]
                N2 = N[1]

                field1_N[0][0] += N1*N1*self.alpha*weight
                field1_N[0][1] += N1*N2*self.alpha*weight
                field1_N[1][0] += N2*N1*self.alpha*weight
                field1_N[1][1] += N2*N2*self.alpha*weight
                for a in range(4):
                    p_1[a] += N[a] * weight * detJ[0] *(-1)* self.alpha * self.t_ambient

            field1 = field1_N*detJ[0]
            fields_matrices.append(field1)
            p_vectors.append(p_1)

        if edgesBC[1] != None:
            #print("Powierzchnia 2")
            field2_N = np.zeros((4,4))
            p_2 = np.zeros(4)
            for integral_point in range(self.scheme):
                ksi = self.integral_points[integral_point][1]
                eta = self.integral_points[integral_point][0]
                weight = self.int_points_weight[integral_point]

                N = self.generateShapeFunction(ksi, eta)

                N2 = N[1]
                N3 = N[2]

                field2_N[1][1] += N2*N2*self.alpha*weight
                field2_N[1][2] += N2*N3*self.alpha*weight
                field2_N[2][1] += N3*N2*self.alpha*weight
                field2_N[2][2] += N3*N3*self.alpha*weight
                for a in range(4):
                    p_2[a] += N[a] * weight * detJ[1] *(-1)* self.alpha * self.t_ambient

            field2 = field2_N * detJ[1]
            fields_matrices.append(field2)
            p_vectors.append(p_2)

        if edgesBC[2] != None:
            #print("Powierzchnia 3")
            field3_N = np.zeros((4,4))
            p_3 = np.zeros(4)
            for integral_point in range(self.scheme):
                ksi = self.integral_points[self.scheme - 1 - integral_point][0]
                eta = self.integral_points[integral_point][1]
                weight = self.int_points_weight[integral_point]

                N = self.generateShapeFunction(ksi, eta)

                N3 = N[2]
                N4 = N[3]

                field3_N[2][2] += N3*N3*self.alpha*weight
                field3_N[2][3] += N3*N4*self.alpha*weight
                field3_N[3][2] += N4*N3*self.alpha*weight
                field3_N[3][3] += N4*N4*self.alpha*weight
                for a in range(4):
                    p_3[a] += N[a] * weight * detJ[2] *(-1)* self.alpha * self.t_ambient

            field3 = field3_N * detJ[2]
            fields_matrices.append(field3)
            p_vectors.append(p_3)

        if edgesBC[3] != None:
            #print("Powierzchnia 4")
            field4_N = np.zeros((4,4))
            p_4 = np.zeros(4)
            for integral_point in range(self.scheme):
                ksi = self.integral_points[integral_point][1]*(-1)
                eta = self.integral_points[self.scheme - 1 - integral_point][0] 
                weight = self.int_points_weight[integral_point]

                N = self.generateShapeFunction(ksi, eta)

                N1 = N[0]
                N4 = N[3]
        
                field4_N[0][0] += N1*N1*self.alpha*weight
                field4_N[0][3] += N1*N4*self.alpha*weight
                field4_N[3][0] += N4*N1*self.alpha*weight
                field4_N[3][3] += N4*N4*self.alpha*weight
                for a in range(4):
                    p_4[a] += N[a] * weight * detJ[3] *(-1)* self.alpha * self.t_ambient

            field4 = field4_N * detJ[3]
            fields_matrices.append(field4)
            p_vectors.append(p_4)

        for matrix in fields_matrices:
            self.Hbc +=  matrix

        for p in p_vectors:
            self.P_vector += p 

    def calculateGlobalP(self, elemID, grid):
        for i in range(4):
            self.globalP[grid.Elem[elemID].ID[i]-1] += grid.Elem[elemID].localP[i]

    def calculateNewHC(self, grid):
        for i in range(grid.global_data.nN):
            for j in range(grid.global_data.nN):
                self.globalH[i][j] = self.globalH[i][j] + (self.globalC[i][j]/self.dtau)

        for i in range(len(self.globalC)):
            row = 0.0
            for j in range(len(self.globalC[0])):
                row+= (self.globalC[i][j]/self.dtau) * grid.ND[j].t0
            self.globalP[i] = -1.0 * self.globalP[i] + row

    def clearAll(self):
        self.globalH = np.zeros((fem.Grid.global_data.nN, fem.Grid.global_data.nN))
        self.globalC = np.zeros((fem.Grid.global_data.nN, fem.Grid.global_data.nN))
        self.globalP = np.zeros((fem.Grid.global_data.nN))

    def displayGlobalH(self,grid):
        print("H Global")
        for i in range(0, fem.Grid.global_data.nN):
            for j in range(0, fem.Grid.global_data.nN):
                if(self.globalH[i][j]==0):
                    print(round(self.globalH[i][j],3), end="\t")
                else:
                    print(round(self.globalH[i][j],3), end="\t")
            print()

    def displayGlobalC(self,grid):
        print("C Global")
        for i in range(0, fem.Grid.global_data.nN):
            for j in range(0, fem.Grid.global_data.nN):
                print(round(self.globalC[i][j],3), end="\t")
            print()

    def displayGlobalP(self,grid):
        print("P Global")
        for i in range(0, fem.Grid.global_data.nN):
            print(round(self.globalP[i],3), end="\t")
        print()
