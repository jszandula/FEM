import numpy as np

class Node:
    x = 0.0
    y = 0.0
    BC = False
    t0 = 0.0
    def __init__(self):
        self.x = 0.0
        self.y = 0.0
        self.BC = False
        self.t0 = 0.0

class Element:
    k = 0.0
    cp = 0.0
    ro = 0.0
    ID = []
    localH = [[],[],[],[]]
    localC = [[],[],[],[]]
    localHbc = [[],[],[],[]]
    localP = [[],[],[],[]]
    def __init__(self):
        self.ID = np.zeros(4)
        self.localH = np.zeros((4,4))
        self.localC = np.zeros((4,4))
        self.localHbc = np.zeros((4,4))

class GlobalData:
    def __init__(self):
        self.W = 0.0
        self.H = 0.0
        self.nW = 0.0
        self.nH = 0.0

        self.nN = 0.0
        self.nE = 0.0
        self.pc = 0
        self.cp = 0.0
        self.ro = 0.0
        self.alfa = 0
        self.t_ambient = 0
        self.dtau = 0.0
        self.simulation_time = 0
        self.t0 = 0

    def getData(self):
        with open("grid.txt", 'r') as file:
            lines=file.readlines()
            data=[]   
        for x in lines:
            data.append(x.split(' ')[0])
        file.close()
        self.W = float(data[0])
        self.H = float(data[1])
        self.nW = int(data[2])
        self.nH = int(data[3])
        self.pc = int(data[4])
        self.alfa = float(data[5])
        self.t_ambient = int(data[6])
        self.dtau = float(data[7])
        self.simulation_time = int(data[8])
        self.t0 = int(data[9])
        self.nN = self.nW * self.nH
        self.nE = (self.nH - 1) * (self.nW - 1)

class Material:
    conductivity = 0.0
    specific_heat = 0.0
    density = 0.0

    def __init__(self,  k, cp, ro):
        self.conductivity = k
        self.specific_heat = cp
        self.density = ro

class Grid(Node, Element, GlobalData, Material):
    global_data = GlobalData()
    global_data.getData()
    ND = []
    Elem = []
    Element 
    x = 0.0
    y = 0.0
    dx = float(global_data.W / (global_data.nW - 1))
    dy = float(global_data.H / (global_data.nH - 1))
    glass = Material(0.8, 840, 2500)
    air = Material(0.025, 1005, 1.27)
    argon = Material(0.0177, 520, 1.784)
    def generateNodes(self):        
        for i in range(0, self.global_data.nW):
            for j in range(0, self.global_data.nH):
                node = Node()
                node.x = i * self.dx
                node.y = j * self.dy
                if(node.x==0):
                    node.BC = True
                else:
                    node.BC = False
                node.t0 = self.global_data.t0
                self.ND.append(node)
        return self.ND

    def printNodes(self):
        print("Nodes:")
        for i in range(0, self.global_data.nN):
            print(i+1, "x = " ,round(self.ND[i].x,3),"\ty = ", round(self.ND[i].y,3),"\tBC = ",self.ND[i].BC)
    
    def generateElements(self):
        for i in range(self.global_data.nE):
            self.Elem.append(Element())
        nE = 0
        for i in range(1, (self.global_data.nE+(self.global_data.nW-1))):
            if i % self.global_data.nH == 0:
                i += 1
                continue
            ID1 = i
            ID2 = i + self.global_data.nH
            ID3 = ID2 + 1
            ID4 = ID1 + 1
            self.Elem[nE].ID = (ID1, ID2, ID3, ID4)
            x_tmp = (self.ND[ID1-1].x + self.ND[ID2-1].x)/2
            if x_tmp < 0.004 or (x_tmp > 0.016 and x_tmp < 0.02) or (x_tmp > 0.032 and x_tmp < 0.036) :
                self.Elem[nE].k = self.glass.conductivity
                self.Elem[nE].cp = self.glass.specific_heat
                self.Elem[nE].ro = self.glass.density
            else:
                self.Elem[nE].k = self.argon.conductivity
                self.Elem[nE].cp = self.argon.specific_heat
                self.Elem[nE].ro = self.argon.density
            nE += 1
        return self.Elem

    def printElements(self):
        print("Elements: ")
        for i in range(0, self.global_data.nE):
            print(i, self.Elem[i].ID, "k = ", self.Elem[i].k, " cp = ", self.Elem[i].cp, " ro = ", self.Elem[i].ro)
    
    def printElement(self, id):
        print("Element[{}] = {}".format(id,self.Elem[id].ID))
        n1, n2, n3, n4 = self.Elem[id].ID
        print("Nodes coordinates:\n{},{},{},{}".format(self.ND[n1-1], 
                self.ND[n2-1], self.ND[n3-1], self.ND[n4-1]))

    def clearAll(self):
        for i in range(self.global_data.nE):
            self.Elem[i].localH = np.zeros((4,4))
            self.Elem[i].localC = np.zeros((4,4))


