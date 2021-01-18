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
        self.dtau = 0
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
        self.k = int(data[5])
        self.cp = int(data[6])
        self.ro = int(data[7])
        self.alfa = int(data[8])
        self.t_ambient = int(data[9])
        self.dtau = int(data[10])
        self.simulation_time = int(data[11])
        self.t0 = int(data[12])
        self.nN = self.nW * self.nH
        self.nE = (self.nH - 1) * (self.nW - 1)

class Grid(Node, Element, GlobalData):
    global_data = GlobalData()
    global_data.getData()
    ND = []
    Elem = []
    Element 
    x = 0.0
    y = 0.0
    dx = float(global_data.W / (global_data.nW - 1))
    dy = float(global_data.H / (global_data.nH - 1))
    def generateNodes(self):        
        for i in range(0, self.global_data.nW):
            for j in range(0, self.global_data.nH):
                node = Node()
                node.x = i * self.dx
                node.y = j * self.dy
                if(node.x==0 or node.y==0 or node.x==self.global_data.W or node.y==self.global_data.H):
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
            nE += 1
        return self.Elem

    def printElements(self):
        print("Elements: ")
        for i in range(0, self.global_data.nE):
            print(i, self.Elem[i].ID)
    
    def printElement(self, id):
        print("Element[{}] = {}".format(id,self.Elem[id].ID))
        n1, n2, n3, n4 = self.Elem[id].ID
        print("Nodes coordinates:\n{},{},{},{}".format(self.ND[n1-1], 
                self.ND[n2-1], self.ND[n3-1], self.ND[n4-1]))

    def clearAll(self):
        for i in range(self.global_data.nE):
            self.Elem[i].localH = np.zeros((4,4))
            self.Elem[i].localC = np.zeros((4,4))
