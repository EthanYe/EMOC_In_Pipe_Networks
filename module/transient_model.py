import sys
sys.path.append(sys.path[0] + '\\Hydraulic Structure')
import numpy as np
from module.readconfig import *
from enum import Enum
from module.components import *
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
figure_root_path = 'figures/'
data_root_path = 'data/'
config_root_path = 'config/'

class MODE(Enum):
    """ This is the boundary set """
    MOC = 0  # Sun的value被设定为0
    SmallSignal = 1
    FiniteDiff1 = 2
    FiniteDiff2 = 3
    FiniteDiff3 = 4


class PipelineSystem(SysParam):
    # 0
    def __init__(self, name,isRecord=True,isTF=False,isKV=False,isLeak=True,id=0,saveFilePath=None):
        self.name = name
        self.id=id
        self.T = 0
        self.steady_state = False
        self.isRecord=isRecord
        self.isTF=isTF
        self.isKV=isKV
        self.isLeak=isLeak
        self.input()
        self.time_step=self.time_step
        self.steps=self.steps
        self.saveFilePath=saveFilePath
    # 1
    def input(self):

        # Read configure
        self.fjs = ReadJson(config_root_path+ self.name)
        fjs=self.fjs
        fjs.read_control_papameters()

        # read pipes
        self.pipes, self.n_pipe = fjs.read_pipes()

        # Initialize nodes, includeing id, nps,list_ps, npe,list_pe
        self.nodes = []
        self.nNodes = self.structure_nodes()
        self.init_topology_matrix()

        # read reservoirs
        self.reservoirs, self.n_reservoir = fjs.read_reservoirs()  # Get reservoirs
        for r in self.reservoirs:
            r.connect(self.nodes)
         # read nonreflectings
        self.nonreflectings, self.n_nonreflecting = fjs.read_nonreflectings()  # Get reservoirs
        for nf in self.nonreflectings:
            nf.connect(self.nodes)

        # read ballvalves
        self.ballvalves, self.n_ballvalve = fjs.read_ballvalves()  # Get ballvalves

        for bv in self.ballvalves:
            bv.connect(self.nodes)
            bv.Time_step = self.time_step
        self.n_endballvalve = BallValve.number_end  # Get ballvalve number
        # read inline valves
        self.inlineValves, self.n_inlineValve = fjs.read_inlineValves()  # Get ballvalves
        for inv in self.inlineValves:
            inv.connect(self.nodes)
        # self.n_endballvalve = BallValve.number_end  # Get ballvalve number
        # read end valves
        self.endValves, self.n_endValve = fjs.read_endValves()  # Get ballvalves
        for env in self.endValves:
            env.connect(self.nodes)
        self.coords, self.n_coord = fjs.read_coords()  # Get ballvalves
        for coor in self.coords:
            coor.connect(self.nodes)
            # tau=fft.ifft(env.complex_f)
            # plt.plot(self.t,tau)
            # plt.show()
        self.demands, self.n_demand = fjs.read_demands()  # Get ballvalves
        for dem in self.demands:
            dem.connect(self.nodes)
         # read leaks
        if self.isLeak:
            leaks,self.n_leak  = fjs.read_leaks()  # Get leaks
            for L in leaks:
                self.pipes[L.pipeID].addLeak(L)
        self.output = fjs.read_output()
        self.ini_dist()
    
    def ini_dist(self):
        # calculate total time steps
        # self.steps = round(self.total_time / self.time_step) + 1
        # read output file   
        self.recorder= []
        for pipe in self.pipes:
            pipe.steps=self.steps
            self.recorder.append(np.zeros((self.steps, pipe.NN * 2 + 3)))
            pipe.leak_flow=np.zeros((pipe.nLeak,self.steps))
            if self.isKV:
                pipe.get_visco()

        # Amount segments of pipes
        self.total_seg = 0
        for pipe in self.pipes:
            self.total_seg = self.total_seg + pipe.NN + 1
        self.N = 3 * self.total_seg   # total var
       

        # time series
        self.t = np.zeros(self.steps)
        for i in range(self.steps):
            self.t[i] = self.time_step * i
        # Get the sequence number and distance of the segments
        count = 0
        ini_distance = 0
        for pipe in self.pipes:
            pipe.count = count
            pipe.ini_distance = ini_distance
            count = count + (pipe.NN + 1) * 3  # Update the counter to the end of the pipe
            pipe.distance = np.zeros(pipe.NN + 1)
            for i in range(pipe.NN + 1):
                pipe.distance[i] = ini_distance + i * pipe.dx
            ini_distance = ini_distance + pipe.length
        self.total_length = ini_distance
    
    def update_pipe(self,pipe_id,D=None,L=None,a=None,f=None):
        pipe=self.pipes[pipe_id]
        if D is not None:
            pipe.D =D
        if L is not None:
            pipe.length=L
        if a is not None:
            pipe.a=a
        if f is not None:
            pipe.f=f
        pipe.get_intermedia_variables()
        self.ini_dist()
    
    def update_pipes(self):
        for pipe in self.pipes:
            pipe.get_intermedia_variables()
        self.ini_dist()
    # 2
    
    def structure_nodes(self):
        '''
        Define nodes that are at the start/ end of pipes. Each node has id, start pipe and end pipe (inlet or outlet).
        The nodes in the input file may be unordered, so they should be ordered by id. ID should be continuous.
        Check if there is a insular node which is not allowed.
        '''
        list_node = []
        # generate nodes
        for pipe in self.pipes:
            # add start node of the pipe
            if pipe.js not in list_node:
                list_node.append(pipe.js)  # Add node id into list
                node = Node(pipe.js)  # Create node
                self.nodes.append(node)  # Add node in nodes
            else:
                # search the start node
                for n in self.nodes:
                    if n.id == pipe.js:
                        node = n
                        break
            node.add_pipe('end pipe', pipe)
            pipe.start_node = node

            # add end node of the pipe
            if pipe.je not in list_node:
                list_node.append(pipe.je)  # Add node id into list
                node = Node(pipe.je)  # Create node
                self.nodes.append(node)  # Add node in nodes
            else:
                # search the end node
                for n in self.nodes:
                    if n.id == pipe.je:
                        node = n
                        break
            node.add_pipe('start pipe', pipe)
            pipe.end_node=node

        # sort nodes
        id_node = [(id, node) for id, node in zip(list_node, self.nodes)]  # Convert to tuple
        id_node.sort()  # Sort by id
        self.nodes = [node for score, node in id_node]
        # check nodes
        for i,node in enumerate(self.nodes):
            # node.id=node.index
            node.np = node.n_sp + node.n_ep
            node.pID = node.epID + node.spID
            node.pipes=node.start_p+node.end_p
            # if node.n_sp == 0 and node.n_ep == 1:
            #     node.position = POS.Upper
            #     node.type=BD.DeadEnd
            # elif node.n_sp != 0 and node.n_ep != 0:
            #     node.position = POS.Internal
            # elif node.n_sp == 1 and node.n_ep == 0:
            #     node.position = POS.Lower
            #     node.type = BD.DeadEnd
            # elif node.np>1:
            #     node.position = POS.Internal
            # else:
            #     error("Isolate node or multi-connected reservoir!")
            if node.np==0:
                error("Isolate node or multi-connected reservoir!")
            elif node.np == 1:
                node.type = BD.DeadEnd
            elif node.np == 2:
                node.type = BD.Series
            elif node.np > 2:
                node.type = BD.Branch

        return len(self.nodes)

    def init_topology_matrix(self):
        nNodes=self.nNodes
        # Topology matrix are used to find the pipe id between two nodes A[node 1, node 2]=pipe id (from node1 to node2)
        # topology_matrix_dire is based on the digraph
        # topology_matrix is based on the undigraph
        self.topology_matrix_dire=np.ones((nNodes, nNodes),dtype='int')*-1
        self.topology_matrix = np.ones((nNodes, nNodes), dtype='int') * -1  
        # adjacency_matrix, A[node 1, node 2]=1
        self.adjacency_matrix=np.zeros((nNodes, nNodes),dtype='int')
        # A[node 1, node 2]= number of dx 
        self.weighted_adjacency_matrix = np.ones((nNodes, nNodes), dtype='int') * np.inf
        for pipe in self.pipes:
            js=pipe.js
            je=pipe.je
            self.topology_matrix_dire[js,je]=pipe.id
            self.topology_matrix[js,je]=pipe.id
            self.topology_matrix[je,js]=pipe.id
            self.adjacency_matrix[js,je]=1
            self.adjacency_matrix[je,js]=1
            self.weighted_adjacency_matrix[js,js]=0
            self.weighted_adjacency_matrix[je,je]=0
            self.weighted_adjacency_matrix[js,je]=pipe.NN
            self.weighted_adjacency_matrix[je,js]=pipe.NN
        np.save("adjacency_matrix.npy",self.adjacency_matrix)
        a=1
    # 3
    def steady_moc(self, Q, leak=0.00):
        # Initial leak
        # self.pipes[0].K0[20] = leak
        #############
        self.T=0
        for pipe in self.pipes:
            pipe.H0[:] = self.reservoirs[0].water_level
            pipe.Q0[:] = Q
            if pipe.id==1:
                pipe.Q0[:]=0
            pipe.H[:] = pipe.H0[:]
            pipe.Q[:] = pipe.Q0[:]
            pipe.K[:] = pipe.K0[:]

        print('constant method:', self.pipes[0].H[0], self.pipes[0].Q[0])
        # # 耦合20 s的特征线法，消除恒定流和非恒定流之间的跳跃
        for i in range(5000):
            self.transient(steady=True)
        # 保存恒定流状态
        for pipe in self.pipes:
            pipe.H0[:] = pipe.H[:]
            pipe.Q0[:] = pipe.Q[:]
        print('after moc:', self.pipes[0].H0[0], self.pipes[0].Q0[0])
        # initialize Xs (X(0))
        self.XsH = self.get_state_var(var='H')
        self.XsQ = self.get_state_var(var='Q')
        self.Xs = self.get_state_var()
        self.steady_state = True
        for pipe in self.pipes:
            pipe.H0[:] = pipe.H[:]
            pipe.Q0[:] = pipe.Q[:]
            pipe.K0[:] = pipe.K[:]
        # initialize Xs (X(0))
        self.Xs = np.zeros(self.N)
        self.get_state_var(self.Xs[:])
        # As
        # self.As = np.zeros((self.N, self.N))
        # self.update_Jocob(self.As)



    def steady_unit(self, iniQ=0.00001, leak=0.00):
        IMAX = 50
        dim = self.nNodes
        QS = np.zeros((dim, dim))
        SY = np.zeros((dim, dim))
        SJ = np.zeros((dim, dim))
        S = np.zeros((dim, dim))
        E = np.zeros(dim)
        for pipe in self.pipes:
            # js = pipe.js
            # je = pipe.je
            js = pipe.start_node.id
            je = pipe.end_node.id
            QS[js, je] = iniQ
            QS[je, js] = -QS[js, je]
            # # 沿程损失流量系数
            SY[js, je] = 8.0 * pipe.f * pipe.length / self.pi / self.pi / pipe.D**5.0 / self.g
            # # on - way loss factor
            SY[je, js] = SY[js, je]
            # # 局部损失流量系数
            SJ[js, je] = pipe.zeta / 2 / self.g / pipe.A / pipe.A
            # # 流量系数
            S[js, je] = SY[js, je] + SJ[js, je]
            S[je, js] = SY[je, js] + SJ[je, js]

        # #初始化访问变量
        visit = np.zeros(dim)
        # #初始化水头，从水库开始遍历
        for res in self.reservoirs:
            id = res.node_id
            E[id] = res.water_level
            visit[id] = True
            self.Q2E(self.nodes[id], E, S, QS, visit)

        # Iteration
        for iter in range(IMAX):
            p = np.zeros((dim, dim))
            q = np.zeros((dim, dim))
            r = np.zeros((dim, dim))
            A = np.zeros((dim, dim))
            B = np.zeros(dim)
            DE = np.zeros(dim)
            # p, q, r
            for pipe in self.pipes:
                # je = pipe.je
                # js = pipe.js
                js = pipe.start_node.id
                je = pipe.end_node.id
                s = S[js, je]
                typeID = self.nodes[je].typeID
                # 计算能量损失（流量的函数）fQ及其导数
                fQ = s * abs(QS[js, je]) * QS[js, je]
                dfQ = 2 * s * abs(QS[js, je])
                if dfQ == 0:
                    dfQ = 1E-11
                if self.nodes[je].type == BD.InlineValve:
                    sj=self.inlineValves[typeID].s
                    fQ = fQ+sj * np.abs(QS[js, je] * QS[js, je])
                    dfQ = dfQ+2 * sj*np.abs(QS[js, je])
                # elif self.nodes[je].type == BD.EndValve:
                #     fQ = fQ+sj * np.abs(QS[js, je] * QS[js, je])
                #     dfQ = dfQ+2 * sj*np.abs(QS[js, je])
                # fQ,dfQ=nodes[je].get_fQ_dfQ()
                p[js, je] = 1 / dfQ
                q[js, je] = -p[js, je]
                r[js, je] = p[js, je] * (E[js] - E[je] - fQ)
                p[je, js] = p[js, je]
                q[je, js] = q[js, je]
                r[je, js] = -r[js, je]
            # 计算总体单元矩阵A和B
            # A(i,i) = 求和 pij ，A(i,j)=q(i, j)
            # B(i)=ci- 求和 Q - 求和 r
            for i in range(dim):
                for j in range(dim):
                    if i == j:
                        for k2 in range(dim):
                            A[i, i] = A[i, i] + p[i, k2]  # total pij of Ei
                    else:
                        A[i, j] = q[i, j]  # qif of Ej
                        B[i] = B[i] - QS[i, j] - r[i, j]
                if self.nodes[i].type == BD.Demand:
                    # if self.nodes[i].obj.burstTime>0:
                    #     self.nodes[i].demand=0
                    # else:
                        if E[i] - self.nodes[i].obj.z<0:
                            self.nodes[i].demand=self.nodes[i].obj.demand = self.nodes[i].obj.cda2g *100
                        else:
                            self.nodes[i].demand=self.nodes[i].obj.demand = self.nodes[i].obj.cda2gt0 *np.sqrt(E[i]- self.nodes[i].obj.z)
                if self.nodes[i].type ==BD.EndValve and self.nodes[i].obj.Q0:
                    B[i] -= self.nodes[i].obj.Q0
                B[i] -= self.nodes[i].demand
                # 水库单元方程 ΔE(i,i)=0
                if (self.nodes[i].type == BD.UpperReservoir or self.nodes[i].type == BD.LowerReservoir):  # reservoir
                    for j in range(dim):
                        if i == j:
                            A[i, j] = 1
                        else:   
                            A[i, j] = 0
                    B[i] = 0
            # for i in range(dim):
            #     if max(A[i]>1E10):
            #         A[i] /= 1E10
            #         B[i] /= 1E10
            # 解单元方程组
            # print(np.linalg.cond(A))
            """ solution 1: QR decomposition"""
            # Perform QR decomposition on A
            Q, R = np.linalg.qr(A)
            # Solve the system of linear equations
            y = np.dot(Q.T, B)
            try:
                DE = np.linalg.solve(R, y)
            except:
                print("singular matrix")
                break
            """ solution 2: normal"""
            # DE = np.linalg.solve(A, B)
            # 更新流量
            for pipe in self.pipes:
                # js = pipe.js
                # je = pipe.je
                js = pipe.start_node.id
                je = pipe.end_node.id
                QS[js, je] = QS[js, je] + p[js, je] * DE[js] + q[js, je] * DE[je] + r[js, je]  # flow of every pipe
                QS[je, js] = -QS[js, je]

            # 更新水头
            # 初始化访问变量
            visit = np.zeros(dim)
            for r in self.reservoirs:
                id = r.node_id
                E[id] = r.water_level
                visit[id] = True
                self.Q2E(self.nodes[id], E, S, QS, visit)

            # 判断结果精度
            if (np.max(np.abs(DE)) < 0.00001 and iter>3):  
                break
        if iter<IMAX-1:
            print("Steady flow completed: " + str(iter))
        else:
            print("Steady flow failed!", iter)
        # print(E,QS)
        # 利用恒定流计算结果，计算管道的水头、流量
        for pipe in self.pipes:
            # js = pipe.js
            # je = pipe.je
            js = pipe.start_node.id
            je = pipe.end_node.id
            # 管道流量
            pipe.Q[0] = QS[js, je]
            # 管道头结点的水头，能量减去速度头
            # pipe.H[0] = E[js] - pipe.Q[0] * pipe.Q[0] / 2 / self.g / pipe.A / pipe.A
            pipe.H[0] = E[js]
            # 从头结点到尾结点，依次计算管道其余节点的水头
            for j in range(1, pipe.NN+1):
                # 管道流量不变
                pipe.Q[j] = pipe.Q[0]
                # 判断流量方向
                if pipe.Q[j] > 0:
                    s = SY[js, je]
                else:
                    s = SY[je, js]
                pipe.H[j] = pipe.H[j - 1] - s * abs(pipe.Q[j]) * pipe.Q[j] / pipe.NN
            delta=E[je]-pipe.H[pipe.NN]
            if abs(delta)>0.01:
                a=1
        for node in self.nodes:
            E=[]
            for pipe in node.pipes:
                if node.id==pipe.js:
                    E.append(pipe.H[0])
                else:
                    E.append(pipe.H[pipe.NN])
            for i in E:
                delta=i-E[0]
                if abs(delta)>0.01:
                    a=1
        # print('unit method:', self.pipes[0].H[0],self.pipes[-1].H[-1], self.pipes[0].Q[0], self.pipes[-1].Q[-1])
        # 保存恒定流状态
        for pipe in self.pipes:
            pipe.H0[:] = pipe.H[:]
            pipe.Q0[:] = pipe.Q[:]
        self.XsH = self.get_state_var(var='H')
        self.XsQ = self.get_state_var(var='Q')
        print(self.XsH)
        # # 耦合20 s的特征线法，消除恒定流和非恒定流之间的跳跃
        # h_s=np.zeros(5000)
        # for i in range(5000):
        #     self.transient(steady=True)
        #     h_s[i]=self.pipes[0].HP[-1]
        #     stdq=0
        #     for pipe in self.pipes:
        #         stdq+=np.var(pipe.Q)
           
            # if stdq<1e-15:
            #     print("Stand variation: ",stdq)
            #     print("Steady MOC setps: ",i)
            #     print('----------------------------------------------------')
            #     break
        # if i>1:
        #     plt.plot(h_s[:i])
        #     plt.show()
        print('MOC steady:', self.pipes[0].H[0], self.pipes[0].Q[0])
        # 保存恒定流状态
        for pipe in self.pipes:
            pipe.H0[:] = pipe.H[:]
            pipe.Q0[:] = pipe.Q[:]
        for pipe in self.pipes:
            for j in range(pipe.NN + 1):
                pipe.b[j] = 0
                pipe.epsi[j] = 0
                for i in range(pipe.NK):  # Update b
                    pipe.bk[i, j] = pipe.km1[i] * (pipe.H[j] - pipe.H0[j]) - pipe.km2k[i] * \
                        pipe.H[j] + pipe.km3k[i] * pipe.epsik[i, j]
                    pipe.b[j] = pipe.b[j] + pipe.bk[i, j]
                    pipe.epsipk[i, j] = pipe.mm1k[i] * pipe.HP[j] - pipe.J[i] * \
                        pipe.k3 * pipe.H0[j] - pipe.tau[i] / pipe.k1 * pipe.bk[i, j]
                    pipe.epsik[i, j] = pipe.epsipk[i, j]
                    pipe.epsi[j] = pipe.epsi[j] + pipe.epsik[i, j]
        # # 保存恒定流状态
        # for pipe in self.pipes:
        #     pipe.H0[:] = pipe.H[:]
        #     pipe.Q0[:] = pipe.Q[:]
        
        # print('after moc:', self.pipes[0].H0, self.pipes[0].Q0[0])
        # initialize Xs (X(0))
        self.XsH = self.get_state_var(var='H')
        self.XsQ = self.get_state_var(var='Q')
        self.Xs = self.get_state_var()
        # np.save(self.name+'_Xs.npy',self.Xs)
        self.steady_state = True
        
        if self.isTF:
            self.x = self.get_transfer_matrix(4.189)
            self.plot_frequency_diagram()
            self.plot_frequency_response()
        
    def steady_unit_emoc(self, iniQ=0.08, leak=0.00,unknowns=[],sensors=[],sensor_data=[]):
        IMAX = 50
        dim = self.nNodes
        QS = np.zeros((dim, dim))
        SY = np.zeros((dim, dim))
        SJ = np.zeros((dim, dim))
        S = np.zeros((dim, dim))
        E = np.zeros(dim)
        for pipe in self.pipes:
            # js = pipe.js
            # je = pipe.je
            js = pipe.start_node.id
            je = pipe.end_node.id
            QS[js, je] = iniQ
            QS[je, js] = -QS[js, je]
            # # 沿程损失流量系数
            SY[js, je] = 8.0 * pipe.f * pipe.length / self.pi / self.pi / pipe.D**5.0 / self.g
            # # on - way loss factor
            SY[je, js] = SY[js, je]
            # # 局部损失流量系数
            SJ[js, je] = pipe.zeta / 2 / self.g / pipe.A / pipe.A
            # # 流量系数
            S[js, je] = SY[js, je] + SJ[js, je]
            S[je, js] = SY[je, js] + SJ[je, js]

        # #初始化访问变量
        visit = np.zeros(dim)
        # 初始化水头，从水库开始遍历
        for i,id in enumerate(sensors):
            E[id] = sensor_data[i]
            visit[id] = True
            self.Q2E(self.nodes[id], E, S, QS, visit,unknowns)

        # Iteration
        for iter in range(IMAX):
            p = np.zeros((dim, dim))
            q = np.zeros((dim, dim))
            r = np.zeros((dim, dim))
            A = np.zeros((dim, dim))
            B = np.zeros(dim)
            DE = np.zeros(dim)
            # p, q, r
            for pipe in self.pipes:
                # je = pipe.je
                # js = pipe.js
                js = pipe.start_node.id
                je = pipe.end_node.id
                s = S[js, je]
                typeID = self.nodes[je].typeID
                # 计算能量损失（流量的函数）fQ及其导数
                fQ = s * abs(QS[js, je]) * QS[js, je]
                dfQ = 2 * s * abs(QS[js, je])
                if dfQ == 0:
                    dfQ = 1E-11
                if self.nodes[je].type == BD.InlineValve:
                    sj=self.inlineValves[typeID].s
                    fQ = fQ+sj * np.abs(QS[js, je] * QS[js, je])
                    dfQ = dfQ+2 * sj*np.abs(QS[js, je])
                # elif self.nodes[je].type == BD.EndValve:
                #     fQ = fQ+sj * np.abs(QS[js, je] * QS[js, je])
                #     dfQ = dfQ+2 * sj*np.abs(QS[js, je])
                # fQ,dfQ=nodes[je].get_fQ_dfQ()
                p[js, je] = 1 / dfQ
                q[js, je] = -p[js, je]
                r[js, je] = p[js, je] * (E[js] - E[je] - fQ)
                p[je, js] = p[js, je]
                q[je, js] = q[js, je]
                r[je, js] = -r[js, je]
            # 计算总体单元矩阵A和B
            # A(i,i) = 求和 pij ，A(i,j)=q(i, j)
            # B(i)=ci- 求和 Q - 求和 r
            for i in range(dim):
                for j in range(dim):
                    if i == j:
                        for k2 in range(dim):
                            A[i, i] = A[i, i] + p[i, k2]  # total pij of Ei
                    else:
                        A[i, j] = q[i, j]  # qif of Ej
                        B[i] = B[i] - QS[i, j] - r[i, j]
                if i in unknowns:
                    for j in range(dim):
                        if j== sensors[unknowns.index(i)]:
                            A[i, j] = 1
                        else:   
                            A[i, j] = 0
                    B[i] = 0
                    continue
                if self.nodes[i].type == BD.Demand:
                    # if self.nodes[i].obj.burstTime>0:
                    #     self.nodes[i].demand=0
                    # else:
                        if E[i] - self.nodes[i].obj.z<0:
                            self.nodes[i].demand=self.nodes[i].obj.demand = self.nodes[i].obj.cda2g *100
                        else:
                            self.nodes[i].demand=self.nodes[i].obj.demand = self.nodes[i].obj.cda2gt0 *np.sqrt(E[i]- self.nodes[i].obj.z)
                if self.nodes[i].type ==BD.EndValve and self.nodes[i].obj.Q0:
                    B[i] -= self.nodes[i].obj.Q0
                B[i] -= self.nodes[i].demand
                # 水库单元方程 ΔE(i,i)=0
                if (self.nodes[i].type == BD.UpperReservoir or self.nodes[i].type == BD.LowerReservoir):  # reservoir
                    for j in range(dim):
                        if i == j:
                            A[i, j] = 1
                        else:   
                            A[i, j] = 0
                    B[i] = 0
                # if i in sensors:
                #     for j in range(dim):
                #         if i == j:
                #             A[i, j] = 1
                #         else:   
                #             A[i, j] = 0
                #     B[i] = 0
            # for i in range(dim):
            #     if max(A[i]>1E10):
            #         A[i] /= 1E10
            #         B[i] /= 1E10
            # 解单元方程组
            # print(np.linalg.cond(A))
            """ solution 1: QR decomposition"""
            # Perform QR decomposition on A
            Q, R = np.linalg.qr(A)
            # Solve the system of linear equations
            y = np.dot(Q.T, B)
            try:
                DE = np.linalg.solve(R, y)
            except:
                print("singular matrix")
                break
            """ solution 2: normal"""
            # DE = np.linalg.solve(A, B)
            # 更新流量
            for pipe in self.pipes:
                # js = pipe.js
                # je = pipe.je
                js = pipe.start_node.id
                je = pipe.end_node.id
                QS[js, je] = QS[js, je] + p[js, je] * DE[js] + q[js, je] * DE[je] + r[js, je]  # flow of every pipe
                QS[je, js] = -QS[js, je]

            # 更新水头
            # 初始化访问变量
            visit = np.zeros(dim)
            for i,id in enumerate(sensors):
                E[id] = sensor_data[i]
                visit[id] = True
                self.Q2E(self.nodes[id], E, S, QS, visit,unknowns)


            # 判断结果精度
            if (np.max(np.abs(DE)) < 0.000000001 and iter>3):  
                break
        if iter<IMAX-1:
            print("Steady flow completed: " + str(iter))
        else:
            print("Steady flow failed!", iter)
        # print(E,QS)
        # 利用恒定流计算结果，计算管道的水头、流量
        for pipe in self.pipes:
            # js = pipe.js
            # je = pipe.je
            js = pipe.start_node.id
            je = pipe.end_node.id
            # 管道流量
            pipe.Q[0] = QS[js, je]
            # 管道头结点的水头，能量减去速度头
            # pipe.H[0] = E[js] - pipe.Q[0] * pipe.Q[0] / 2 / self.g / pipe.A / pipe.A
            pipe.H[0] = E[js]
            # 从头结点到尾结点，依次计算管道其余节点的水头
            for j in range(1, pipe.NN+1):
                # 管道流量不变
                pipe.Q[j] = pipe.Q[0]
                # 判断流量方向
                if pipe.Q[j] > 0:
                    s = SY[js, je]
                else:
                    s = SY[je, js]
                pipe.H[j] = pipe.H[j - 1] - s * abs(pipe.Q[j]) * pipe.Q[j] / pipe.NN
            delta=E[je]-pipe.H[pipe.NN]
            if abs(delta)>0.01:
                a=1
        for node in self.nodes:
            E=[]
            for pipe in node.pipes:
                if node.id==pipe.js:
                    E.append(pipe.H[0])
                else:
                    E.append(pipe.H[pipe.NN])
            for i in E:
                delta=i-E[0]
                if abs(delta)>0.01:
                    a=1
        # print('unit method:', self.pipes[0].H[0],self.pipes[-1].H[-1], self.pipes[0].Q[0], self.pipes[-1].Q[-1])
        # 保存恒定流状态
        for pipe in self.pipes:
            pipe.H0[:] = pipe.H[:]
            pipe.Q0[:] = pipe.Q[:]
        for pipe in self.pipes:
            for j in range(pipe.NN + 1):
                pipe.b[j] = 0
                pipe.epsi[j] = 0
                for i in range(pipe.NK):  # Update b
                    pipe.bk[i, j] = pipe.km1[i] * (pipe.H[j] - pipe.H0[j]) - pipe.km2k[i] * \
                        pipe.H[j] + pipe.km3k[i] * pipe.epsik[i, j]
                    pipe.b[j] = pipe.b[j] + pipe.bk[i, j]
                    pipe.epsipk[i, j] = pipe.mm1k[i] * pipe.HP[j] - pipe.J[i] * \
                        pipe.k3 * pipe.H0[j] - pipe.tau[i] / pipe.k1 * pipe.bk[i, j]
                    pipe.epsik[i, j] = pipe.epsipk[i, j]
                    pipe.epsi[j] = pipe.epsi[j] + pipe.epsik[i, j]
        # # 保存恒定流状态
        # for pipe in self.pipes:
        #     pipe.H0[:] = pipe.H[:]
        #     pipe.Q0[:] = pipe.Q[:]
        
        # print('after moc:', self.pipes[0].H0, self.pipes[0].Q0[0])
        # initialize Xs (X(0))
        self.XsH = self.get_state_var(var='H')
        self.XsQ = self.get_state_var(var='Q')
        self.Xs = self.get_state_var()
        np.save(self.name+'_Xs.npy',self.Xs)
        self.steady_state = True
        
        if self.isTF:
            self.x = self.get_transfer_matrix(4.189)
            self.plot_frequency_diagram()
            self.plot_frequency_response()
        
     
        
    def Q2E(self, node, E, S, QS, visit,unknowns=[]):
        #     # # 目的：通过流量计算各节点的能量
        # # # 输入：流量，流量系数，访问标记序列
        # # # 输出：各节点的能量
        # #当前节点ID
        n1 = node.id
        # 能量损失
        # 遍历节点处的管道
        for p2 in node.pID:
            # 判断下个节点
            if self.pipes[p2].js == node.id:
                # n2 = self.pipes[p2].je
                n2 = self.pipes[p2].end_node.id
            else:
                # n2 = self.pipes[p2].js
                n2 = self.pipes[p2].start_node.id
            # 判断是否已初始化
            if visit[n2]:
                continue
            typeID = self.nodes[n2].typeID
            # 判断流量方向
            s = S[n1, n2]
            # 计算单元能量损失
            type = self.nodes[n2].type
            if n2 in unknowns:
                DH = 0
            elif type == BD.UpperReservoir:
                # 上游节点
                E[n2] = self.reservoirs[typeID].water_level
                visit[n2] = True
                continue
            elif type == BD.LowerReservoir:
                # 下游节点
                E[n2] = self.reservoirs[typeID].water_level
                visit[n2] = True
                continue
            elif type == BD.EndValve:
                # dead end
                obj = self.endValves[typeID]
                if obj.motion == 'sinusoidal' or obj.motion == 'sudden' or obj.motion == 'phase' or obj.motion == "linear" or obj.motion == 'closed'or obj.motion == 'udf':
                    Q = obj.Q0
                    DH = s * abs(Q) * Q - s * abs(QS[n1, n2]) * QS[n1, n2]
                    E[n2] = E[n1] - s * abs(QS[n1, n2]) * QS[n1, n2] - DH
                    obj.cda0 = np.sqrt(Q**2 / 2 / self.g / E[n2])
                    a=1
                elif obj.motion == 'static':
                    if obj.Q0:
                        Q = obj.Q0
                        DH = s * abs(Q) * Q - s * abs(QS[n1, n2]) * QS[n1, n2]
                        E[n2] = E[n1] - s * abs(QS[n1, n2]) * QS[n1, n2] - DH
                        obj.cda0 = np.sqrt(Q**2 / 2 / self.g / E[n2])
                    else:
                        E[n2] = E[n1] - s * abs(QS[n1, n2]) * QS[n1, n2]
                        obj.Q0=Q = obj.cda0*np.sqrt(2 * self.g * E[n2])
                        DH = s * abs(Q) * Q - s * abs(QS[n1, n2]) * QS[n1, n2]
                        E[n2] = E[n1] - s * abs(QS[n1, n2]) * QS[n1, n2] - DH
                        
                else:
                    error('Undefined dead valve')
                visit[n2] = True
            elif type==BD.InlineValve:
                DH = np.abs(QS[n1, n2]) * QS[n1, n2] *self.inlineValves[typeID].s
            else:
                DH = 0
            E[n2] = E[n1] - s * abs(QS[n1, n2]) * QS[n1, n2] - DH
            visit[n2] = True
            # 遍历下一个节点
            self.Q2E(self.nodes[n2], E, S, QS, visit,unknowns)
    # 3

    def transient(self, steady=False ,timeup=None   ):
        '''
        Moc: in-node, boundary node, update state variable and time
        '''
        if timeup is False:
            pass
        else:
            if not steady:
                self.T = self.T + self.time_step
                if self.isRecord:
                    self.record()
        # Update internal node flow and head
        for pipe in self.pipes:
            pipe.moc(self.T)

        # Update boundary node flow and head
        for node in self.nodes:
            # boundary can inherite  node 
            # self.pipes[0].Q=100
            typeID = node.typeID
            if (node.type == BD.UpperReservoir) or (node.type == BD.LowerReservoir):  # reservoir
                r = self.reservoirs[typeID]
                r.moc(self.T)
            elif node.type == BD.Series:
                node.series_moc()
            elif node.type == BD.Branch:
                node.branch_moc()
            elif node.type==BD.DeadEnd:
                node.deadend_moc()
            elif node.type == BD.EndBallValve:
                bv = self.ballvalves[typeID]
                bv.moc(self.T)
            elif node.type == BD.InlineValve:
                inv = self.inlineValves[typeID]
                inv.moc(self.T)
            elif node.type == BD.EndValve:
                env = self.endValves[typeID]
                env.moc(self.T)
            elif node.type == BD.NonReflecting:
                nf = self.nonreflectings[typeID]
                nf.moc()
            elif node.type == BD.Demand:
                dem = self.demands[typeID]
                dem.moc(self.T)
            else:
                error('Unknow node type!')

        # Update flow and head at previous time step
        for pipe in self.pipes:
            for j in range(pipe.NN + 1):
                pipe.H[j] = pipe.HP[j]
                pipe.Q[j] = pipe.QP[j]
                pipe.K[j] = pipe.KP[j]
                pipe.b[j] = 0
                pipe.epsi[j] = 0
                for i in range(pipe.NK):  # Update b
                    pipe.bk[i, j] = pipe.km1[i] * (pipe.H[j] - pipe.H0[j]) - pipe.km2k[i] * \
                        pipe.H[j] + pipe.km3k[i] * pipe.epsik[i, j]
                    pipe.b[j] = pipe.b[j] + pipe.bk[i, j]
                    pipe.epsipk[i, j] = pipe.mm1k[i] * pipe.HP[j] - pipe.J[i] * \
                        pipe.k3 * pipe.H0[j] - pipe.tau[i] / pipe.k1 * pipe.bk[i, j]
                    pipe.epsik[i, j] = pipe.epsipk[i, j]
                    pipe.epsi[j] = pipe.epsi[j] + pipe.epsik[i, j]
    
    def simple_steady(self,Q,H1,H2):
        for pipe in self.pipes:
            h1=H1-pipe.distance[0]/self.total_length*(H1-H2)
            h2=H1-pipe.distance[-1]/self.total_length*(H1-H2)
            pipe.H0[:] = np.linspace(h1,h2,pipe.NN+1)
            pipe.Q0[:] =Q
            pipe.H[:] = np.linspace(h1,h2,pipe.NN+1)
            pipe.Q[:] =Q
            self.Xs=self.get_state_var()
    # 4
    def run(self, Q=0.001, leak=0.00, mode=MODE.MOC):
        # distribute mode
        if not self.steady_state:
            self.steady_unit()
        self.T = 0
        self.inv_state_var(self.Xs[:])
        if mode == MODE.MOC:
            for i in range(self.steps):
                self.transient()
        elif mode == MODE.SmallSignal:
            self.small_signal()
        elif mode == MODE.FiniteDiff1:
            for i in range(self.steps):
                self.finite_diff(format='lax')
        elif mode == MODE.FiniteDiff2:
            for i in range(self.steps):
                self.finite_diff(format='mac')
        elif mode == MODE.FiniteDiff3:
            for i in range(self.steps):
                self.finite_diff(format='Henu')
        else:
            print('Unknow mode in run!')
            exit(0)
        # for i in self.demands:
        #     i.plot_demands()

    # 5
    def get_state_var(self, Xs=None,var=None):
        '''
        pipe(H,Q,K) -> Xs
        '''
        
        if var is None:
            if Xs is None:
                Xs = np.zeros(self.N)
            for pipe in self.pipes:  # Get Xs=X0
                for j in range(pipe.NN + 1):
                    Xs[pipe.count + 3 * j] = pipe.H[j]
                    Xs[pipe.count + 3 * j + 1] = pipe.Q[j]
                    Xs[pipe.count + 3 * j + 2] = pipe.K[j]
        elif var=='H':
            Xs = np.zeros(self.total_seg)
            for pipe in self.pipes:  # Get Xs=X0
                for j in range(pipe.NN + 1):
                    Xs[round(pipe.count/3) + j] = pipe.H[j]
        elif var == 'Q':
            Xs = np.zeros(self.total_seg)
            for pipe in self.pipes:  # Get Xs=X0
                for j in range(pipe.NN + 1):
                    Xs[round(pipe.count/3) +  j] = pipe.Q[j]

        # if self.n_ballvalve:
        #     Xs[-1] = self.bv_end.T2CDA(self.T)
        return Xs
        

    # 6
    def inv_state_var(self, Xs, p=False):
        '''
        Xs -> pipe(H,Q,K), if p, Xs -> pipe(HP, QP, KP)
        '''
        if not p:  # for second round, practice proofed it is useless
            for pipe in self.pipes:  # Get Xs=X0
                for j in range(pipe.NN + 1):
                    pipe.H[j] = Xs[pipe.count + 3 * j]
                    pipe.Q[j] = Xs[pipe.count + 3 * j + 1]
                    if Xs[pipe.count + 3 * j + 2] < 0:
                        pipe.K[j] = 0
                    pipe.K[j] = Xs[pipe.count + 3 * j + 2]
        else:
            for pipe in self.pipes:  # Get Xs=X0
                for j in range(pipe.NN + 1):
                    pipe.HP[j] = Xs[pipe.count + 3 * j]
                    pipe.QP[j] = Xs[pipe.count + 3 * j + 1]
                    if Xs[pipe.count + 3 * j + 2] < 0:
                        pipe.KP[j] = 0
                    pipe.KP[j] = Xs[pipe.count + 3 * j + 2]

    # 7
    def inv_state_QH(self, Xs):
        '''
        Xs -> pipe(H,Q)
        '''
        for pipe in self.pipes:  # Get Xs=X0
            for j in range(pipe.NN + 1):
                pipe.H[j] = Xs[pipe.count + 3 * j]
                pipe.Q[j] = Xs[pipe.count + 3 * j + 1]
                if Xs[pipe.count + 3 * j + 2] < 0:  # leak < 0, K=0
                    pipe.K[j] = 0
                # pipe.K[j] = Xs[pipe.count + 3 * j + 2]

    # 8
    def pick_data(self, var=None, h=None,q=None,pipenum=0):
        if var==None:
            if h==None:
                return self.recorder[pipenum][:, 2*q]
            elif q==None:
                if h>=0:
                    return self.recorder[pipenum][:, 2 * h + 1]
                else:
                    return self.recorder[pipenum][:, 2 * h]
        if var == 'time':
            return self.recorder[pipenum][ :, 0]
        elif var == 'H0':
           return self.recorder[pipenum][:, 1]
        elif var == 'Q0':
           return self.recorder[pipenum][:, 2]
        elif var == 'HNN':
           return self.recorder[pipenum][ :, -2]
        elif var == 'QNN':
            return self.recorder[pipenum][ :, -1]

    # 9
    def record(self):
        '''
        Record the pressure and flow at the start/end of pipes that are listed in the output (H(0), H(NN), Q(0), Q(NN))
        index(pipe number, steps, order)

        '''
        step = round(self.T / self.time_step) - 1
        for i, pipe in enumerate(self.pipes):
            NN = pipe.NN
            self.recorder[i][step,0] = self.T
            for j in range(0, NN + 1):
                self.recorder[i][step,2 * j + 1] = pipe.H[j]
                self.recorder[i][step,2 * j + 2] = pipe.Q[j]

    # 10
    def write_recorder(self, filename=None,T=0):
        startStep=int(T/self.time_step)
        if filename is None:
            filename=self.name
        else:
            filename=filename
        self.pipesData = {'pipe'+str(i): pipeData[startStep:] for i,pipeData in enumerate(self.recorder)}
        self.nodesData={}
        for node in self.nodes:
            id=node.id
            if node.id==node.pipes[0].js:
                nodeData = self.recorder[node.pipes[0].id][startStep:,1]
            else:
                nodeData = self.recorder[node.pipes[0].id][startStep:,-2]
            self.nodesData['node'+str(node.id)]=nodeData
        # self.nodesData['demandNode']=self.demands[0].node.id
        # self.nodesData['demandSize']=self.demands[0].cda2g0
        np.savez(filename+'_pipes_data.npz',**self.pipesData)
        np.savez(filename+'_nodes_data.npz',**self.nodesData)

    
