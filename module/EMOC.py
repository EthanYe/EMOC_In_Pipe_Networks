import pickle
from module.transient_model import *
from module.components import *
from module.myfun import *
from module.Floyd import *
import os.path as osp

class EMOC():
    def __init__(self, sensors, unknowns, name='EMOC_model', dataname=None,  model=None,mode=None,noise=0, parameter_changed_pipes=None, sensor_data=None,ini_data=None,ini_data_name=None,
                 sensor_data_name=None,extra_sensor=None,extra_sensor_data=None):
        self.name = name # name of pipe network model
        self.sensors=sensors # sensor locations
        self.unknowns=unknowns # unknown boundary conditions
        self.dataname = dataname # name of the input data, default is None
        self.mode = mode # initial condition method, default is None (use MOC results)
        self.sensor_data = sensor_data  # input data, default is [] (use MOC results)
        self.ini_data=ini_data # input data, default is [] (use MOC results)
        self.ini_data_name =ini_data_name
        self.sensor_data_name = sensor_data_name
        self.extra_sensor = extra_sensor #
        self.extra_sensor_data = extra_sensor_data
        self.isSteady=False
        # load pipe network model
        if model is None:
            a_file = open(data_root_path + self.name + "_class.pkl", "rb")
            self.model = model = pickle.load(a_file)
        else:
            self.model = model
        # get the shorest paths from sensors to unknown boundaries
        self.node_routes=node_routes=shortest_paths(model.weighted_adjacency_matrix,sensors,unknowns)
        self.sensor_name=''
        for i in self.sensors:
            self.sensor_name = self.sensor_name+str(i)
        # create uncertainty by changing the parameters of specific pipes (default is None)
        # if len(parameter_changed_pipes):
           
        #     model.pipes[3].D = parameter_changed_pipes[0]
        #     model.pipes[3].get_intermedia_variables()
        #     model.pipes[4].D = parameter_changed_pipes[1]
        #     model.pipes[4].get_intermedia_variables()
        self.noise = noise # the noise added into sensor data, default is 0
        # number of pipes
        self.nPipes = self.model.n_pipe
        # get pipe object
        self.pipes = model.pipes[:]
        # get node object
        self.nodes = model.nodes[:]
        # total time steps
        self.steps = model.steps
        # the maximum time to calculte the first front, not exceed the total pipe reaches of the system
        self.max_init_time = model.total_seg
        # number of sensors
        self.nsensor = len(sensors)
        # get pipe routes by node routes
        self.node_routes_to_pipe_routes()
        # define route direction, pipe.dire=1 if the path is from start node to end node, else pipe.dire=-1
        self.define_route_direction()
        # input sensor data
        # if sensor_data is None:
        self.input_data()

    def node_routes_to_pipe_routes(self):
        '''find the passing pipes for each node route'''
        self.pipe_routes = [] # list of pipe_routes
        self.first_pipe_routes=[] # list of the first pipes of pipe_routes
        for node_route in self.node_routes:
            pipe_route=[]
            for i in range(len(node_route)-1):
                node=self.nodes[node_route[i]]
                for pipe in node.pipes:  
                    if (node_route[i+1]==pipe.js or node_route[i+1]==pipe.je) and (node_route[i+1]!=node_route[i]):
                        # add the pipe between node_route[i] and node_route[i+1] into pipe route
                        pipe_route.append(pipe.id)
                        # define the pipe as nextPipe of node_route[i]
                        node.nextPipe=pipe 
            self.pipe_routes.append(pipe_route)
            if len(pipe_route)>0:
                self.first_pipe_routes.append(pipe_route[0])
            else:
                self.first_pipe_routes.append([])
    
    def define_route_direction(self):
        ''' define the route direction of the pipe'''
        '''go from the start to the end of the pipe, pipe.dire=1, otherwise, pipe.dire=-1'''
        pipes=self.pipes
        for i, pipe_route in enumerate(self.pipe_routes):
            for j in range(len(pipe_route)):
                if self.node_routes[i][j] == pipes[pipe_route[j]].js:
                    # from the start to the end of the pipe
                    pipes[pipe_route[j]].dire = 1
                    # if self.node_routes[i][j]==self.node_routes[i][j+1]:
                    #     pipes[pipe_route[j]].dire = -1
                else:
                    
                    pipes[pipe_route[j]].dire = -1
                    # if self.node_routes[i][j]==self.node_routes[i][j+1]:
                    #     pipes[pipe_route[j]].dire = 1
  
    def input_data(self):
        '''Read simulated MOC data into self.data. Store the sensor data in self.sensor_data'''
        # read simulated data
        steps = self.steps
        self.data = []
        nsensor = self.nsensor
        # fetch sensor data
        steps = self.steps
        
        # The data has the list type of shape [pipe number][2*pipe reaches][time steps]. The 2nd dim of length 2*pipe reaches represents [head 1, flow 1, head 2, flow 2, ..., ...]
        
        if self.dataname != None:
            
            loaded_data=np.load(data_root_path + self.dataname + '_pipes_data.npz')
            for i in range(self.model.n_pipe):
                array_name = f'pipe{i}'
                loaded_array = loaded_data[array_name][:,1:]
                self.data.append(loaded_array)
        elif self.sensor_data is not None:
            return
        else:
            loaded_data=np.load(data_root_path + self.name+'_pipes_data.npz')
            for i in range(self.model.n_pipe):
                array_name = f'pipe{i}'
                loaded_array = loaded_data[array_name][:,1:]
                self.data.append(loaded_array)
        if self.sensor_data is not None:
        #     self.steps=self.model.steps=self.sensor_data.shape[1]
            if self.mode == 'delta':
                for i in range(nsensor):
                    self.sensor_data[i]-=self.sensor_data[i][0]
            return
 
        self.sensor_data = np.zeros((nsensor, steps))
        for i in range(nsensor):
            if len(self.pipe_routes[i])==0:
                for pipe in self.pipes:
                    if pipe.js==self.node_routes[i][0]:
                        self.sensor_data[i] = self.data[pipe.id
                                                        ][:, 0][:steps] + self.noise[i]
                    elif pipe.je==self.node_routes[i][0]:
                        self.sensor_data[i] = self.data[pipe.id
                                                        ][:, -2][:steps] +  self.noise[i]
            
            elif self.sensors[i] == self.pipes[self.pipe_routes[i][0]].js:
                self.sensor_data[i] = self.data[self.pipe_routes[i][0]
                                                ][:, 0][:steps] + np.random.random(steps) * self.noise
            else:
                self.sensor_data[i] = self.data[self.pipe_routes[i][0]][:, -
                                                                        2][:steps] + np.random.random(steps) * self.noise
            if self.mode == 'delta':
                self.sensor_data[i][0]*=1.05
                self.sensor_data[i]-=self.sensor_data[i][0]
    
    def emoc(self):
        '''The EMOC method is used to reconstruct the flow and head based on the sensor data
            The reconstructed state is stored in pipe.qi and pipe.hi, which have the shape [total time steps, pipe reaches]'''
        # get the front by using Nan at unknown boundaries as tracers
        self.delta_new_front = np.zeros(self.nsensor, dtype="int")
        self.find_front()
        self.print_front("1st find front")
        # update front if the sensor nodes and branch nodes are failed to meet the requirements
        self.update_front()
        
        # input the data of initial condition
        self.input_initial_condition()
        
        self.front_checking()
        # exit(0)
        # initialize the head and flow in the known zone (beneath the first front)
        
        # get the maximum and minimum levels of the front
        self.min_front, self.max_front = self.min_max_front()
        # the maximum time steps for plotting
        self.plotSteps = self.steps - (self.max_front - self.min_front)
        '''EMOC loop body, including Step 1, Step 2, and Step 3'''
        if (not self.isSteady) and (self.mode==''):
            # self.sensor_data[0][0]=0
            # self.sensor_data[1][0]=0
            # stdh=0
            # stdq=0
            # for pipe in self.pipes:
            #         stdh+=np.var(pipe.hi[0])
            #         stdq+=np.var(pipe.qi[0])
            # print("Stand variation: ", stdh,stdq)
            # print('----------------------------------------------------')
            # if stdq>1e-25:  
            for i in range(30):
                self.initialize_known_zone(steady=True)
                self.emoc_loop(steady=True)
                stdh=0
                stdq=0
                for pipe in self.pipes:
                    stdh+=np.var(pipe.hi[self.plotSteps-1])
                    stdq+=np.var(pipe.qi[self.plotSteps-1])
                    pipe.qi[0, :] =pipe.qi[self.plotSteps-1, :]
                    pipe.hi[0, :] = pipe.hi[self.plotSteps-1, :] # ::2 is the index of head
                print('The '+str(i+1)+'st steady state calculation completed')
                print("Stand variation: ", stdh,stdq)
                print('----------------------------------------------------')
                if stdq<1e-9:
                    break
            if i>1:
                self.saveInitialConditions()
            self.isSteady = True
        
        self.initialize_known_zone()
        self.emoc_loop()
    
    def saveInitialConditions(self):
        self.pipesDataH = {}
        self.pipesDataQ = {}
        
        for pipe in self.pipes:

            pipeDataQ = pipe.qi[0]
            pipeDataH = pipe.hi[0]

            self.pipesDataH['pipeH'+str(pipe.id)]=pipeDataH
            self.pipesDataQ['pipeQ'+str(pipe.id)]=pipeDataQ
        # self.nodesData['demandNode']=self.demands[0].node.id
        # self.nodesData['demandSize']=self.demands[0].cda2g0
        np.savez(data_root_path+self.name+self.sensor_name+'_pipes_steady_dataH.npz',**self.pipesDataH)
        np.savez(data_root_path+self.name+self.sensor_name+'_pipes_steady_dataQ.npz',**self.pipesDataQ)

    def emoc_loop(self,steady=False):  
        model = self.model
        pipes = self.pipes
        nPipes = self.nPipes
        sensor_data = self.sensor_data
        for ti in range(1, self.steps - self.max_front):
            for kk, pipe_route in enumerate(self.pipe_routes):
                '''Step 1: update sensor node'''
                # sensor node
                node = self.nodes[self.sensors[kk]] 
                # if the sensor node is at a boundary
                if len(self.pipe_routes[kk]) == 0:
                    continue
                else:
                    # the first SB pipe
                    pipe = pipes[pipe_route[0]]
                # get the next front at sensors
                if pipe.dire == 1:  # the sensor is at the start of the pipe
                    # the next front
                    front_t = pipe.front_t[0] + ti
                else:  # the sensor is at the end of the pipe
                    front_t = pipe.front_t[-1] + ti
                # calculate the flow at the sensor node
                if steady:
                    sensor_h = sensor_data[kk][0]
                else:
                    sensor_h = sensor_data[kk][front_t]
                if node.type == BD.Series:
                    
                    if node.pipes[0]==pipe:
                        previousPipe = node.pipes[1]
                    else:
                        previousPipe = node.pipes[0]
                    if pipe.dire == 1:  # sensor -> unknown pipe 
                        pipe.hi[front_t, 0] = sensor_h
                        if previousPipe.je==node.id:   # known pipe (previousPipe) -> sensor -> unknown pipe (pipe2),node.n_sp=1,node.n_p=1
                            previousPipe.hi[front_t, previousPipe.NN] = pipe.hi[front_t, 0]
                            NN = previousPipe.NN
                            # C+ line
                            previousPipe.qi[front_t, NN] = previousPipe.QCPi(
                                front_t, NN) - previousPipe.CQPi(front_t, NN) * sensor_h
                            pipe.qi[front_t, 0] = previousPipe.qi[front_t,NN]
                        else: # known pipe (previousPipe) <- sensor -> unknown pipe (pipe2),node.n_ep=2
                            previousPipe.hi[front_t, 0] = pipe.hi[front_t, 0]
                            # C- line
                            previousPipe.qi[front_t, 0] = previousPipe.QCMi(
                                front_t, 0) + previousPipe.CQMi(front_t, 0) * sensor_h
                            pipe.qi[front_t, 0] = -previousPipe.qi[front_t,0]
                    elif pipe.dire == -1:  # sensor <- unknown pipe
                        pipe.hi[front_t, pipe.NN] = sensor_h
                        if previousPipe.js==node.id:      # known pipe (previousPipe) <- sensor <- unknown pipe (pipe2),node.n_sp=2
                            previousPipe.hi[front_t, 0] = sensor_h
                            # C- line
                            previousPipe.qi[front_t, 0] = previousPipe.QCMi(
                                front_t, 0) + previousPipe.CQMi(front_t, 0) * sensor_h
                            pipe.qi[front_t, pipe.NN] = previousPipe.qi[front_t, 0]
                        else: # known pipe (previousPipe) -> sensor <- unknown pipe (pipe2)
                            previousPipe.hi[front_t, previousPipe.NN] = sensor_h
                            # C- line
                            previousPipe.qi[front_t,  previousPipe.NN] = previousPipe.QCPi(
                                front_t,  previousPipe.NN) - previousPipe.CQPi(front_t,  previousPipe.NN) * sensor_h
                            pipe.qi[front_t, pipe.NN] = -previousPipe.qi[front_t, previousPipe.NN]
                elif node.type == BD.Branch:
                    node.cross_branch_i(front_t, sensor_h)
                elif node.type == BD.UpperReservoir:
                    pipe = node.end_p[0]
                    pipe.hi[front_t, 0] = sensor_h
                    pipe.qi[front_t, 0] = pipe.QCMi(front_t, 0) + pipe.CQMi(front_t, 0) * sensor_h
                elif node.type == BD.LowerReservoir:
                    pipe = node.start_p[0]
                    NN = pipe.NN
                    pipe.hi[front_t, NN] = sensor_h
                    pipe.qi[front_t, NN] = pipe.QCPi(front_t, NN) - pipe.CQPi(front_t, NN) * sensor_h
                else:
                    print('Undefined sensor node!')

                '''Step 2: update sensor-boundary pipes'''
                for j in range(len(pipe_route)):  # from pipe[sensors[0]-1]to pipe 0
                    pipe = pipes[pipe_route[j]]
                    self.update_pipe(pipe, ti)
                    for xx in range(pipe.NN):
                        if np.isnan(pipe.hi[pipe.front_t[xx]+ti,xx]) and pipe.front_t[xx]+ti>=0:
                            print("Illegal values in Step 2 for pipe "+ str(pipe.id)+' at steps '+str(ti))
                            exit(0)
                            pass
                        
            '''Step 3: update front in non-sensor-boundary pipes using MOC'''
            for i in range(self.min_front, self.max_front + 1):
                front_t = ti + i
                # skip the points where the front <=0
                if front_t < 1:
                    continue
                # update boundary nodes
                for node in model.nodes:
                    # update sensor node, excluding the nodes in the node routes
                    # because the sensor node has different front for each pipe. Only the front in the node routes has been updated. Now the fronts in the other pipes are updated here.
                    if node.id in self.sensors:
                        for pipe in node.pipes:
                            # skip the pipe in pipe routes
                            if pipe.id in self.first_pipe_routes:
                                continue
                            # sensor node
                            sensor_ind = self.sensors.index(node.id) # index of this sensor node
                            if pipe.front_t.max() == i:
                                if node.id == pipe.je:
                                    sp = pipe
                                    if steady:
                                        sensor_h = self.sensor_data[sensor_ind][0]
                                    else:
                                        sensor_h = self.sensor_data[sensor_ind][front_t]
                                    # sensor_h = self.sensor_data[sensor_ind][front_t]
                                    sp.hi[front_t, sp.NN] = sensor_h
                                    sp.qi[front_t, sp.NN] = sp.QCPi(front_t, sp.NN) - \
                                        sp.CQPi(front_t, sp.NN) * sp.hi[front_t, sp.NN]
                                    if np.isnan(pipe.hi[front_t,  sp.NN] ) or np.isnan(pipe.qi[front_t,  sp.NN] ):
                                        pass
                                else:
                                    ep = pipe
                                    if steady:
                                        sensor_h = self.sensor_data[sensor_ind][0]
                                    else:
                                        sensor_h = self.sensor_data[sensor_ind][front_t]
                                    # sensor_h = self.sensor_data[sensor_ind][front_t]
                                    ep.hi[front_t, 0] = sensor_h
                                    ep.qi[front_t, 0] = ep.QCMi(front_t, 0) + ep.CQMi(front_t, 0) * ep.hi[front_t, 0]
                                    if np.isnan(pipe.hi[front_t, 0] ) or np.isnan(pipe.qi[front_t, 0] ):
                                          a=1
                                          
                        a=1
                    elif node.id in self.unknowns:
                        pass
                    else:
                        # if node is not Sensor nodes:
                        if node.front_t.min() == i:
                            if node.type == BD.Series:
                                node.series_i(front_t)
                            elif node.type == BD.Branch:
                                node.branch_i(front_t)
                            elif node.type == BD.DeadEnd:
                                node.deadend_i(front_t)
                            elif node.type == BD.Demand:
                                if node.id in self.unknowns:
                                    node.unknown_demand_i(front_t)
                                else:
                                    node.demand_i(front_t)
                            elif node.type == BD.InlineValve:
                                node.inlineValve_i(front_t)
                            # elif node.type == BD.EndValve:
                            #     node.endValve_i(front_t)
                            else:
                                error("Undefined node type in non-SB pipe")
                            # update known reservoir boundaries
                            if node.id not in self.unknowns:
                                if node.type == BD.UpperReservoir or node.type == BD.LowerReservoir:
                                    typeID = node.typeID
                                    model.reservoirs[typeID].moc_i(front_t)
                                    if self.mode == "delta":
                                        model.reservoirs[typeID].moc_i(front_t, waterLevel=0)
                            # elif node.type == BD.Demand:
                            #     node.known_demand_i(front_t)
                            # elif node.type == BD.UpperReservoir:

                        a=1
                    # if node.id in self.extra_sensor:
                    #     if node.front_t.min() == i:
                    #         for pipe in node.pipes:
                    #             if node.id == pipe.je:
                    #                 sp = pipe
                    #                 delta=pipe.hi[front_t, sp.NN]  -self.extra_sensor_data[front_t]
                    #                 pipe.hi[front_t, sp.NN] =pipe.hi[front_t, sp.NN] -delta
                    #             elif node.id == pipe.je:
                    #                 ep=pipe
                    #                 delta=pipe.hi[front_t, 0]  -self.extra_sensor_data[front_t]
                    #                 pipe.hi[front_t, 0] =pipe.hi[front_t, 0] -delta
                a=1
                # update the flow and head at interal points
                for k in range(nPipes):
                    pipe = pipes[k]
                    for j in range(1, pipe.NN):
                        if pipe.front_t[j] == i:
                            
                            pipe.hi[front_t, j] = (pipe.QCPi(front_t, j) - pipe.QCMi(front_t, j)) / \
                                (pipe.CQPi(front_t, j) + pipe.CQMi(front_t, j))
                            pipe.qi[front_t, j] = pipe.QCPi(front_t, j) - pipe.CQPi(front_t, j) * pipe.hi[front_t, j]
                            if np.isnan(pipe.hi[front_t, j] ) or np.isnan(pipe.qi[front_t, j] ):
                                print("Illegal values in Step 3 for pipe "+ str(pipe.id))
                                exit(0)
                                pass
            a=1
            for pipe in self.pipes:
                for xx in range(pipe.NN+1):
                    if np.isnan(pipe.hi[pipe.front_t[xx]+ti,xx]) and pipe.front_t[xx]+ti>=0:
                        print("Illegal values in Step 3 for pipe "+ str(pipe.id))
                        exit(0)
                        pass
    
    def input_initial_condition(self):
        
        '''input initial condition'''
        if self.ini_data is not None:
            # self.model.inv_state_var(self.ini_data)
            for i in range(self.nPipes):
                # i is the pipe index, 0 is the time step, 1::2 is the index of flow
                self.pipes[i].qi[0, :] = self.ini_data[i][1::2]
                self.pipes[i].hi[0, :] = self.ini_data[i][::2] # ::2 is the index of head
                # self.pipes[i].qi[0, :] = self.model.pipes[i].Q[ :]
                # self.pipes[i].hi[0, :] =  self.model.pipes[i].H[ :] # ::2 is the index of head
        elif self.ini_data_name is not None:
            init_dataset=[]
            loaded_data=np.load(data_root_path + self.ini_data_name + '_pipes_data.npz')
            for i in range(self.model.n_pipe):
                array_name = f'pipe{i}'
                loaded_array = loaded_data[array_name][:,1:]
                init_dataset.append(loaded_array)
                self.pipes[i].qi[0, :] = loaded_array[0, 1::2]
                self.pipes[i].hi[0, :] = loaded_array[0, ::2]  # ::2 is the index of head
            # for i in range(self.nPipes):
            #     init_dataset += [np.load(data_root_path + self.ini_data_name + '_pipe' + str(i) + '.npy')[:, 1:]]
            #     self.pipes[i].qi[0, :] = init_dataset[i][0, 1::2]
            #     self.pipes[i].hi[0, :] = init_dataset[i][0, ::2]  # ::2 is the index of head
        elif self.mode == 'delta':
            for i in range(self.nPipes):
                self.pipes[i].qi[0, :] = 0
                self.pipes[i].hi[0, :] = 0 # ::2 is the index of head
            return
        elif osp.exists(data_root_path+self.name+self.sensor_name+'_pipes_steady_dataQ.npz'):
                save_file=data_root_path+self.name+self.sensor_name+'_pipes_steady_dataQ.npz'
                print("Using Cached file: {}".format(save_file))
                loaded_inidataQ=np.load(save_file)
                loaded_inidataH=np.load(save_file[:-5]+'H.npz')
                for i,pipe in enumerate(self.pipes):
                    array_nameH = f'pipeH{pipe.id}'
                    array_nameQ = f'pipeQ{pipe.id}'
                    loaded_arrayH = loaded_inidataH[array_nameH]
                    loaded_arrayQ = loaded_inidataQ[array_nameQ]
                    pipe.qi[0, :] = loaded_arrayQ
                    pipe.hi[0, :] = loaded_arrayH # ::2 is the index of head

        else: 
             
            for i in range(self.nPipes):
                # i is the pipe index, 0 is the time step, 1::2 is the index of flow
                self.pipes[i].qi[0, :] = self.data[i][0, 1::2]
                self.pipes[i].hi[0, :] = self.data[i][0, ::2]  # ::2 is the index of head
            
                # self.pipes[i].qi[0, :] = 0.001
                # self.pipes[i].hi[0, :] =30
            
    def initialize_known_zone(self,steady=False):
        ''' This function is used to calculate the state until the time steps of max_init_time'''
        model = self.model
        pipes = self.pipes
        nPipes = self.nPipes
        
        for in_t in range(1, self.max_init_time + 1):
            # calculate the head and flow at internal points
            for k in range(nPipes):
                pipe = pipes[k]
                for j in range(1, pipe.NN ):
                    pipe.hi[in_t, j] = (pipe.QCPi(in_t, j) - pipe.QCMi(in_t, j)) / \
                                        (pipe.CQPi(in_t, j) + pipe.CQMi(in_t, j))
                    pipe.qi[in_t, j] = pipe.QCPi(in_t, j) - pipe.CQPi(in_t, j) * pipe.hi[in_t, j]


            # # calculate the head and flow at boundary nodes
            for node in model.nodes:
                front_t=in_t
                typeID = node.typeID
                # if not node.isSensor:
                if node.id in self.sensors:
                    i=self.sensors.index(node.id)
                    if steady:
                        sensor_h = self.sensor_data[i][0]
                    else:
                        sensor_h = self.sensor_data[i][front_t]
                    for sp in node.start_p:
                        sp.hi[front_t, sp.NN] = sensor_h
                        sp.qi[front_t, sp.NN] = sp.QCPi(front_t, sp.NN) - sp.CQPi(front_t, sp.NN) * sp.hi[front_t, sp.NN]
                    for ep in node.end_p:
                        ep.hi[front_t, 0] = sensor_h
                        ep.qi[front_t, 0] = ep.QCMi(front_t, 0) + ep.CQMi(front_t, 0) * ep.hi[front_t, 0]
                else:
                    '''Interior boundary'''
                    if node.id not in self.unknowns:
                        if node.type == BD.Series:
                            node.series_i(in_t)
                        elif node.type == BD.Branch:
                            node.branch_i(in_t)
                        elif node.type == BD.Demand:
                            node.demand_i(in_t)
                        elif node.type == BD.UpperReservoir or node.type == BD.LowerReservoir:
                            model.reservoirs[typeID].moc_i(in_t)
                            if self.mode =="delta":
                                model.reservoirs[typeID].moc_i(in_t,waterLevel=0)
                        elif node.type == BD.DeadEnd:
                            node.deadend_i(in_t)    
                        elif node.type == BD.EndValve:
                            node.endvalve_i(in_t)  
                        elif node.type == BD.InlineValve:
                            node.inlineValve_i(in_t)  
                        else:
                            a=1
                            print('''Unknown interior boundary''')             
                    elif node.id in self.unknowns:
                        '''Exterior boundary'''
                        '''unknown exterior boundary'''
                        if node.type == BD.Demand:
                            node.demand_i(in_t)
                        else:
                            pass
                            # node.unknown_interior_node_i(in_t)

            for pipe in self.pipes:
                for xx in range(pipe.NN+1):
                    for i in range(pipe.front_t[xx]):
                        if np.isnan(pipe.hi[i,xx]) or np.isnan(pipe.qi[i,xx]) :
                            print("Invalid value in initializating the known zone for pipe " +str(pipe.id))          

    def update_pipe(self,pipe,ti):
        dire=pipe.dire
        if dire == -1:
            '''EMOC from the end node to the start node of the pipe'''
            k = pipe.NN - 1
            while k > -1:
                k2 = k
                t1 = pipe.front_t[k2] +ti
                if t1<1:
                    k -= 1
                    continue
                hA1 = pipe.hi[t1 - 1, k2 + 1]
                qA1 = pipe.qi[t1 - 1, k2 + 1]
                hA2 = pipe.hi[t1 + 1, k2 + 1]
                qA2 = pipe.qi[t1 + 1, k2 + 1]
                bA2 = pipe.bi[t1 +1, k2 + 1]
                bB = pipe.bi[t1, k2]
                pipe.qi[t1, k2] = (-hA1 + (1 + pipe.km2)**2 * hA2 + pipe.m1 * ((1 + pipe.km2)*qA2 + qA1) + (1+pipe.km2) * bA2 +bB) / \
                    ((2 +pipe.km2) * pipe.m1 - pipe.m2 * ((1 +pipe.km2) *qA2 - abs(qA1)))
                if pipe.qi[t1, k2]<0:
                    kkk=pipe.qi[t1, k2]
                    pipe.qi[t1, k2] = (-hA1 + (1 + pipe.km2)**2 * hA2 + pipe.m1 * ((1 + pipe.km2) * qA2 + qA1) + pipe.km2 * bA2 + bB) / \
                                        ((2 + pipe.km2) * pipe.m1 - pipe.m2 * (-(1 + pipe.km2) * qA2 - abs(qA1)))
                    if  pipe.qi[t1, k2]>0:
                        print('error in flow direction -1')
                        # exit(0)  
                pipe.hi[t1, k2] = (hA1 + pipe.m1 * (pipe.qi[t1, k2] - qA1) + pipe.m2 * pipe.qi[t1, k2] * abs(qA1) -bB) /(1 +pipe.km2)
                # go across the last node
                if k2 == 0  and (pipe.js not in self.unknowns): 
                    if pipe.start_node.type == BD.Series:
                        pipe.start_node.cross_series_i(t1,pipe.id)
                    elif pipe.start_node.type == BD.Branch:
                        pipe.start_node.cross_branch_i(t1, pipe.hi[t1, k2])
                    elif pipe.start_node.type == BD.Demand:
                        pipe.start_node.cross_demand_i(t1,dire)
                # if k2==0 and (pipe.js in self.extra_sensor):
                #     delta=pipe.hi[t1, k2] -self.extra_sensor_data[t1]
                #     pipe.hi[t1, k2]=pipe.hi[t1, k2]-delta
                k -= 1
            
        elif dire == 1:
            '''EMOC from the start node to the end node of the pipe'''
            k = 1
            while k < pipe.NN + 1:
                # if i < j:
                #     break
                k2 = k
                t1 = pipe.front_t[k2] +ti
                if t1<1:
                    k += 1
                    continue
                hA1 = pipe.hi[t1 + 1, k2 - 1]
                qA1 = pipe.qi[t1 + 1, k2 - 1]
                hA2 = pipe.hi[t1 - 1, k2 - 1]
                qA2 = pipe.qi[t1 - 1, k2 - 1]
                bA1 = pipe.bi[t1 +1, k2 - 1]
                bB = pipe.bi[t1, k2]
                bB=0
                bA1=0
                # go upwards
                pipe.qi[t1, k2] = (hA2 - (1 + pipe.km2)**2 * hA1 + pipe.m1 * ((1 + pipe.km2) * qA1 + qA2) - (1 + pipe.km2) * bA1 - bB) / \
                    ((2 + pipe.km2) * pipe.m1 - pipe.m2 * ((1 + pipe.km2) * qA1 - abs(qA2)))
                if pipe.qi[t1, k2]<0:
                    pipe.qi[t1, k2] = (hA2 - (1 + pipe.km2)**2 * hA1 + pipe.m1 * ((1 + pipe.km2) * qA1 + qA2) - (1 + pipe.km2) * bA1 - bB) / \
                    ((2 + pipe.km2) * pipe.m1 - pipe.m2 * (-(1 + pipe.km2) * qA1 - abs(qA2)))
                    if  pipe.qi[t1, k2]>0:
                        A=1
                        print('error in flow direction 1')
                        print('ti',ti)
                        # exit(0) 
                pipe.hi[t1, k2] = hA1 + pipe.m1 * (pipe.qi[t1, k2] - qA1) - pipe.m2 * pipe.qi[t1, k2] * abs(qA1)+pipe.km2 * hA1+bA1
                # go across the last node
                if (k2 - pipe.NN == 0) and (pipe.je not in self.unknowns):
                    if pipe.end_node.type == BD.Series:
                        pipe.end_node.cross_series_i(t1,pipe.id)
                    elif pipe.end_node.type == BD.Branch:
                        pipe.end_node.cross_branch_i(t1, pipe.hi[t1, k2])
                    elif pipe.end_node.type == BD.Demand:
                        pipe.end_node.cross_demand_i(t1,dire)
                # if k2==pipe.NN  and (pipe.je in self.extra_sensor):
                #     delta=pipe.hi[t1, k2] -self.extra_sensor_data[t1]
                #     pipe.hi[t1, k2]=pipe.hi[t1, k2]-delta
                k += 1
        
    def min_max_front(self):
        '''Returns the minimum and maximum values of the front'''
        max_front = 0
        min_front = 0
        for k in range(self.nPipes):
            pipe = self.pipes[k]
            if max(pipe.front_t) >max_front:
                max_front = max(pipe.front_t)
            if min(pipe.front_t) <min_front:
                min_front = min(pipe.front_t)
        return min_front,max_front
    
    def front_checking(self):
        '''Check if the branch node and sensor node meet the requirements'''
        flag=True
        for node_route in self.node_routes:
            for j in node_route[:-1]:
                node=self.nodes[j]
                if node.type == BD.Branch:
                    flag=node.check_front()
                    if flag== False:
                        break
        if flag==True:
            print("Front checking succeed!")

    def update_front(self):
        min_level=0
        # update the front in sensor-boundary pipes
        for i in range(self.nsensor):
            for k,pipe_id in enumerate( self.pipe_routes[i]):
                pipe=self.pipes[pipe_id]
                if pipe.dire==1:
                    # update front at the sensor node       
                    if k==0:
                        node=self.nodes[pipe.js] # sensor node
                        if pipe.front_t[0]>node.front_t.min() and node.np>1:
                            # If the sensor node failed to meet the requirements, let the front at sensor node be the (minimum front - 1)
                            pipe.front_t[0]=node.front_t.min()-1
                        # pipe.front_t[-1] -= i
                    # update front of other nodes, decreasing linearly
                    elif k>0:
                        # previous pipe
                        pipe0=self.pipes[self.pipe_routes[i][k-1]]
                        # use the minimum front of the previous pipe to update the front of the next pipe
                        for pipe1 in self.nodes[pipe.js].pipes:
                            if pipe1.js == pipe.js:
                                pipe1.front_t[0] = pipe0.front_t.min()
                            else:
                                pipe1.front_t[-1] = pipe0.front_t.min()
                    # update front of internal points
                    for j in range(pipe.NN):
                        pipe.front_t[j + 1] = pipe.front_t[j]-1
                    node2=pipe.end_node
                else: # from the end node to the start node of the pipe
                    if k==0:
                        node=self.nodes[pipe.je]   
                        # min_front,count=np.unique(node.front_t)
                        if  pipe.front_t[-1]>node.front_t.min(): # BUG:if [min, min, max]
                            pipe.front_t[-1]=node.front_t.min()-1
                        # pipe.front_t[-1]-=i
                    elif k>0:
                        pipe0=self.pipes[self.pipe_routes[i][k-1]]
                        for pipe1 in self.nodes[pipe.je].pipes:
                            if pipe1.je == pipe.je:
                                pipe1.front_t[-1] = pipe0.front_t.min()                   
                            else:
                                pipe1.front_t[0] = pipe0.front_t.min()
                    for j in range(pipe.NN):
                        pipe.front_t[-1-j-1] = pipe.front_t[-1-j]-1
                    node2=pipe.start_node
                self.update_node_front(node)
                self.update_node_front(node2)
                if min_level >node2.front_t.min():
                    min_level = node2.front_t.min()

        # Negative values of fronts exists after updating the front in sensor-boundary pipes
        # Move up the new front in sensor-boundary pipes until the minimum front is 0
        # Let the state under the new front at boundaries nodes be knowns (set them 0) to find the new front for all the pipe
        for i in range(self.nsensor):
            node=self.nodes[self.node_routes[i][-1]]
            self.delta_new_front[i]= node.front_t.min()-min_level #
            if len(self.pipe_routes[i])==0:
                self.delta_new_front[i]=self.model.total_seg
        self.print_front('Modified front')
        self.find_front()
        self.print_front("2nd find front")
        max_f=0
        for i,j in enumerate(self.delta_new_front):
            if len(self.pipe_routes[i])==0:
                continue
            elif j>max_f:
                max_f=j
                
        for pipe in self.pipes:
            pipe.front_t-=(max_f)
        # update the front at nodes
        for node in self.nodes:
            self.update_node_front(node) 
        self.print_front("Maximum font at unknown nodes is 0")
    
    def update_node_front(self,node):
        for i, pipe in enumerate(node.pipes):
                    if pipe.js == node.id:
                        node.front_t[i] = pipe.front_t[0]
                    else:
                        node.front_t[i] = pipe.front_t[-1]
    
    def plot_front(self ):
        plt.figure(10,figsize=(10, 3))
        num=0
        for i,pipe in enumerate(self.pipes):
            
            if i > 4:
                break
            plt.plot(np.linspace(num, num + pipe.NN, pipe.NN + 1, dtype='int'), pipe.front_t, 'k')
            plt.plot(np.linspace(num,num+pipe.NN,pipe.NN+1,dtype='int'),pipe.front_t,'ok',markersize=4)
            num += pipe.NN
            
        # plt.plot(np.linspace(4, 6,3, dtype='int'), self.pipes[5].front_t, 'r')
        # plt.plot(np.linspace(6, 8, 3, dtype='int'), self.pipes[6].front_t, 'r')
        plt.plot(np.linspace(4, 6, 3, dtype='int'), self.pipes[5].front_t, 'or',markersize=4)
        plt.plot(np.linspace(6, 8, 3, dtype='int'), self.pipes[6].front_t, 'or', markersize=4)
        # plt.grid(1)
        # plt.xticks(np.arange(0,20,1))
        # plt.yticks(np.arange(0,13,1))
        # plt.axis([0,20,0,13])
        plt.show()
    
    def find_front(self):
        '''Finds the front'''
        # # # # # # # # # # # # # # # # # #
        '''initiallize qi and hi'''
        for pipe in self.pipes:
            pipe.front_t = np.zeros(pipe.NN+1,dtype='int')
            pipe.hi[1:] = np.nan
            pipe.qi[1:] = np.nan
        
        # set 0 at sensors and known reservoirs, set nan at unknown boundaries
        for i, node_ids in enumerate(self.node_routes):
            # set the head at sensors to 0
            node = self.nodes[node_ids[0]]
            for pipei in node.start_p:
                pipei.hi[1:,-1] = 0
            for pipei in node.end_p:
                pipei.hi[1:,0] = 0
            # set the head at unknown boundaries under the new front t0 0
            node = self.nodes[node_ids[-1]]
            for pipei in node.start_p:
                pipei.hi[:1+self.delta_new_front[i], -1] = 0
            for pipei in node.end_p:
                pipei.hi[:1 + self.delta_new_front[i], 0] = 0
            if node.type == BD.Demand:
                node.demand = np.nan*np.ones(self.steps)
                node.demand[:1 + self.delta_new_front[i]] = 0
                node.obj.cda2g = np.nan*np.ones(self.steps)
                node.obj.cda2g[:1 + self.delta_new_front[i]] = 0
                # node.obj.cda2g = np.nan
        # set the head of known reservoirs to 0
        for node in self.nodes:
            if node.id not in self.unknowns:
                if node.type == BD.LowerReservoir :     
                    node.start_p[0].hi[:,-1]=0
                elif node.type == BD.UpperReservoir:
                    node.end_p[0].hi[:,0]=0
                elif node.type == BD.DeadEnd:
                    if node.n_sp:
                        node.pipes[0].qi[:,-1]=0
                    else:
                        node.pipes[0].qi[:,0]=0
        '''run MOC'''
        model = self.model
        pipes = self.pipes
        self.initialize_known_zone()
        ''' Find nan
            The last point which the flow and head are not Nan is the front
            It represents the last point where the information from known IC/BC can reach'''
        max_front=0
        second_max_front=0
        for ti in range(self.max_init_time+1):
            for pipe in pipes:
                for i in range(pipe.NN +1):
                    if not ((np.isnan(pipe.hi[ti, i]) or np.isnan(pipe.qi[ti,i]))):
                        pipe.front_t[i] = ti  
                        if pipe.front_t[i]>max_front:
                            max_front = pipe.front_t[i]
                        
        for pipe in pipes:
            for i in range(pipe.NN +1):
                if pipe.front_t[i] >second_max_front and pipe.front_t[i]<max_front:
                    second_max_front = pipe.front_t[i]
        for pipe in pipes:
            for i in range(pipe.NN +1):
                if pipe.front_t[i] ==max_front:
                    pipe.front_t[i]=second_max_front+1
        self.max_init_time=second_max_front+2
        for node in model.nodes:
            node.front_t=np.zeros(node.np,dtype="int")
            self.update_node_front(node)
        # self.plot_front()
        # self.print_front()

    def pick_data(self, h, pipnum):
        return self.pipes[pipnum].hi[:, h]
    
    def print_front(self, hint=''):
        print("Front: " + hint+" ###############")
        for pipe in self.pipes:
            print("pipe" + str(pipe.id) + ": ", pipe.front_t)

    def save_results(self,addname=''):
        self.nodesDataH={}
        self.nodesDataQ={}
        for node in self.nodes:
            for pipe in self.pipes:
                if pipe.js==node.id:
                    # plt.plot(pipe.hi[:self.plotSteps, 0], '--', label='Predicted')
                    nodeDataH=pipe.hi[:self.plotSteps, 0]
                    nodeDataQ=pipe.qi[:self.plotSteps, 0]
                    break
                elif pipe.je==node.id:
                    # plt.plot(pipe.hi[:self.plotSteps, -1], '--', label='Predicted')
                    nodeDataH=pipe.hi[:self.plotSteps, -1]
                    nodeDataQ=pipe.qi[:self.plotSteps, -1]
                    break
            self.nodesDataH['node'+str(node.id)]=nodeDataH
            self.nodesDataQ['node'+str(node.id)]=nodeDataQ

        np.savez(data_root_path+'emoc_results_'+self.name+'_'+addname+'_nodes_dataH.npz',**self.nodesDataH)
        np.savez(data_root_path+'emoc_results_'+self.name+'_'+addname+'_nodes_dataQ.npz',**self.nodesDataQ)

    def loadResults(self,filename):
        loaded_H=np.load(data_root_path+filename+"_nodes_dataH.npz")
        for node in self.nodes:
            for pipe in self.pipes:
                array_nameH = f'node{node.id}'
                loaded_arrayH = loaded_H[array_nameH]
                self.plotSteps=loaded_arrayH.shape[0]
                if pipe.js==node.id:
                    # plt.plot(pipe.hi[:self.plotSteps, 0], '--', label='Predicted')
                    pipe.hi[:self.plotSteps, 0]=loaded_arrayH
                    break
                elif pipe.je==node.id:
                    # plt.plot(pipe.hi[:self.plotSteps, -1], '--', label='Predicted')
                    nodeDataH=pipe.hi[:self.plotSteps, -1]=loaded_arrayH
                    break

