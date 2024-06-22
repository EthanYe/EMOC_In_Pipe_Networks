import configparser
import os
import sys
# sys.path.append(sys.path[0] + '\\Hydraulic Structure')
import numpy as np
from module.components import *

import json


class ReadIni:
    def __init__(self,name):
        os.chdir(sys.path[0])  # Change the dir to the path of current .py file
        self.cf = configparser.ConfigParser()  # Creat a config handle
        self.cf.read(name+"_config.ini", encoding="utf-8")  # Read .ini file


    def read_control_papameters(self):
        # Read control papameters
        cf=self.cf
        if cf.has_section('CONTROL PARAMETERS'):
            SysParam.time_step = cf.getfloat('CONTROL PARAMETERS', 'TimeStep')
            SysParam.total_time = cf.getfloat('CONTROL PARAMETERS', 'TotalTimes')
            SysParam.g = cf.getfloat('CONTROL PARAMETERS', 'Gravity')
            SysParam.factor = cf.get('CONTROL PARAMETERS', 'Factor')

    def read_pipes(self):  # Read pipes
        cf = self.cf
        pipes = []
        for i in range(Pipe.MAX_PIPE):
            name = 'L' + str(i + 1)
            if cf.has_section(name):
                js = cf.getint(name, 'J_s')
                je = cf.getint(name, 'J_e')
                wavespeed = cf.getfloat(name, 'wavespeed')
                diameter = cf.getfloat(name, 'diameter')
                length = cf.getfloat(name, 'length')
                f = cf.getfloat(name, 'f')
                pipes.append(Pipe(i, js, je, wavespeed, diameter, length, f))
            else:
                Pipe.number = i  # Get pipe number
                break  # Reading pipes ends
        return pipes

    def read_reservoirs(self):  # Read pipes
        cf = self.cf
        reservoirs = []
        #
        for i in range(Reservoir.MAX_RESERVOIR):
            name = 'R' + str(i + 1)
            if cf.has_section(name):
                node = cf.getint(name, 'J')
                level = cf.getfloat(name, 'level')
                reservoirs.append(Reservoir(i, node, level))
            else:
                Reservoir.number = i  # Get pipe number
                break  # Reading pipes ends
        return reservoirs

    def read_ballvalves(self):  # Read valves
        cf = self.cf
        ballvalves = []
        #
        for i in range(BallValve.MAX_BALLVALVE):
            name = 'BV' + str(i + 1)
            if cf.has_section(name):
                node = cf.getint(name, 'J')
                ballvalves.append(BallValve(i, node))
            else:
                BallValve.number = i  # Get pipe number
                break  # Reading pipes ends
        return ballvalves

    def read_output(self):  # Read valves
        cf = self.cf
        output = []
        #
        for i in range(Pipe.MAX_PIPE):
            name = 'OUTPUT' + str(i + 1)
            if cf.has_section(name):
                pipe_id = cf.getint(name, 'pip')
                output.append(pipe_id)
            else:
                break  # Reading pipes ends
        return output

# def 
class ReadJson:
    def __init__(self,name=None):
        self.name=name
        if name is None:
            try:    
                with open('option.json','r') as f:
                    option=json.load(f)
                if 'option' in option:
                    self.name=dict['option']
            except FileNotFoundError:
                error("option.json is not found.")
        try:
            with open(name + '_config.json', 'r') as f:
                self.js = json.load(f)
        except FileNotFoundError:
            error(name + '_config.json is not found.')
        # self.sort_nodes()
    def read_control_papameters(self):
        js = self.js
        if 'control parameters' in js:
            cp = js['control parameters']
            try:   
                SysParam.time_step = cp['time step']
                SysParam.total_time = cp['total time']
                SysParam.g = cp['gravity']
                SysParam.factor = cp['friction factor']
                SysParam.steps = int(round(SysParam.total_time / SysParam.time_step))
                SysParam.tlen = SysParam.steps
                SysParam.tt = SysParam.time_step
                SysParam.t = np.arange(0, SysParam.total_time, SysParam.time_step)
                SysParam.freqs = fft.fftfreq(SysParam.steps, SysParam.time_step)  # f
                SysParam.wfs = SysParam.freqs * 2 * np.pi # w=2pif
                if 'wn' in cp:
                    SysParam.wn = cp['wn']
                if 'unit' in cp:
                    SysParam.unit = cp['unit']
                else:
                    SysParam.unit = 'SI'
            except KeyError:
                print('Incomplete control parameters in config file!')
        else:
            print('control parameters are not defined!')
    
    def default_dump(self,obj):
        """Convert numpy classes to JSON serializable objects."""
        if isinstance(obj, (np.integer, np.floating, np.bool_)):
            return obj.item()
        elif isinstance(obj, np.ndarray):
            return obj.tolist()
        else:
            return obj
    
    def sort_nodes(self):
        js = self.js
        if 'pipes' not in js:
            error('No pipe is defined!')
        pdict = js['pipes']
        pipes = []
        Pipe.number = len(pdict)
        nodes = []
        for i, p in enumerate(pdict):
            p['id'] = i
            js = p['js']
            je = p['je']
            if js not in nodes:
                nodes.append(js)
            if je not in nodes:
                nodes.append(je)
        nodes_id_sorted = np.sort(np.array(list(set(nodes)), dtype='int'))

        for p in pdict:
            js = p['js']
            je = p['je']
            n, m = np.where(nodes_id_sorted == js)[0][0], np.where(nodes_id_sorted == je)[0][0]
            p['js'], p['je'] = n, m
        for part in self.js:
            try:
                for p in self.js[str(part)]:
                    if 'node' in p:
                        node= p['node']
                        p['node'] = np.where(nodes_id_sorted == node)[0][0]
            except:
                continue
        with open(self.name + '_config.json', 'w') as f:
            json.dump(self.js, f, ensure_ascii=False, default=self.default_dump)
        return pipes, Pipe.number
    
    def read_pipes(self):  # Read pipes
        js = self.js    
        pdict = js['pipes']
        Pipe.number=len(pdict)
        pipes = []
        for p in pdict:
            id = p['id']
            js =p['js'] 
            je = p['je']
            wavespeed = p['wavespeed']
            if SysParam.unit == 'US':
                diameter = p['diameter']/1000*25.4
                length = p['length']/100*30.48
            else :
                diameter = p['diameter']
                length = p['length']
            if 'alpha' in p:
                alpha = p['alpha']
            else:
                alpha = 0
            if 'f' in p:
                f = p['f']
                pipes.append(Pipe(id, js, je, wavespeed, diameter, length, f, alpha=alpha))
            elif 'roughness' in p:
                n = p['roughness']
                pipes.append(Pipe(id, js, je, wavespeed, diameter, length,  alpha=alpha, n=n))

        return pipes, Pipe.number

    def read_reservoirs(self):  # Read pipes
        js = self.js
        reservoirs = []
        Reservoir.number=0
        if 'reservoirs' not in js:
            print('No reservoir is defined!')
            return reservoirs,Reservoir.number
        else:
            rdict = js['reservoirs']

        Reservoir.number=len(rdict)
        #
        for r in rdict:
            id = r['id']
            node = r['node']
            if "water level file" in r:
                water_levels = np.load(r['water level file']+'.npy')
                reservoirs.append(Reservoir(id, node, water_levels=water_levels))
                continue
            if "mode" in r:
                mode = r['mode']
                if mode=='sinusoidal':
                    water_level = r['water level']
                    A = r['A']
                    w = r['w']
                reservoirs.append(Reservoir(id, node, water_level,A=A, w=w))
                continue
            water_level = r['water level']
            reservoirs.append(Reservoir(id,node,water_level))
        return reservoirs,Reservoir.number

    def read_output(self):  # Read valves
        js = self.js
        for i in range(Pipe.MAX_PIPE):
            if 'output' in js:
                return js['output']
            else:
                return []

    def read_ballvalves(self):  # Read valves
        js = self.js
        ballvalves = []
        if 'ballvalves' not in js:
            return ballvalves,0
        bdict = js['ballvalves']
        BallValve.number=len(bdict)
        for b in bdict:
            id=b['id']
            node = b['node']
            ballvalves.append(BallValve(id, node))
        return ballvalves,BallValve.number

    def read_inlineValves(self):  # Read pipes
        js = self.js
        inlineValves = []
        if 'inline valve' not in js:
            return inlineValves, 0
        else:
            inv_dict = js['inline valve']
        
        InlineValve.number = len(inv_dict)
        #
        Cd, A, closingTime, tau_end=0,0,0,0
        for inv in inv_dict:
            id = inv['id']
            node = inv['node']
            if "status" in inv:
                status=inv['status']
                inlineValves.append(InlineValve(id, node, Cd, A, closingTime, tau_end,status=status))
            else:
                inlineValves.append(InlineValve(id, node, Cd, A, closingTime, tau_end))
        return inlineValves, InlineValve.number

    def read_endValves(self):  # Read pipes
        js = self.js
        endValves = []
        if 'end valve' not in js:
            return endValves, 0
        else:
            env_dict = js['end valve']
        
        EndValve.number = len(env_dict)
        #
        for env in env_dict:
            id = env['id']
            node = env['node']
            if 'Q0' in env:
                        Q0 = env['Q0']
            if 'closingTime' in env:
                        closingTime = env['closingTime']
            # if 'motion' in env:
            motion = env['motion']
            if motion=="sinusoidal":
                amplitude = env['amplitude']
                tau0 = env['tau0']
                wf = env['wf'] 
                closingTime = env['closingTime']
                Q0 = env['Q0']
                endValves.append(EndValve(id, node, tau0=tau0, closingTime=closingTime,
                                    amplitude=amplitude, Q0=Q0, motion=motion, wf=wf))
            elif motion == "closed":
                endValves.append(EndValve(id, node, motion=motion))

            elif motion == "static":
                if 'Q0' in env:
                    Q0 = env['Q0']
                    endValves.append(EndValve(id, node, motion=motion, Q0=Q0))
                elif 'iscda' in env:
                    cda=env['cda']
                    endValves.append(EndValve(id, node, motion=motion, cda=cda,iscda=True))
            elif motion == "udf":
                Q0 = env['Q0']
                endValves.append(EndValve(id, node, motion=motion, Q0=Q0))
            elif motion == "linear":
                closingTime = env['closingTime']
               
                Q0 = env['Q0']
                if 'duration' in env:
                    duration=env['duration']
                    endValves.append(EndValve(id, node, closingTime=closingTime, motion=motion,duration=duration, Q0=Q0))
                else:
                    endValves.append(EndValve(id, node, closingTime=closingTime, motion=motion, Q0=Q0))
                # else:
                #     closingTime = env['closingTime']
                #     Q0 = env['Q0']
                    # endValves.append(EndValve(id, node, closingTime=closingTime, motion=motion, Q0=Q0))
            elif motion == "sudden":
                Q0 = env['Q0']
                endValves.append(EndValve(id, node, motion=motion, Q0=Q0))
            else:
                endValves.append(EndValve(id, node, closingTime=closingTime, Q0=Q0))
        return endValves, EndValve.number

    def read_leaks(self):  # Read leaks
        js = self.js
        leaks = []
        if 'leaks' not in js:
            return leaks, 0
        ldict = js['leaks']
        Leak.number = len(ldict)
        for l in ldict:
            id = l['id']
            node = l['node']
            pipeID = l['pipe']
            cda = l['cda']
            burstTime = l['burstTime']
            duration = l['duration']
            elevation = l['elevation']
            leaks.append(Leak(id, pipeID, node, cda, burstTime, duration,elevation))
        return leaks, Leak.number

    def read_nonreflectings(self):  # Read valves
        js = self.js
        nonreflectings = []
        if 'non-reflectings' not in js:
            return nonreflectings, 0
        nfdict = js['non-reflectings']
        Nonreflecting.number = len(nfdict)
        for nf in nfdict:
            id = nf['id']
            node = nf['node']
            nonreflectings.append(Nonreflecting(id, node))
        return nonreflectings, Nonreflecting.number

    def read_demands(self):  # Read leaks
        js = self.js
        demands = []
        if 'demands' not in js:
            return demands, 0
        ldict = js['demands']
        Demand.number = len(ldict)
        for l in ldict:
            id = l['id']
            node = l['node']
            demand = 0
            burstTime=0
            duration=0
            elevation=0
            udf=False
            cda=0
            mode=None
            if 'demand' in l:
                demand = l['demand']
            if 'burstTime' in l:
                burstTime = l['burstTime']
            if 'duration' in l:
                duration = l['duration']
            if 'elevation' in l:
                elevation = l['elevation']
            if "udf" in l:
                udf=l['udf']
            if 'cda' in l:
                cda=l['cda']
            if 'mode' in l:
                mode=l['mode']
            demands.append(Demand(id,  node, demand=demand,burstTime=burstTime, duration=duration,z=elevation,udf=udf,cda=cda,mode=mode))
        return demands, Demand.number

    def read_coords(self):  # Read leaks
        js = self.js
        coords = []
        if 'coordinates' not in js:
            return coords, 0
        cdict = js['coordinates']
        Coordinate.number = len(cdict)
        for c in cdict:
            node = c['node']
            X_coord= c['X-coord']
            Y_coord= c['Y-coord']
            epanode=c["epa_node"]
            coords.append(Coordinate( node,X_coord,Y_coord,epa_node=epanode))
        return coords, Coordinate.number
