import numpy as np
from enum import Enum, unique
import os
import numpy.fft as fft
import matplotlib.pyplot as plt


@unique  # @unique装饰器可以帮助我们检查保证没有重复值
class BD(Enum):
    """ This is the boundary set """
    Series = 0  # Sun的value被设定为0
    UpperReservoir = 1
    LowerReservoir = 2
    Branch = 3
    EndValve = 4
    EndBallValve = 5
    InlineValve = 6,
    DeadEnd = 7,
    NonReflecting = 8,
    Demand = 9


class BV(Enum):
    """ This is the boundary set """
    Start = 0  # Sun的value被设定为0
    Medium = 1
    End = 2


class RE(Enum):
    """ This is the boundary set """
    Upper = 0  # Sun的value被设定为0
    Lower = 1


class POS(Enum):
    Upper = 0
    Internal = 1
    Lower = 2


class SysParam:
    # default system parameters
    # time_step   # TimeStep
    # total_time
    pi = 3.14159  # Circular constant
    g = 9.81  # Gravity acceleration


# This is PIpe


class Pipe(SysParam):  # Pipe model
    # Cylinder
    # Uniform loss factor
    MAX_PIPE = 100
    number = 0

    def __init__(self, id=0, start_node=0, end_node=0, a=0, D=0, length=0, f=0, sigam=0, n=0, z=0, e=0.063, alpha=0, rou=998.1):
        self.id = id  # ID
        self.js = start_node  # Prenode
        self.je = end_node  # Rear node
        self.a = a  # Wave speed
        self.D = D  # Diameter
        self.length = length  # Pipe length
        self.f = f
        self.n = n
        self.time_step = SysParam.time_step
        self.sigam = sigam  # Local water loss factor
        self.dt = SysParam.time_step
        self.zeta = 0
        self.z = z  # Elevation
        # self.cdaP = np.zeros(self.NN + 1)  # Leak of previous time step
        self.leaks = []
        self.nLeak = 0
        # Internal variables
        self.count = 0
        # self.partial_ep = np.zeros(self.NN + 1)
        # self.tau = [1.05, 0.5, 1.5, 5]
        # self.J = [1.048e-5, 1.029e-5, 1.134e-5, 8.083e-7]
        self.alpha = alpha
        self.e = e
        self.rou = rou
        self.get_intermedia_variables()

    def get_intermedia_variables(self):
        # if SysParam.factor == 'n':
        #     self.f = 8 * self.g * self.n * self.n * (self.D / 4)**(-1.0 / 3.0)
        self.A = SysParam.pi * self.D * self.D / 4  # Area
        self.NN = round(self.length / self.time_step / self.a)  # Get segments
        if self.NN == 0:
            self.NN = 1  # reaches>=1

        # self.a = self.length / self.time_step / self.NN  # Adjuct wave speed

        self.dx = self.length / self.NN  # Space step
        self.m1 = self.a / self.g / self.A  # a/gA
        self.m2 = self.time_step * self.f * self.a / 2.0 / self.g / self.D / self.A / self.A  # fadt/(2gDA2)
        self.m3 = self.g * self.A * self.time_step / self.dx  # gAdt/dx
        self.m4 = self.time_step * self.f / 2 / self.D / self.A  # dtf/(2DA)
        self.m5 = self.a ** 2 * self.time_step / self.g / self.A / self.dx  # a^2dt/gAdx
        self.m6 = self.g * self.A * self.time_step   # gAdt
        self.m7 = self.a ** 2 * self.time_step / self.g / self.A  # a^2dt/gA
        self.m8 = 1 / 2.0 / self.g / self.A / self.A  # 1/2gA^2
        # FD3
        self.direction = 'left'
        self.dz = np.ones(self.NN) * self.dx  # Leak of current time step
        self.cda = np.zeros(self.NN + 1)  # Leak of current time step
        # Hydraulic factor
        self.H0 = np.zeros(self.NN + 1)  # Initial head
        self.Q0 = np.zeros(self.NN + 1)  # Initial flow
        self.K0 = np.zeros(self.NN + 1)  # Initial leak
        self.H = np.zeros(self.NN + 1)  # Head of previous time step
        self.Q = np.zeros(self.NN + 1)  # Flow of previous time step
        self.HP = np.zeros(self.NN + 1)  # Head of current time step
        self.QP = np.zeros(self.NN + 1)  # Flow of current time step
        self.HE = np.zeros(self.NN + 1)  # Estimated head of current time step
        self.QE = np.zeros(self.NN + 1)  # Estimated flow of current time step
        self.K = np.zeros(self.NN + 1)  # Leak of previous time step
        self.KP = np.zeros(self.NN + 1)  # Leak of current time step
        self.cda = np.zeros(self.NN + 1)  # Leak of previous time step

        # KVmodel
        # KVmodel damping is larger when tau is smaller or J is larger
        self.tau = [1.05]
        self.J = [4.029e-5]
        self.NK = len(self.tau)
        self.bk = np.zeros((self.NK, self.NN + 1))
        self.b = np.zeros(self.NN + 1)
        self.epsik = np.zeros((self.NK, self.NN + 1))
        self.epsipk = np.zeros((self.NK, self.NN + 1))
        self.epsip = np.zeros(self.NN + 1)
        self.epsi = np.zeros(self.NN + 1)
        # create matrices of q and h
        self.qi = np.zeros((self.steps, self.NN + 1))
        self.hi = np.zeros((self.steps, self.NN + 1))
        self.bi = np.zeros((self.steps, self.NN + 1))
        self.epsipki = np.zeros((self.steps, self.NK, self.NN + 1))
        
        tau = self.tau
        J = self.J
        self.k2k = np.zeros(self.NK)
        self.km1 = np.zeros(self.NK)
        self.km2k = np.zeros(self.NK)
        self.km3k = np.zeros(self.NK)
        self.mm1k = np.zeros(self.NK)
        self.mm2k = np.zeros(self.NK)
        self.mm3k = np.zeros(self.NK)
        self.mm4k = np.zeros(self.NK)
        self.k1 = 2 * self.a * self.a * self.dt / self.g
        self.k2 = 0
        self.k3 = 0
        self.km2 = 0
        self.km3 = 0
        self.mm1 = 0
        self.mm2 = 0
        self.mm3 = 0
        self.mm4 = 0
        # self.k3 = self.rou*SysParam.g* self.alpha /1000* self.D / 2 / self.e
        if self.alpha == 0:
            self.update_b = self.update_b0
            return
        else:
            self.update_b = self.update_bi
        self.k3 = self.alpha * self.D / 2 / self.e
        for i in range(self.NK):
            self.k2k[i] = pow(np.e, -self.dt / tau[i])
            self.km1[i] = self.k1 * self.k2k[i] * self.k3 * J[i] / tau[i]
            self.km2k[i] = self.k1 * J[i] * (1 - self.k2k[i]) * self.k3 / self.dt
            self.km3k[i] = -self.k1 * self.k2k[i] / tau[i]
            self.mm1k[i] = (J[i] * self.k3 - tau[i] / self.k1 *
                            self.km2k[i])
            self.mm2k[i] = self.km1[i] - self.km2k[i]
            self.mm3k[i] = -self.tau[i] / self.k1 * self.mm2k[i]
            self.mm4k[i] = -self.tau[i] / self.k1 * self.km3k[i]
            self.k2 = self.k2 + self.k2k[i]
            self.km2 = self.km2 + self.km2k[i]
            self.km3 = self.km3 + self.km3k[i]
            self.mm1 = self.mm1 + self.mm1k[i]
            self.mm2 = self.mm2 + self.mm2k[i]
            self.mm3 = self.mm3 + self.mm3k[i]
            self.mm4 = self.mm4 + self.mm4k[i]

    def moc(self, T=0, step=None):
        if step is None:
            step = int(T / self.time_step) - 1

        for i, leak in enumerate(self.leaks):
            if T > leak.burstTime - 0.00000001 and T < leak.burstTime + 0.00000001 + leak.duration:
                self.cda[leak.node] = leak.cda * (T - leak.burstTime) / leak.duration
            elif T > leak.burstTime + 0.00000001 + leak.duration:
                self.cda[leak.node] = leak.cda
            if leak.z > self.HP[leak.node]:
                self.KP[leak.node] = 0
            else:
                self.KP[leak.node] = abs(self.cda[leak.node]) * np.sqrt((self.HP[leak.node] - leak.z))
            self.leak_flow[i][step] = self.KP[leak.node]
            # print(self.leak_flow[i][step])
        for j in range(1, self.NN):
          
            self.HP[j] = (self.QCP(j) - self.QCM(j)) / (self.CQP(j) + self.CQM(j))
            # print(self.HP[j])
            self.QP[j] = self.QCP(j) - self.CQP(j) * self.HP[j]
        if len(self.leaks) == 0:
            for j in range(1, self.NN):
                self.KP[j] = self.K[j]
                if self.KP[j] < 0:
                    self.KP[j] = 0

    def update_b0(self, ti, j):
        return

    def update_bi(self, ti, j):

        self.epsi[j] = 0
        self.bi[ti, j] = 0
        if ti == 2 and j == 1 and self.js == 0:
            a = 1
        for i in range(self.NK):  # Update b
            self.bk[i, j] = self.km1[i] * (self.hi[ti - 1, j] - self.hi[0, j]) - self.km2k[i] * \
                self.hi[ti - 1, j] + self.km3k[i] * self.epsipki[ti - 1, i, j]
            self.bi[ti, j] = self.bi[ti, j] + self.bk[i, j]
            if self.js == 0 and j == 1:
                a = 1
            self.epsipki[ti, i, j] = self.mm1k[i] * self.hi[ti, j] - self.J[i] * \
                self.k3 * self.hi[0, j] - self.tau[i] / self.k1 * self.bk[i, j]

    def QCP(self, J):
        qcp = ((self.Q[J - 1] - self.K[J - 1]) * self.m1 + self.H[J - 1] - self.b[J]) / \
            (self.m1 + self.m2 * abs(self.Q[J - 1] - self.K[J - 1]))
        return qcp

    def CQP(self, J):
        cqp = (1 + self.km2) / (self.m1 + self.m2 * abs(self.Q[J - 1] - self.K[J - 1]))
        return cqp

    def QCM(self, J):
        qcm = ((self.Q[J + 1] + self.KP[J]) * self.m1 + self.m2 * self.KP[J] * self.Q[J + 1] -
               self.H[J + 1] + self.b[J]) / (self.m1 + self.m2 * abs(self.Q[J + 1]))
        return qcm

    def CQM(self, J):
        cqm = (1 + self.km2) / (self.m1 + self.m2 * abs(self.Q[J + 1]))
        return cqm

    def QCPi(self, i, j, k=1):
        QA = self.qi[i - 1, j - 1]
        HA = self.hi[i - 1, j - 1]
        if np.isnan(QA ) or np.isnan(HA ):
            a=1
        qcp = ((QA) * self.m1 + HA - self.bi[i, j]) / (self.m1 + k * self.m2 * abs(QA))
        return qcp

    def CQPi(self, i, j, k=1):
        QA = self.qi[i - 1, j - 1]
        if np.isnan(QA ) :
            a=1
        cqp = (1 + self.km2) / (self.m1 + k * self.m2 * abs(QA))
        return cqp

    def QCMi(self, i, j, k=1):
        QB = self.qi[i - 1, j + 1]
        HB = self.hi[i - 1, j + 1]
        if np.isnan(QB ) or np.isnan(HB ):
            a=1
        qcm = ((QB) * self.m1 - HB + self.bi[i, j]) / (self.m1 + k * self.m2 * abs(QB))
        return qcm

    def CQMi(self, i, j, k=1):
        QB = self.qi[i - 1, j + 1]
        if np.isnan(QB ) :
            a=1
        cqm = (1 + self.km2) / (self.m1 + k * self.m2 * abs(QB))
        return cqm

    def FD_lax(self, j):
        '''
        return HP,QP
        '''
        HA = self.H[j - 1]
        HB = self.H[j + 1]
        QA = self.Q[j - 1]
        QB = self.Q[j + 1]
        Qave = (QA + QB) / 2
        QP = Qave - self.m3 * (HB - HA) / 2 - self.m4 * Qave * abs(Qave)
        HP = (HA + HB) / 2 - self.m5 * (QB - QA) / 2
        return HP, QP

    def FD_continuousTime1(self, j):
        '''
        Get HE,QE
        '''
        HC = self.H[j]
        QC = self.Q[j]
        NN = self.NN
        if j == 0:
            self.HE[0] = self.H[0]
            self.QE[0] = QC - self.m6 / self.dz[0] * (self.H[1] - HC) - self.m4 * QC * abs(QC)
            return
        elif j == self.NN:
            self.HE[NN] = self.H[NN]
            # self.QE[NN] = QC - self.m6 / self.dz[NN-1] * (self.H[NN-1] - HC) - self.m4 * QC * abs(QC)
            return
        else:
            HA = self.H[j - 1]
            QA = self.Q[j - 1]
            HB = self.H[j + 1]
            QB = self.Q[j + 1]
            self.HE[j] = HC - self.m7 / self.dz[j] * (QC - QA)
            self.QE[j] = QC - self.m6 / self.dz[j] * (HB - HC) - self.m4 * QC * abs(QC)

    def FD_continuousTime2(self, j):
        '''
        Get HP,QP
        '''
        HC = self.H[j]
        QC = self.Q[j]
        HEC = self.HE[j]
        QEC = self.QE[j]
        NN = self.NN
        if j == 0:
            self.QP[0] = QC - 0.5 * (self.m6 / self.dz[0] * (self.H[1] - HC + self.HE[1] -
                                                             HEC) + self.m4 * QC * abs(QC) + self.m4 * QEC * abs(QEC))
            self.HP[0] = self.H[0]
            return
        elif j == self.NN:
            self.HP[NN] = self.H[NN]
            # self.QP[j] = QC - 0.5 * (self.m6 / self.dz[NN - 1] * (self.H[NN - 1] - HC +
            #                                                 self.HE[NN - 1] - HEC) + self.m4 * QC * abs(QC) + self.m4 * QEC * abs(QEC))
            return
        else:
            HA = self.H[j - 1]
            QA = self.Q[j - 1]
            HB = self.H[j + 1]
            QB = self.Q[j + 1]
            HEA = self.HE[j - 1]
            QEA = self.QE[j - 1]
            HEB = self.HE[j + 1]
            QEB = self.QE[j + 1]
            self.HP[j] = HC - 0.5 * self.m7 / self.dz[j] * (QC - QA + QEC - QEA)
            self.QP[j] = QC - 0.5 * (self.m6 / self.dz[j] * (HB - HC + HEB -
                                                             HEC) + self.m4 * QC * abs(QC) + self.m4 * QEC * abs(QEC))

    def FD_mac1(self, j):
        '''
        get HE, QE
        '''
        HC = self.H[j]
        QC = self.Q[j]
        if self.direction == 'right':  # 1->n
            HB = self.H[j + 1]
            QB = self.Q[j + 1]
            self.HE[j] = HC - self.m5 * (QB - QC)
            self.QE[j] = QC - self.m3 * (HB - HC) - self.m4 * QC * abs(QC)
        elif self.direction == 'left':  # 2->n+1
            HA = self.H[j - 1]
            QA = self.Q[j - 1]
            self.HE[j] = HC - self.m5 * (QC - QA)
            self.QE[j] = QC - self.m3 * (HC - HA) - self.m4 * QC * abs(QC)
        else:
            print('Unknown dire in Fa_mac:', self.direction)
        if j == self.NN:
            self.HE[j] = HC

    def FD_mac2(self, j):
        '''
        return HP,QP
        '''
        HC = self.H[j]
        QC = self.Q[j]
        HEC = self.HE[j]
        QEC = self.QE[j]
        if self.direction == 'right':  # 2->n+1
            HEA = self.HE[j - 1]
            QEA = self.QE[j - 1]
            self.HP[j] = (HC + HEC - self.m5 * (QEC - QEA)) * 0.5
            self.QP[j] = (QC + QEC - self.m3 * (HEC - HEA) - self.m4 * QEC * abs(QEC)) * 0.5
        elif self.direction == 'left':  # 1->n
            HEB = self.HE[j + 1]
            QEB = self.QE[j + 1]
            self.HP[j] = (HC + HEC - self.m5 * (QEB - QEC)) * 0.5
            self.QP[j] = (QC + QEC - self.m3 * (HEB - HEC) - self.m4 * QEC * abs(QEC)) * 0.5
        else:
            print('Unknown dire in Fa_mac:', self.direction)

    def finite_diff(self, format='lax'):
        '''
        Lax(in-node) or max
        '''
        if format == 'lax':
            for j in range(1, self.NN):
                self.HP[j], self.QP[j] = self.FD_lax(j)
        elif format == 'mac':
            if self.direction == 'right':
                for j in range(0, self.NN):
                    self.FD_mac1(j)
                for j in range(1, self.NN):
                    self.FD_mac2(j)
                # self.direction = 'left'
            elif self.direction == 'left':
                for j in range(1, self.NN + 1):
                    self.FD_mac1(j)
                for j in range(0, self.NN):
                    self.FD_mac2(j)
                # self.direction = 'right'  # change the direction
        elif format == 'Henu':
            for j in range(0, self.NN + 1):
                self.FD_continuousTime1(j)
            for j in range(1, self.NN):
                self.FD_continuousTime2(j)
        else:
            print('Unknown format in finite_diff:', format)

    def get_field_matrix(self, omega):
        self.R = (2 * self.f * self.Q0[0]) / (2 * self.g * self.D * self.A**2)
        self.u2 = complex(-omega**2 / self.a**2, self.g * self.A * omega * self.R / self.a**2)
        self.u = self.u2**(0.5)
        self.Zc = complex(0, -self.u * self.a**2 / omega / self.g / self.A)
        self.field_matrix = np.zeros((2, 2), dtype=complex)
        self.field_matrix[1, 1] = self.field_matrix[0, 0] = np.cosh(self.u * self.length)
        temp = np.sinh(self.u * self.length)
        self.field_matrix[0, 1] = -temp / self.Zc
        self.field_matrix[1, 0] = -temp * self.Zc
        # self.field_matrix_inv=np.linalg.inv(self.field_matrix)
        # print(self.field_matrix_inv@self.field_matrix)
        # self.point_matrix=np.eye(2)
        # self.point_matrix[0,0]=0
        # self.field_matrix=self.field_matrix@self.point_matrix
        a = 1

    def addLeak(self, leak):
        self.leaks.append(leak)
        self.nLeak += 1

class Node:
    MAX_NODE = 100
    number = 0

    def __init__(self, id=0, type='', typeID=0):
        self.id = id
        self.type = BD.Series
        self.typeID = typeID
        self.n_sp = 0
        self.n_ep = 0
        self.np = 0
        self.pID = []
        self.spID = []
        self.epID = []
        self.start_p = []
        self.end_p = []
        self.pipes = []
        self.demand=self.demand0 = 0
        self.is_transfer_q = False

    def deadend_moc(self):
        if self.n_sp == 0:
            ep = self.end_p[0]
            NN = ep.NN
            ep.QP[0] = 0
            ep.HP[0] = -ep.QCM(0) / ep.CQM(0)
        elif self.n_ep == 0:
            sp = self.start_p[0]
            NN = sp.NN
            sp.QP[NN] = 0
            sp.HP[NN] = sp.QCP(NN) / sp.CQP(NN)

    def series_moc(self):
        if self.n_sp == 1 and self.n_ep == 1:
            sp = self.start_p[0]
            ep = self.end_p[0]
            # a1 = sp.m8 - ep.m8
            NN = sp.NN
            # b0 = 1 / sp.CQP(NN) + 1 / ep.CQM(0)
            # c0 = -sp.QCP(NN) / sp.CQP(NN) - ep.QCM(0) / ep.CQM(0)
            # a0 = 0
            # if (abs(a0) < 0.000001):
            #     sp.QP[NN] = -c0 / b0
            # else:
            #     y1 = (-b0 - np.sqrt(b0 * b0 - 4.0 * a0 * c0)) / (2.0 * a0)
            #     y2 = (-b0 + np.sqrt(b0 * b0 - 4.0 * a0 * c0)) / (2.0 * a0)
            #     if abs(y1 - sp.Q[NN]) < abs(y2 - sp.Q[NN]):
            #         sp.QP[NN] = y1
            #     else:
            #         sp.QP[NN] = y2
                        # ep.QP[0] = sp.QP[NN]
            # # print(ep.HP[0])
            # sp.HP[NN] = (sp.QCP(NN) - sp.QP[NN]) / sp.CQP(NN)
            # ep.HP[0] = (ep.QP[0] - ep.QCM(0)) / ep.CQM(0)
            sp.HP[ NN] = (sp.QCP( NN) - ep.QCM( 0)) / (sp.CQP( NN)+ep.CQM( 0))
            ep.HP[ 0] = sp.HP[ NN]
            ep.QP[ 0] = sp.QP[ NN]=sp.QCP(NN) - sp.CQP(NN) * sp.HP[NN]

        elif self.n_sp == 2:
            sp1, sp2 = self.start_p[0], self.start_p[1]
            NN1, NN2 = sp1.NN, sp2.NN
            sp1.HP[NN1] = sp2.HP[NN2] = (sp1.QCP(NN1) + sp2.QCP(NN2)) / (sp1.CQP(NN1) + sp2.CQP(NN2))
            sp1.QP[NN1] = sp1.QCP(NN1) - sp1.CQP(NN1) * sp1.HP[NN1]
            sp2.QP[NN2] = -sp1.QP[NN1]
        elif self.n_ep == 2:
            ep1, ep2 = self.end_p[0], self.end_p[1]
            NN1, NN2 = 0, 0
            ep1.HP[NN1] = ep2.HP[NN2] = (-ep1.QCM(NN1) - ep2.QCM(NN2)) / (ep1.CQM(NN1) + ep2.CQM(NN2))
            ep1.QP[NN1] = ep1.QCM(NN1) + ep1.CQM(NN1) * ep1.HP[NN1]
            ep2.QP[NN2] = -ep1.QP[NN1]

    def series_i(self, i):  # for RMOC
        if self.n_sp == 1 and self.n_ep == 1:
            sp = self.start_p[0]
            ep = self.end_p[0]
            NN = sp.NN
            sp.hi[i, NN] = (sp.QCPi(i, NN) - ep.QCMi(i, 0)) / (sp.CQPi(i, NN)+ep.CQMi(i, 0))
            ep.hi[i, 0] = sp.hi[i, NN]
            ep.qi[i, 0] = sp.qi[i, NN]=sp.QCPi(i,NN) - sp.CQPi(i,NN) * sp.hi[i,NN]
        elif self.n_sp == 2:
            sp1, sp2 = self.start_p[0], self.start_p[1]
            NN1, NN2 = sp1.NN, sp2.NN
            sp1.hi[i, NN1] = sp2.hi[i, NN2] = (sp1.QCPi(i, NN1) + sp2.QCPi(i, NN2)) / \
                (sp1.CQPi(i, NN1) + sp2.CQPi(i, NN2))
            sp1.qi[i, NN1] = sp1.QCPi(i, NN1) - sp1.CQPi(i, NN1) * sp1.hi[i, NN1]
            sp2.qi[i, NN2] = -sp1.qi[i, NN1]
        elif self.n_ep == 2:
            ep1, ep2 = self.end_p[0], self.end_p[1]
            NN1, NN2 = 0, 0
            ep1.hi[i, NN1] = ep2.hi[i, NN2] = (-ep1.QCMi(i, NN1) - ep2.QCMi(i, NN2)) / \
                (ep1.CQMi(i, NN1) + ep2.CQMi(i, NN2))
            ep1.qi[i, NN1] = ep1.QCMi(i, NN1) + ep1.CQMi(i, NN1) * ep1.hi[i, NN1]
            ep2.qi[i, NN2] = -ep1.qi[i, NN1]
    
    def deadend_i(self,i):
        if self.n_sp == 0:
            ep = self.end_p[0]
            NN = ep.NN
            ep.qi[i, 0] = 0
            ep.hi[i, 0] = - ep.QCMi(i, 0) / ep.CQMi(i, 0)
        elif self.n_ep == 0:
            sp = self.start_p[0]
            NN = sp.NN
            sp.qi[i, NN] = 0
            sp.hi[i, NN] = sp.QCPi(i, NN) / sp.CQPi(i, NN)
    
    def endvalve_i(self,i):
        if self.n_sp == 0:
            ep = self.end_p[0]
            NN = ep.NN
            ep.qi[i, 0] = self.obj.Q0
            ep.hi[i, 0] = - ep.QCMi(i, 0) / ep.CQMi(i, 0)
        elif self.n_ep == 0:
            sp = self.start_p[0]
            NN = sp.NN
            sp.qi[i, NN] =  self.obj.Q0
            sp.hi[i, NN] = sp.QCPi(i, NN) / sp.CQPi(i, NN)

    def sensor_i(self, i, sensor_h=0):
        for sp in self.start_p:
            sp.hi[i, sp.NN] = sensor_h
            sp.qi[i, sp.NN] = sp.QCPi(i, sp.NN) - sp.CQPi(i, sp.NN) * sp.hi[i, sp.NN]
        for ep in self.end_p:
            ep.hi[i, 0] = sensor_h
            ep.qi[i, 0] = ep.QCMi(i, 0) + ep.CQMi(i, 0) * ep.hi[i, 0]

    def cross_series_i(self, i,pipe_id):  # for RMOC
        

        if self.n_sp ==2 :
            sp1=self.start_p[0]
            sp2=self.start_p[1]
            if sp1.id==pipe_id:
                sp2.qi[i, sp2.NN]=-sp1.qi[i, sp1.NN]
                sp2.hi[i, sp2.NN]= sp1.hi[i, sp1.NN]
            else:
                sp1.qi[i, sp1.NN]=-sp2.qi[i, sp2.NN]
                sp1.hi[i, sp1.NN]=sp2.hi[i, sp2.NN]
        elif self.n_ep ==2:
            ep1=self.end_p[0]
            ep2=self.end_p[1]
            if ep1.id==pipe_id:
                ep2.qi[i, 0]=-ep1.qi[i, 0]
                ep2.hi[i, 0]=ep1.hi[i, 0]
            else:
                ep1.qi[i, 0] = -ep2.qi[i, 0]
                ep1.hi[i, 0] = ep2.hi[i, 0]
        else:
            sp = self.start_p[0]
            ep = self.end_p[0]
            NN = sp.NN
            if sp.id==pipe_id:
                ep.qi[i, 0] = sp.qi[i, NN]
                # sp.hi[i, NN] = (sp.QCPi(i, NN) - sp.qi[i, NN]) / sp.CQPi(i, NN)
                ep.hi[i, 0] = sp.hi[i, NN]
            else:
                sp.qi[i, NN]= ep.qi[i, 0]
                sp.hi[i, NN]=ep.hi[i, 0] 
    
    def cross_branch_i(self, ti, local_h):  # for RMOC
        
        other_pipe_q = 0
        # get h and q in start pipes (exclude forward pipe)
        for ep in self.end_p:
            if ep.id==self.nextPipe.id:
                continue
            ep.hi[ti, 0] = local_h
            ep.qi[ti, 0] = ep.QCMi(ti, 0) + ep.CQMi(ti, 0) * ep.hi[ti, 0]
            other_pipe_q += ep.qi[ti, 0]
        # get h and q in end pipes (exclude forward pipe)
        for sp in self.start_p:
            if sp.id==self.nextPipe.id:
                continue
            sp.hi[ti, sp.NN] = local_h
            sp.qi[ti, sp.NN] = sp.QCPi(ti, sp.NN) - sp.CQPi(ti, sp.NN) * sp.hi[ti, sp.NN]
            other_pipe_q -= sp.qi[ti, sp.NN]

        # get h and q in forward pipes
        pipe = self.nextPipe
        if pipe.id in self.spID:
            pipe.hi[ti, pipe.NN] = local_h
            pipe.qi[ti, pipe.NN] = other_pipe_q
            if np.isnan(pipe.qi[ti, pipe.NN] ):
                a=1
        else:
            pipe.hi[ti, 0] = local_h
            pipe.qi[ti, 0] = -other_pipe_q
            if np.isnan(pipe.qi[ti, 0] ):
                a=1
        return True
    
    def check_front(self):
        ''' check sovability'''
        # The second front in other pipe should be greater than the the sencond ront in the next SB pipe 
        if self.nextPipe.dire==1:
            min_front=self.nextPipe.front_t[1]
        else:
            min_front=self.nextPipe.front_t[-2]
        for pipe in self.pipes:
            NN=pipe.NN
            if pipe.id== self.nextPipe.id:
                continue
            if pipe.id in self.spID:
                front =pipe.front_t[NN-1]
            else:
                front =pipe.front_t[1]
            if front<=min_front:
                print('Warning: Node '+str(self.id)+' is not solvable')
                return False
        return True

    """ old functions for EMOC2"""
    '''
    # def cross_branch_i(self, ti, local_h, pipe0):  # for RMOC
    #     pipe = self.forward_pipe
    #     if pipe.id in self.spID:
    #         front_1 = pipe.front_t[pipe.NN]
    #         front_2 = pipe.front_t[pipe.NN - 1]
    #     else:
    #         front_1 = pipe.front_t[0]
    #         front_2 = pipe.front_t[1]
    #     if front_1 < front_2:
    #         if pipe.id in self.spID:
    #             pipe.hi[ti, pipe.NN] = local_h
    #             pipe.qi[ti, pipe.NN] = pipe.qi[ti, pipe.NN] = pipe.QCPi(
    #                 ti, pipe.NN) - pipe.CQPi(ti, pipe.NN) * pipe.hi[ti, pipe.NN]
    #         else:
    #             pipe.hi[ti, 0] = local_h
    #             pipe.qi[ti, 0] = pipe.QCMi(ti, 0) + pipe.CQMi(ti, 0) * pipe.hi[ti, 0]
    #     else:
    #         # get h and q in start pipes (exclude forward pipe)
    #         for ep in self.other_ep:
    #             ep.hi[ti, 0] = local_h
    #             ep.qi[ti, 0] = ep.QCMi(ti, 0) + ep.CQMi(ti, 0) * ep.hi[ti, 0]
    #         # get h and q in end pipes (exclude forward pipe)
    #         for sp in self.other_sp:
    #             sp.hi[ti, sp.NN] = local_h
    #             sp.qi[ti, sp.NN] = sp.QCPi(ti, sp.NN) - sp.CQPi(ti, sp.NN) * sp.hi[ti, sp.NN]
    #         # get q in other pipes
    #         other_pipe_q = 0
    #         for osp in self.other_sp:
    #             other_pipe_q -= osp.qi[ti, osp.NN]
    #         for oep in self.other_ep:
    #             other_pipe_q += oep.qi[ti, 0]
    #         # get h and q in forward pipes
    #         pipe = self.forward_pipe
    #         if pipe.id in self.spID:
    #             pipe.hi[ti, pipe.NN] = local_h
    #             pipe.qi[ti, pipe.NN] = other_pipe_q
    #         else:
    #             pipe.hi[ti, 0] = local_h
    #             pipe.qi[ti, 0] = -other_pipe_q
    '''
    
    def branch_moc(self):
        l = 0.0
        m = 0.0
        B = 0.0
        for sp in self.start_p:
            l = l + sp.QCP(sp.NN)
            B = B + sp.CQP(sp.NN)
        for ep in self.end_p:
            m = m - ep.QCM(0)
            B = B + ep.CQM(0)
        if self.demand!=0:
            print("Demand is not 0 at branch node")
        hp = (l + m-self.demand) / B
        for sp in self.start_p:
            sp.HP[sp.NN] = hp
            sp.QP[sp.NN] = sp.QCP(sp.NN) - sp.CQP(sp.NN) * sp.HP[sp.NN]
        for ep in self.end_p:
            ep.HP[0] = hp
            ep.QP[0] = ep.QCM(0) + ep.CQM(0) * ep.HP[0]

    def branch_i(self, i):
        l = 0.0
        m = 0.0
        B = 0.0
        for sp in self.start_p:
            l = l + sp.QCPi(i, sp.NN)
            B = B + sp.CQPi(i, sp.NN)
        for ep in self.end_p:
            m = m - ep.QCMi(i, 0)
            B = B + ep.CQMi(i, 0)
        hp = (l + m) / B
        for sp in self.start_p:
            sp.hi[i, sp.NN] = hp
            sp.qi[i, sp.NN] = sp.QCPi(i, sp.NN) - sp.CQPi(i, sp.NN) * sp.hi[i, sp.NN]
            if np.isnan( sp.qi[i, sp.NN]):
                a=1
        for ep in self.end_p:
            ep.hi[i, 0] = hp
            ep.qi[i, 0] = ep.QCMi(i, 0) + ep.CQMi(i, 0) * ep.hi[i, 0]
            if np.isnan(ep.qi[i, 0]):
                a = 1

    def add_pipe(self, position, pipe):
        if position == 'end pipe':
            self.n_ep = self.n_ep + 1  # Update the number of ending pipe
            self.epID.append(pipe.id)  # Add the pipe id in the ending pipe list
            self.end_p.append(pipe)  # Add the pipe in the ending pipe list
        elif position == 'start pipe':
            self.n_sp = self.n_sp + 1  # Update the number of ending pipe
            self.spID.append(pipe.id)  # Add the pipe id in the ending pipe list
            self.start_p.append(pipe)  # Add the pipe id in the ending pipe list
        else:
            print('Pipe position is wrong!')

    def demand_i(self, i):
        sp = self.start_p[0]
        ep = self.end_p[0]
        NN = sp.NN
        # if knownPipe == 'sp':
        #     ep.hi[i, 0] = sp.hi[i, NN]
        #     ep.qi[i, 0] = ep.QCMi(i, 0) + ep.CQMi(i, 0) * ep.hi[i, 0]
        # else:
        #     sp.hi[i, NN]=ep.hi[i, 0]
        #     sp.qi[i, NN] = sp.QCPi(i, NN) - sp.CQPi(i, NN) * sp.hi[i, NN]
        # self.demandi=self.demand
        # if sp.hi[i-1,NN]-self.obj.z<0:
        #     self.demandi=0
        # else:
        #     self.demandi = self.obj.cda2g*np.sqrt(sp.hi[i-1,NN]-self.obj.z)
        self.demandi = self.obj.cda2g*np.sqrt(sp.hi[i-1,NN]-self.obj.z)
        sp.hi[i,NN] = (sp.QCPi(i,NN) - ep.QCMi(i,0) - self.demandi[i]) / (sp.CQPi(i,NN) + ep.CQMi(i,0))
        sp.qi[i,NN] = sp.QCPi(i,NN) - sp.CQPi(i,NN) * sp.hi[i,NN]
        ep.qi[i,0] = sp.qi[i,NN] - self.demandi[i]
        ep.hi[i,0] = sp.hi[i,NN]
    
    def unknown_demand_i(self, i):
        sp = self.start_p[0]
        ep = self.end_p[0]
        NN = sp.NN
        try:
            if sp.dire:
                ep.hi[i, 0] = sp.hi[i, NN]
                ep.qi[i, 0] = ep.QCMi(i, 0) + ep.CQMi(i, 0) * ep.hi[i, 0]
        except:
                sp.hi[i, NN] = ep.hi[i, 0]
                sp.qi[i, NN] = sp.QCPi(i, NN) - sp.CQPi(i, NN) * sp.hi[i, NN]
        a=1
    
    def cross_demand_i(self,i,dire):
        sp = self.start_p[0]
        ep = self.end_p[0]
        NN = sp.NN
            
        if dire == 1:
            demand = self.obj.cda2g * np.sqrt(sp.hi[i - 1, NN] - self.obj.z)
            ep.hi[i, 0] = sp.hi[i, NN]
            ep.qi[i, 0] = sp.qi[i, NN]- demand
            a=1
        else:
            print('Warning: undefined process')
    
    def inlineValve_i(self,i):
        if self.obj.status=='closed':
            if self.n_sp == 1 and self.n_ep == 1:
                sp = self.start_p[0]
                ep = self.end_p[0]
                NN = sp.NN
                sp.qi[i,NN]=0
                ep.qi[i,0]=0
                ep.hi[i,0] = (ep.qi[i,0] - ep.QCMi(i,0)) / ep.CQMi(i,0)
                sp.hi[i,NN] = (sp.QCPi(i,NN) - sp.qi[i,NN]) / sp.CQPi(i,NN)
            elif self.n_sp == 2 and self.n_ep == 0:
                sp1, sp2 = self.start_p[0], self.start_p[1]
                NN1, NN2 = sp1.NN, sp2.NN
                sp2.qi[i,NN2] = sp1.qi[i,NN1]=0
                sp1.hi[i,NN1] = (sp1.QCPi(i,NN1) - sp1.qi[i,NN1]) / sp1.CQPi(i,NN1)
                sp2.hi[i,NN2] = (sp2.QCPi(i,NN2) - sp2.qi[i,NN2]) / sp2.CQPi(i,NN2)
            elif self.n_ep == 2 and self.n_sp == 0:
                ep1, ep2 = self.end_p[0], self.end_p[1]
                NN1, NN2 = 0, 0
                ep2.qi[i,NN2] = ep1.qi[i,NN1]=0
                ep1.hi[i,NN1] =(ep1.qi[i,0] - ep1.QCMi(i,0)) / ep1.CQMi(i,0)
                ep2.hi[i,NN2] = (ep2.qi[i,0] - ep2.QCMi(i,0)) / ep2.CQMi(i,0)

        elif self.obj.status=="open":
            if self.n_sp == 1 and self.n_ep == 1:
                sp = self.start_p[0]
                ep = self.end_p[0]
                NN = sp.NN
                sp.hi[i,NN] = (sp.QCPi(i, NN) - ep.QCMi(i, 0)) / (sp.CQPi( i,NN)+ep.CQMi(i, 0))
                ep.hi[i, 0] = sp.hi[i, NN]
                ep.qi[ i,0] = sp.qi[ i,NN]=sp.QCPi(i,NN) - sp.CQPi(i,NN) * sp.hi[i,NN]

            elif self.n_sp == 2:
                sp1, sp2 = self.start_p[0], self.start_p[1]
                NN1, NN2 = sp1.NN, sp2.NN
                sp1.hi[i,NN1] = sp2.hi[i,NN2] = (sp1.QCPi(i,NN1) + sp2.QCPi(i,NN2)) / (sp1.CQPi(i,NN1) + sp2.CQPi(i,NN2))
                sp1.qi[i,NN1] = sp1.QCPi(i,NN1) - sp1.CQPi(i,NN1) * sp1.hi[i,NN1]
                sp2.qi[i,NN2] = -sp1.qi[i,NN1]
            elif self.n_ep == 2:
                ep1, ep2 = self.end_p[0], self.end_p[1]
                NN1, NN2 = 0, 0
                ep1.hi[i,NN1] = ep2.hi[i,NN2] = (-ep1.QCMi(i,NN1) - ep2.QCMi(i,NN2)) / (ep1.CQMi(i,NN1) + ep2.CQMi(i,NN2))
                ep1.qi[i,NN1] = ep1.QCMi(i,NN1) + ep1.CQMi(i,NN1) * ep1.hi[i,NN1]
                ep2.qi[i,NN2] = -ep1.qi[i,NN1]

    def unknown_interior_node_i(self, i):
        for pipe in self.end_p:
            pipe.qi[i, 0] = pipe.QCMi(i, 0) + pipe.CQMi(i, 0) * pipe.hi[i, 0]
        for pipe in self.start_p:
            NN = pipe.NN
            pipe.qi[i, NN] = pipe.QCPi(i, NN) - pipe.CQPi(i, NN) * pipe.hi[i, NN]
    

    def get_point_matrix(self, omega):
        self.point_matrix = np.eye(2)
        if self.type == BD.DeadEnd:
            self.point_matrix[0, 0] = 0
        elif self.type == BD.Branch:
            self.point_matrix[0, 1] = 1
        elif self.type == BD.EndValve:
            self.point_matrix = np.eye(3)
            if self.id == self.pipes[0].js:
                k = 0
            else:
                k = self.pipes[0].NN
            self.point_matrix[1, 0] = -2.0 * self.pipes[0].H0[k] / self.pipes[0].Q0[k]
            self.point_matrix[1, 2] = 2.0 * self.pipes[0].H0[k] * self.obj.amplitude
    # def get_fQ_dfQ(self,Q):
    #     if self.type==BD.InlineValve:
    #         self.fQ=
    #     return self.fQ,self.dfQ

    def get_other_pipe(self, pipe):
        other_pipe = self.pipes.copy()
        other_pipe.remove(pipe)
        if len(other_pipe) == 1:
            return other_pipe[0]
        return other_pipe

class Leak():
    MAX_LEAK = 100
    number = 0

    def __init__(self, id, pipeID, node, cda, burstTime=0, duration=0, z=0):
        self.id = id
        self.pipeID = pipeID
        self.node = node
        self.cda = cda
        self.burstTime = burstTime
        self.duration = duration
        self.z = z

class Reservoir(SysParam):
    MAX_RESERVOIR = 100
    number = 0

    def __init__(self, id: 'identity', node, water_level=0, water_levels=None, A=0, w=0):
        self.id = id
        self.node_id = node
        self.water_levels = water_levels
        if water_levels is not None:
            self.water_level = water_levels[0]
            self.water_levels = water_levels
        else:
            self.water_levels=water_level*np.ones(self.steps)
            self.water_level = water_level
        self.A = A
        self.w = w

    def connect(self, nodes):
        # Recognize resevoirs
        self.node = node = nodes[self.node_id]
        self.pipes = []
        self.end_p = []
        self.start_p = []
        if node.n_sp == 0 and node.n_ep == 1:
            self.type = RE.Upper
            self.end_p += [node.end_p[0]]
            node.type = BD.UpperReservoir
            self.pipes += [node.end_p[0]]
        elif node.n_ep == 0 and node.n_sp == 1:
            self.type = RE.Lower
            self.start_p += [node.start_p[0]]
            node.type = BD.LowerReservoir
            self.pipes += [node.start_p[0]]
        else:
            print('Reservoir definition error! Please check it!')
            exit(0)
        node.typeID = self.id
        node.obj=self
    # def input_water_level(self):
        
    def moc_i(self, i,waterLevel=None):
        if waterLevel is None:
            waterLevel = self.water_level
        if self.type == RE.Upper:
            pipe = self.end_p[0]
            pipe.hi[i, 0] = waterLevel
            pipe.qi[i, 0] = pipe.QCMi(i, 0) + pipe.CQMi(i, 0) * pipe.hi[i, 0]
        elif self.type == RE.Lower:
            pipe = self.start_p[0]
            NN = pipe.NN
            pipe.hi[i, NN] = waterLevel
            pipe.qi[i, NN] = pipe.QCPi(i, NN) - pipe.CQPi(i, NN) * pipe.hi[i, NN]
        else:
            error('undefined reservoir!')
    

    def get_water_levle(self,step=0,time=0):
        if step == 0:
            step = int(round(time / self.time_step))
        if step == self.steps:
            step -= 1
        self.water_level = self.water_levels[step]
        return self.water_levels[step]

    def moc(self, T):
        self.water_level = self.get_water_levle(time=T)
        # if T<0.000001:
        #     self.water_level = 34.79527
        # lower=int(steps/10)
        # frac=steps/10.0-lower
        # self.water_level = self.water_levels[lower] + (self.water_levels[lower]-self.water_levels[lower+1])*frac
        if self.type == RE.Upper:
            pipe = self.end_p[0]
            # /*pipe.HP[0] = waterLevel
            # pipe.QP[0] = pipe.QCM(0) + pipe.CQM(0) * pipe.HP[0]
            # a2 =  pipe.m8 * pipe.CQM(0)
            # pipe.QP[0] = (np.sqrt(1 + 4 * a2 * (pipe.QCM(0) + pipe.CQM(0)*self.water_level)) - 1) / 2 / a2
            # pipe.HP[0] = self.water_level - abs(pipe.QP[0]) * pipe.QP[0] * pipe.m8
            end_p = self.end_p[0]
            self.noise = np.random.randn(1) * 2
            end_p.KP[0] = 0
            end_p.HP[0] = self.water_level * 1.0 + self.A * np.sin(self.w * np.pi * T) + 0 * self.noise
            end_p.QP[0] = end_p.QCM(0) + end_p.CQM(0) * end_p.HP[0]

        elif self.type == RE.Lower:
            pipe = self.start_p[0]
            NN = pipe.NN
            pipe.HP[NN] = self.water_level * 1.0 + self.A * np.sin(self.w * np.pi * T)
            pipe.QP[NN] = pipe.QCP(NN) - pipe.CQP(NN) * pipe.HP[NN]

        # if self.type == RE.Upper:  # Upstream reservoir
        #     end_p = self.end_p  # Get end pipe object
        #     end_p.KP[0] = 0
        #     end_p.HP[0] = self.water_level
        #     end_p.QP[0] = end_p.QCM(0) + end_p.CQM(0) * end_p.HP[0]
        # elif self.type == RE.Lower:  # Downstream reservoir
        #     start_p = self.start_p  # Get start pipe object
        #     NN = start_p.NN  # Get segments
        #     start_p.KP[0] = 0
        #     start_p.HP[NN] = self.water_level
        #     start_p.QP[NN] = start_p.QCP(NN) - start_p.CQP(NN) * start_p.HP[NN]

    def solve(self, format):
        if format == 'moc' or format == 'lax':
            if self.type == RE.Upper:  # Upstream reservoir
                end_p = self.end_p  # Get end pipe object
                end_p.KP[0] = 0
                end_p.HP[0] = self.water_level
                end_p.QP[0] = end_p.QCM(0) + end_p.CQM(0) * end_p.HP[0]
            elif self.type == RE.Lower:  # Downstream reservoir
                start_p = self.start_p  # Get start pipe object
                NN = start_p.NN  # Get segments
                start_p.KP[0] = 0
                start_p.HP[NN] = self.water_level
                start_p.QP[NN] = start_p.QCP(NN) - start_p.CQP(NN) * start_p.HP[NN]
        elif format == 'Henu':
            # pass

            if self.type == RE.Upper:  # Upstream reservoir
                end_p = self.end_p  # Get end pipe object
                end_p.KP[0] = 0
                end_p.HP[0] = self.water_level
                end_p.QP[0] = end_p.QCM(0) + end_p.CQM(0) * end_p.HP[0]
            elif self.type == RE.Lower:  # Downstream reservoir
                start_p = self.start_p  # Get start pipe object
                NN = start_p.NN  # Get segments
                start_p.KP[0] = 0
                start_p.HP[NN] = self.water_level
                start_p.QP[NN] = start_p.QCP(NN) - start_p.CQP(NN) * start_p.HP[NN]

            # if self.type == RE.Upper:  # Upstream reservoir
            #     end_p = self.end_p  # Get end pipe object
            #     end_p.KP[0] = 0
            #     end_p.HP[0] = self.water_level
            #     end_p.QP[0] =end_p.Q[0]  - 0.5 * (end_p.m6 / end_p.dz[0] * (end_p.H[1] - end_p.H[0]+ end_p.HE[1] -
            #                                                                 end_p.HE[0]) + end_p.m4 * end_p.Q[0] * abs(end_p.Q[0]) + end_p.m4 * end_p.QE[0] * abs(end_p.QE[0]))
            # elif self.type == RE.Lower:  # Downstream reservoir
            #     start_p = self.start_p  # Get start pipe object
            #     NN = start_p.NN  # Get segments
            #     start_p.KP[0] = 0
            #     start_p.HP[NN] = self.water_level
            #     start_p.QP[NN] = start_p.Q[NN] - 0.5 * (-start_p.m6 / start_p.dz[NN-1] * (start_p.H[NN-1] - start_p.H[NN] + start_p.HE[NN-1] -
            #                                                        start_p.HE[NN]) + start_p.m4 * start_p.Q[NN] * abs(start_p.Q[NN]) + start_p.m4 * start_p.QE[NN] * abs(start_p.QE[NN]))
        elif format == 'mac':
            if self.type == RE.Upper:  # Upstream reservoir
                end_p = self.end_p  # Get end pipe object
                end_p.KP[0] = 0
                end_p.HP[0] = self.water_level
                # end_p.QP[0] = end_p.QCM(0) + end_p.CQM(0) * end_p.HP[0]
            elif self.type == RE.Lower:  # Downstream reservoir
                start_p = self.start_p  # Get start pipe object
                NN = start_p.NN  # Get segments
                start_p.KP[0] = 0
                start_p.HP[NN] = self.water_level
                # start_p.QP[NN] = start_p.QCP(NN) - start_p.CQP(NN) * start_p.HP[NN]

    def finite_diff(self):
        if self.type == RE.Upper:  # Upstream reservoir
            end_p = self.end_p  # Get end pipe object
            end_p.HP[0] = self.water_level
            end_p.QP[0] = end_p.QCM(0) + end_p.CQM(0) * end_p.HP[0]
        elif self.type == RE.Lower:  # Downstream reservoir
            start_p = self.start_p  # Get start pipe object
            NN = start_p.NN  # Get segments

            start_p.HP[NN] = self.water_level
            start_p.QP[NN] = start_p.QCP(NN) - start_p.CQP(NN) * start_p.HP[NN]

class EndValve(SysParam):
    '''
    Sudden close
    This valve has two status flow = Q0 or flow = 0
    '''
    MAX_ENDVALVE = 100
    number = 0

    def __init__(self, id: 'identity', node_id, tau0=0, amplitude=0,duration=0,z=0,cda=0, closingTime=0, Q0=0, motion='sudden', wf=0, cda_end=0, iscda=False):
        self.id = id
        self.node_id = node_id
        self.closing_time = closingTime
        self.Q0 = Q0
        self.motion = motion
        self.cda_end = cda_end
        
        self.amplitude = amplitude
        self.wf = wf
        self.duration=duration
        self.tv = np.linspace(0, self.total_time, self.steps)
        self.taus = np.zeros(self.steps)
        self.z=z
        if motion == 'sudden':
            self.moc = self.moc_sudden
        elif motion == 'sinusoidal':
            self.moc = self.moc_sinusoidal
            for i in range(self.steps):
                T = i * self.time_step
                self.taus[i] = self.amplitude * np.sin(self.wf * 2 * T)
        elif motion == 'linear':
            self.moc = self.moc_linear
        elif motion == 'closed':
            self.moc = self.moc_closed
        elif motion == 'static':
            self.iscda = iscda
            if iscda:
                self.cda0=cda
            self.moc = self.moc_static
        elif motion == 'udf':
            self.moc = self.moc_udf
        else:
            error("Undefined valve motion")
            # for i in range(self.steps):
            #     T=i*self.time_step
            #     if T > closingTime:  # T > ct, Q=0
            #         self.taus[i]=0
            #     else:
            #         self.taus[i] = (1 - T / closingTime)

        self.complex_f = fft.fft(self.taus)
        self.amplitudes = np.abs(self.complex_f)
        self.real_f = self.complex_f.real
        # self.amplitudes=np.abs(self.complex_f)/self.complex_f.shape[0]
        # self.real_f = self.complex_f.real / self.complex_f.shape[0]
        a = 1

    def connect(self, nodes):
        node = nodes[self.node_id]
        node.typeID = self.id
        node.type = BD.EndValve
        self.n_sp=self.n_ep=0
        if node.n_sp:
            self.start_p = node.pipes[0]
            self.n_sp=node.n_sp
        else:
            self.end_p = node.pipes[0]
            self.n_ep=node.n_ep
        node.demand = self.Q0
        node.obj = self

    def moc_sudden(self, T):   # have not consider the start valve
        # update Q
        if T > self.closing_time:  # T > ct, Q=0
            self.moc_Q(0)
        else:
            self.moc_Q(self.Q0)

    def moc_closed(self, T):   # have not consider the start valve
        self.moc_Q(0)

    def moc_static(self, T):   # have not consider the start valve
        if not self.iscda:
            self.moc_Q(self.Q0)
        else:
            self.moc_cda(self.cda0)

    def moc_sinusoidal(self, T):   # have not consider the start valve

        # tau=self.amplitude*np.exp(complex(0,self.wf*2*T))
        # tau=tau.real+1
        tau = 1. + self.amplitude * np.sin(self.wf * 2 * T)
        cda = tau * self.cda0
        self.moc_cda(cda)

    def moc_udf(self, T):
        T=T
        cda = (1. + 0.2*(0.2 * np.sin(T * 8) - 0.3 * np.cos(4 * T + np.pi / 2) + 0.2 * np.sin(T * 5)**2 - T * 0.05)) * self.cda0
        self.moc_cda(cda)

    def moc_linear(self, T):
        if T<self.closing_time:
            cda=self.cda0
        elif (T >= self.closing_time) and (T < self.closing_time+self.duration):  # T > ct, Q=0
            cda = self.cda0 - (T-self.closing_time) / self.duration * (self.cda0 - self.cda_end)     
        else:
            cda = self.cda_end

        # cda = 0
        # cda =(1. + 0.2 * np.sin(T * 8) - 0.3 * np.cos(4 * T + np.pi / 2) + 0.2 * np.sin(T * 5)**2 - T * 0.05) * self.cda0
        self.moc_cda(cda)
        # cda = (1. + 0.2 * np.sin(T * 0.8) - 0.3 * np.cos(0.4 * T + np.pi / 2) + 0.2 * np.sin(T * 0.5)**2 - T * 0.005) * self.cda0
        # if T < self.closing_time:  # T > ct, Q=0
        #     cda=self.cda0-0.2*T/self.closing_time*(self.cda0-self.cda_end)
        # elif  T < 2*self.closing_time:zz
        #     cda=self.cda0*0.8-0.8*(T-self.closing_time)/self.closing_time*(self.cda0-self.cda_end)

        # else:
        #     cda=self.cda_end
        # if T < self.closing_time:  # T > ct, Q=0
        #     Q = self.Q0 * (1 - T / self.closing_time)
        # else:
        #     Q=0
        # self.moc_Q(Q)

    def moc_Q(self, Q):   # have not consider the start valve
        if self.n_sp:
            start_p = self.start_p
            # get segments
            NN = start_p.NN
            # update Q
            start_p.QP[NN] = Q  # T <= ct, Q=Q0
            # update H
            start_p.HP[NN] = start_p.QCP(NN) / start_p.CQP(NN) - \
                1.0 / start_p.CQP(NN) * start_p.QP[NN]
            # update KP
            start_p.KP[NN] = start_p.K[0]
        else:
            end_p = self.end_p
            # get segments
            
            # update Q
            end_p.QP[0] = Q  # T <= ct, Q=Q0
            # update H
            end_p.HP[0] = -end_p.QCM(0) / end_p.CQM(0) + \
                1.0 / end_p.CQM(0) * end_p.QP[0]
            # update KP
            end_p.KP[0] = end_p.K[0]
            

    def moc_cda(self, cda):
        if cda == 0:
            self.moc_Q(0)
            return
        if self.n_sp:
            sp = self.start_p
            # get segments
            NN = sp.NN
            if sp.Q[NN] > 0:
                A = 1.0 / 2.0 / self.g / cda / cda
            else:
                sp.QP[NN] =0
                sp.HP[NN] = sp.QCP(NN) / sp.CQP(NN) 
                # update KP
                sp.KP[NN] = sp.K[0]
                # error("Q<0 at end valve")
                return
            b = 1 / sp.CQP(NN)
            c = -sp.QCP(NN) / sp.CQP(NN)
            delta = b * b - 4 * A * c
            if delta < 0:
                A = -A
                delta = -delta
            y1 = (-b + np.sqrt(delta)) / 2 / A
            y2 = (-b - np.sqrt(delta)) / 2 / A
            if np.abs(y1 - sp.Q[NN]) < np.abs(y2 - sp.Q[NN]):
                sp.QP[NN] = y1
            else:
                sp.QP[NN] = y2
            # update H
            sp.HP[NN] = sp.QCP(NN) / sp.CQP(NN) - \
                1.0 / sp.CQP(NN) * sp.QP[NN]
            # update KP
            sp.KP[NN] = sp.K[0]
        else:
            ep = self.end_p
            # get segments
            if ep.Q[0] < 0:
                A = 1.0 / 2.0 / self.g / cda / cda
            else:
                ep.QP[0] =0
                ep.HP[0] = -ep.QCM(0) / ep.CQM(0) 
                # update KP
                ep.KP[0] = ep.K[0]
                # error("Q<0 at end valve")
                return
            b = -1 / ep.CQM(0)
            c = ep.QCM(0) / ep.CQM(0)
            delta = b * b - 4 * A * c
            if delta < 0:
                A = -A
                delta = -delta
            y1 = (-b + np.sqrt(delta)) / 2 / A
            y2 = (-b - np.sqrt(delta)) / 2 / A
            if np.abs(y1 - ep.Q[0]) < np.abs(y2 - ep.Q[0]):
                ep.QP[0] = y1
            else:
                ep.QP[0] = y2
            # update H
            ep.HP[0] = -ep.QCM(0) / ep.CQM(0) + \
                1.0 / ep.CQM(0) * ep.QP[0]
            # update KP
            ep.KP[0] = ep.K[0]

class BallValve():
    # The maximum of valves
    MAX_BALLVALVE = 100
    # The existing number
    number = 0
    number_end = 0

    def __init__(self, typeID: 'identity', node_id, name='bv'):
        self.id = typeID
        # The node where the valve exists
        self.node_id = node_id
        # Set prefix name
        self.__name = name
        self.read_param()

    def read_param(self):
        # Get filename: prefix name + id
        fname1 = self.__name + str(self.id + 1) + "_ang_CDA.txt"
        # read Ang_CDA file, eg: bv1_ang_CDA
        self.ang_CDA = np.loadtxt(fname1, usecols=(1, 2))
        # read Time_Ang file, eg: bv1_time_ang
        fname2 = self.__name + str(self.id + 1) + "_time_ang.txt"
        self.t_ang = np.loadtxt(fname2, usecols=(1, 2))

    def connect(self, nodes):
        node = nodes[self.node_id]
        node.typeID = self.id

        if node.n_ep == 0 and node.n_sp == 1:
            self.type = BV.End
            self.start_p = node.start_p[0]
            node.type = BD.EndBallValve
            BallValve.number_end = BallValve.number_end + 1
        elif node.n_ep == 1 and node.n_sp == 1:
            self.type = BV.Medium
            self.start_p = node.start_p[0]
            self.end_p = node.end_p[0]
        else:
            print('Reservoir definition error! Please check it!')
            exit(0)

    def moc(self, T):
        if self.type == BV.End:
            start_p = self.start_p
            # get segments
            NN = start_p.NN
            # get CDA by T
            CDA = self.T2CDA(T)
            if CDA > 0.0001:
                # update Q
                a = start_p.CQP(NN) / CDA / CDA / 2 / SysParam.g
                b = 1
                c = -start_p.QCP(NN)
                delta = b * b - 4 * a * c
                if delta < 0:
                    print('delta<0')
                    print('a=', a, ', b=', b, ', c=', c)
                    print('Time=', T)
                    exit(0)

                start_p.QP[NN] = (-b + np.sqrt(b * b - 4 * a * c)) / 2 / a
            else:
                # Ballvalve full closing
                start_p.QP[NN] = 0
            # update H
            # update H
            start_p.HP[NN] = (start_p.QCP(NN) - start_p.QP[NN]) / start_p.CQP(NN)
            # update KP
            start_p.KP[NN] = start_p.K[0]
        else:
            print("Unkown ballvalve type!")

    def T2CDA(self, T):
        # interpolate ang by T
        ang = np.interp(T, self.t_ang[:, 0], self.t_ang[:, 1])
        CDA = np.interp(ang, self.ang_CDA[:, 0], self.ang_CDA[:, 1])
        return CDA

    def dCDA_dt(self, T):
        return (self.T2CDA(T) - self.T2CDA(T - self.Time_step)) / self.Time_step

    def dCDA(self, T):
        return (self.T2CDA(T) - self.T2CDA(T - self.Time_step))

class InlineValve(SysParam):
    # The maximum of valves
    MAX_BALLVALVE = 100
    # The existing number
    number = 0
    number_end = 0

    def __init__(self, typeID: 'identity', node_id, Cd, A, closingTime, tau_end, name='bv',status=None):
        self.id = typeID
        # The node where the valve exists
        self.node_id = node_id
        self.Cd = Cd
        self.A = A
        self.closingTime = closingTime
        self.tau_end = tau_end
        # Set prefix name
        self.__name = name
        self.cda0 = Cd * A
        
        self.status=status
        if status is not None:
            self.moc=self.moc_status
            if status=='open':
                self.s = 0
            else:
                self.s=1e16
        else:
            self.moc=self.moc_cda

    def connect(self, nodes):
        node = nodes[self.node_id]
        node.typeID = self.id
        node.type = BD.InlineValve
        self.start_p = node.start_p
        self.end_p = node.end_p
        self.n_sp=node.n_sp
        self.n_ep=node.n_ep
        node.obj=self

    def moc_status(self,T):
        if self.status=='closed':
            if self.n_sp == 1 and self.n_ep == 1:
                sp = self.start_p[0]
                ep = self.end_p[0]
                NN = sp.NN
                sp.QP[NN]=0
                ep.QP[0]=0
                ep.HP[0] = (ep.QP[0] - ep.QCM(0)) / ep.CQM(0)
                sp.HP[NN] = (sp.QCP(NN) - sp.QP[NN]) / sp.CQP(NN)
            elif self.n_sp == 2 and self.n_ep == 0:
                sp1, sp2 = self.start_p[0], self.start_p[1]
                NN1, NN2 = sp1.NN, sp2.NN
                sp2.QP[NN2] = sp1.QP[NN1]=0
                sp1.HP[NN1] = (sp1.QCP(NN1) - sp1.QP[NN1]) / sp1.CQP(NN1)
                sp2.HP[NN2] = (sp2.QCP(NN2) - sp2.QP[NN2]) / sp2.CQP(NN2)
            elif self.n_ep == 2 and self.n_sp == 0:
                ep1, ep2 = self.end_p[0], self.end_p[1]
                NN1, NN2 = 0, 0
                ep2.QP[NN2] = ep1.QP[NN1]=0
                ep1.HP[NN1] =(ep1.QP[0] - ep1.QCM(0)) / ep1.CQM(0)
                ep2.HP[NN2] = (ep2.QP[0] - ep2.QCM(0)) / ep2.CQM(0)

        elif self.status=="open":
            if self.n_sp == 1 and self.n_ep == 1:
                sp = self.start_p[0]
                ep = self.end_p[0]
                NN = sp.NN
                sp.HP[ NN] = (sp.QCP( NN) - ep.QCM( 0)) / (sp.CQP( NN)+ep.CQM( 0))
                ep.HP[ 0] = sp.HP[ NN]
                ep.QP[ 0] = sp.QP[ NN]=sp.QCP(NN) - sp.CQP(NN) * sp.HP[NN]

            elif self.n_sp == 2:
                sp1, sp2 = self.start_p[0], self.start_p[1]
                NN1, NN2 = sp1.NN, sp2.NN
                sp1.HP[NN1] = sp2.HP[NN2] = (sp1.QCP(NN1) + sp2.QCP(NN2)) / (sp1.CQP(NN1) + sp2.CQP(NN2))
                sp1.QP[NN1] = sp1.QCP(NN1) - sp1.CQP(NN1) * sp1.HP[NN1]
                sp2.QP[NN2] = -sp1.QP[NN1]
            elif self.n_ep == 2:
                ep1, ep2 = self.end_p[0], self.end_p[1]
                NN1, NN2 = 0, 0
                ep1.HP[NN1] = ep2.HP[NN2] = (-ep1.QCM(NN1) - ep2.QCM(NN2)) / (ep1.CQM(NN1) + ep2.CQM(NN2))
                ep1.QP[NN1] = ep1.QCM(NN1) + ep1.CQM(NN1) * ep1.HP[NN1]
                ep2.QP[NN2] = -ep1.QP[NN1]
        else:
            error("undefined inline valve")

    def moc_cda(self, T):
        sp = self.start_p
        ep = self.end_p
        NN = sp.NN
        if T < self.closingTime:
            tau = 1 - T / self.closingTime * (1 - self.tau_end)
        else:
            tau = self.tau_end
        cda = self.cda0 * tau
        if sp.Q[NN] > 0:
            A = 1.0 / 2.0 / self.g / cda / cda
        else:
            A = -1.0 / 2.0 / self.g / cda / cda

        b = 1 / sp.CQP(NN) + 1 / ep.CQM(0)
        c = -sp.QCP(NN) / sp.CQP(NN) - ep.QCM(0) / ep.CQM(0)
        delta = b * b - 4 * A * c
        if delta < 0:
            A = -A
            delta = -delta
        y1 = (-b + np.sqrt(delta)) / 2 / A
        y2 = (-b - np.sqrt(delta)) / 2 / A
        if np.abs(y1 - sp.Q[NN]) < np.abs(y2 - sp.Q[NN]):
            sp.QP[NN] = y1
        else:
            sp.QP[NN] = y2
        sp.QP[NN] = (-b + np.sqrt(delta)) / 2 / A
        sp.HP[NN] = (sp.QCP(NN) - sp.QP[NN]) / sp.CQP(NN)
        ep.QP[0] = sp.QP[NN]
        ep.HP[0] = (ep.QP[0] - ep.QCM(0)) / ep.CQM(0)

class Nonreflecting(SysParam):
    MAX_NONREFLECTING = 100
    number = 0

    def __init__(self, id: 'identity', node):
        self.id = id
        self.node_id = node

    def connect(self, nodes):
        # Recognize Nonreflecting
        node = nodes[self.node_id]
        # if node.n_sp == 0 and node.n_ep == 1:
        #     self.end_p = node.end_p[0]
        #     # node.type = BD.UpperReservoir
        #     pass
        if node.n_ep == 0 and node.n_sp == 1:
            self.start_p = node.start_p[0]
            node.type = BD.NonReflecting
        else:
            print('Nonreflecting definition error! Please check it!')
            exit(0)
        node.typeID = self.id

    def moc(self):
        start_p = self.start_p
        NN = start_p.NN
        start_p .QP[NN] = start_p.Q[NN - 1]
        start_p .HP[NN] = (start_p.QCP(NN) - start_p.QP[NN]) / start_p.CQP(NN)

class Demand(SysParam):
    MAX_LEAK = 100
    number = 0

    def __init__(self, id, node, demand=0, cda=None, burstTime=0, duration=0, z=0, udf=False,mode=None):
        self.id = id
        self.node_id = node
        self.burstTime = burstTime
        self.duration = duration
        self.duration2 = self.duration * 1.2
        self.z = z

        self.mode=mode
        if cda is not None:
            self.cda2g0=cda*np.sqrt(2*self.g)
            self.cda2gt0=0
        # elif demand:
        #     self.cda=demand/50
        self.demands = []
        self.udf=udf
        # if self.udf:
        #     self.demand = demand

    def connect(self, nodes):
        node = nodes[self.node_id]
        node.typeID = self.id
        node.type = BD.Demand
        self.start_p = node.start_p
        self.end_p = node.end_p
        node.obj=self
        self.node=node
    def update(self, nodes,node,cda):
        self.node_id=node
        self.cda2g0=cda*np.sqrt(2*self.g)
        if self.node.np<3:
            self.node.type=BD.Series
        else: self.node.type=BD.Branch
        self.node.obj=None
        self.connect(nodes)
        self.demands=[]

    def get_demand(self,T):
        sp = self.node.pipes[0]
        NN = sp.NN
        if self.udf:
            if T> self.burstTime+self.duration:
                dt=(T-self.burstTime-self.duration)*1.6
                cda = self.cda * \
                    (1. - 0.5 * (0.2 * np.sin(dt * 8) -
                                                       0.3 * np.cos(4 * dt + np.pi / 2) + 0.2 * np.sin(dt * 5)**2 - dt * 0.05))
                #     (1. + 0.5 * (0.2 * np.sin(dt * 8) - 0.3 * np.cos(4 * dt + np.pi / 2) + 0.2 * np.sin(dt * 5)**2 - dt * 0.05))
                # self.demand = self.cda2g*np.sqrt(sp.H[NN]-self.z)
            elif T> self.burstTime:
                x11 = (T-self.burstTime) / self.duration
                cda =  self.cda * 0.5 * (1 - np.cos(x11 * np.pi))
            else:
                cda = 0
            if sp.H[NN]-self.z<=0:
                    self.demand = 0
            else:
                self.demand = cda*np.sqrt(sp.H[NN]-self.z)
        else:      
            if T <self.burstTime:
                self.demand = 0
            elif T > self.burstTime and T <self.burstTime +self.duration:
                x11 = (T - self.burstTime) / self.duration
                self.cda2g=self.cda2g0*0.5 * (1 - np.cos(x11 * np.pi))
                # if sp.H[NN]-self.z<=0:
                #     self.demand = 0
                # else:
                self.demand = self.cda2g*np.sqrt(sp.H[NN]-self.z)
            else:
                # if sp.H[NN]-self.z<=0:
                #     self.demand = 0
                # else:
                    self.demand = self.cda2g0*np.sqrt(sp.H[NN]-self.z)
        return self.demand
    def moc(self, T):
        # sp = self.start_p[0]
        # ep = self.end_p[0]
        # NN = sp.NN
        # self.get_demand(T)
        # sp.HP[NN] = (sp.QCP(NN) - ep.QCM(0) - self.demand) / (sp.CQP(NN) + ep.CQM(0))
        # sp.QP[NN] = sp.QCP(NN) - sp.CQP(NN) * sp.HP[NN]
        # ep.QP[0] = sp.QP[NN] - self.demand
        # ep.HP[0] = sp.HP[NN]
        self.get_demand(T)
        self.demands.append(self.demand)
        l = 0.0
        m = 0.0
        B = 0.0
        for sp in self.start_p:
            l = l + sp.QCP(sp.NN)
            B = B + sp.CQP(sp.NN)
        for ep in self.end_p:
            m = m - ep.QCM(0)
            B = B + ep.CQM(0)
        hp = (l + m-self.demand) / B
        for sp in self.start_p:
            sp.HP[sp.NN] = hp
            sp.QP[sp.NN] = sp.QCP(sp.NN) - sp.CQP(sp.NN) * sp.HP[sp.NN]
        for ep in self.end_p:
            ep.HP[0] = hp
            ep.QP[0] = ep.QCM(0) + ep.CQM(0) * ep.HP[0]

    def plot_demands(self):
        plt.plot(self.demands)
        plt.show()

class Coordinate(SysParam):
    number = 0

    def __init__(self, node_id, X_coord, Y_coord,epa_node=None):
        self.id = id
        self.node_id = node_id
        self.X_coord = X_coord
        self.Y_coord = Y_coord
        self.epa_node=epa_node

    def connect(self, nodes):
        node = nodes[self.node_id]
        node.X_coord = self.X_coord
        node.Y_coord = self.Y_coord
        node.epa_node=self.epa_node

def error(message):
    print(message)
    os.system('pause')
    exit(0)
