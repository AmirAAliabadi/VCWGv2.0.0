import numpy
import math
import copy

"""
Formulate linear system of equations to solve; make matrix of coefficients A and the right hand side RHS vector
Developed by Mohsen Moradi and Amir A. Aliabadi
Atmospheric Innovations Research (AIR) Laboratory, University of Guelph, Guelph, Canada
Last update: May 2021
Originally developed by Alberto Martilli, Scott Krayenhoff
"""

class Diff:
    def __init__(self,nz,dt,sf,vol,dz,rho):
        self.nz = nz
        self.dt = dt
        self.sf = sf
        self.vol = vol
        self.dz = dz
        self.rho = rho

    def Solver(self,nzm,nz,iz1,izf,dt,rho,co,cd,aa,bb,sf,vl,dz):

        def invert(nzm,nn,A,C):

            # C = copy.copy(c)
            # A = copy.copy(a)
            for In in range(nn-2,-1,-1):
                C[In] = C[In]-A[In][2]*C[In+1]/A[In+1][1]
                A[In][1] = A[In][1]-A[In][2]*A[In+1][0]/A[In+1][1]

            for In in range(1,nn):
                C[In] = C[In]-A[In][0]*C[In-1]/A[In-1][1]

            x = numpy.zeros(nn)
            for In in range(0,nn):
                x[In] = C[In]/A[In][1]
            return x

        rhoz = numpy.zeros(nz)
        rhoz[0] = rho[0]
        for iz in range(1,nz):
            rhoz[iz] = (rho[iz]*dz+rho[iz-1]*dz)/(dz+dz)

        cddz = numpy.zeros(nz+1)
        cddz[0] = sf[0]*rho[0]*cd[0]/dz
        for iz in range(1,nz):
            cddz[iz] = rhoz[iz]*sf[iz]*cd[iz]/((dz+dz)/2)
        if izf > 1:
            cddz[nz] = sf[nz]*rho[nz-1]*cd[nz]/dz
        else:
            cddz[nz] = 0

        a = numpy.zeros((nz, 3))
        c = numpy.zeros(nz)
        if iz1 > 1:
            a[0][0] = 0
            a[0][1] = 1
            a[0][2] = 0
            c[0] = co[0]

        dzv = numpy.zeros(nz)
        for iz in range(iz1-1,nz-izf+1):
            dzv[iz] = vl[iz]*dz
            a[iz][0] = -(1/rho[iz])*cddz[iz]*dt/dzv[iz]
            a[iz][1] = (1/rho[iz])*dt*((1/dzv[iz])*(cddz[iz]+cddz[iz+1]))+1-aa[iz]*dt
            a[iz][2] = -(1/rho[iz])*cddz[iz+1]*dt/dzv[iz]
            c[iz] = co[iz]+bb[iz]*dt

        if izf == 1:
            dzv[nz-1] = vl[nz-1]*dz
            a[nz-1][0] = -cddz[nz-1]*dt/dzv[nz-1]
            a[nz-1][1] = 1+dt*(cddz[nz-1])/dzv[nz-1]-aa[nz-1]*dt
            a[nz-1][2] = 0.
            c[nz-1] = co[nz-1]+bb[nz-1]*dt
        else:
            a[nz - 1][0] = 0
            a[nz - 1][1] = 1.
            a[nz - 1][2] = 0.
            c[nz - 1] = co[nz - 1]

        co_new = invert(nzm,nz,a,c)

        fc = numpy.zeros(nz)
        for iz in range(0,iz1):
            fc[iz] = 0
        for iz in range(iz1,nz):
            fc[iz] = -(cddz[iz] * (co[iz] - co[iz - 1])) / rho[iz]

        df = numpy.zeros(nz)
        for iz in range(0,iz1):
            df[iz] = 0
        for iz in range(iz1,nz-izf+1):
            dzv = vl[iz]*dz
            if iz < nz-1:
                df[iz] = +(co[iz-1]*cddz[iz]-co[iz]*(cddz[iz]+cddz[iz+1])+co[iz+1]*cddz[iz+1])/dzv/rho[iz]
            else:
                df[iz] = +(co[iz-1]*cddz[iz]-co[iz]*(cddz[iz]+cddz[iz+1]))/dzv/rho[iz]

        for iz in range(nz-izf-1,nz):
            df[iz] = 0

        return co_new,fc,df


