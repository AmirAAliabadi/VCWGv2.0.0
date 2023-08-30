import numpy
import math

"""
Calculate sink and source terms associated with the presence of buildings in the 1D model for momentum, heat, and TKE  
Developed by Mohsen Moradi and Amir A. Aliabadi
Atmospheric Innovations Research (AIR) Laboratory, University of Guelph, Guelph, Canada
Last update: February 2020
Originally developed by Alberto Martilli, Scott Krayenhoff, and Negin Nazarian
"""

# explicit and implicit terms for the building
class BuildingCol:

    # number of street direction (assumed to be one)
    nd = 1.0
    roadfrac = 1
    def __init__(self,nz,dz,dt,vol,lambdap,lambdaf,hmean,Ck,Cp,th_ref,vx,vy,th,Cdrag,rho,nz_u,pb,ss
                 ,z0_road,z0_roof,SensHt_HVAC,HVAC_street_frac,HVAC_atm_frac,windMin):
        self.nz = nz                     # Number of grid points in vertical column
        self.dz = dz                     # Grid resolution [m]
        self.dt = dt                     # Time step [s]
        self.vol = vol                   # Fraction of air to total volume in each urban unit cell
        self.lambdap = lambdap           # Plan area fraction of buildings
        self.lambdaf = lambdaf           # Ratio of wall area facing ambient wind to plan area
        self.hmean = hmean               # Average building height [m]
        self.Ck = Ck                     # Coefficient used in the equation of diffusion coefficient (kappa)
        self.Cp = Cp                     # Heat capacity of dry air [J kg^-1 K^-1]
        self.th_ref = th_ref             # Reference potential temperature [K]
        self.vx = vx                     # x component of horizontal wind speed [m s^-1]
        self.vy = vy                     # y component of horizontal wind speed [m s^-1]
        self.th = th                     # Potential temperature [K]
        self.Cdrag = Cdrag               # Drag coefficient due to buildings (sectional drag coefficient)
        self.rho = rho                   # density profile [kg m^-3]
        self.nz_u = nz_u                 # Number of grid points within the canyon
        self.pb = pb                     # Probability distribution of building height (assumed to be one
                                         # within the canyon and zero above the canyon)
        self.ss = ss                     # Probability that a building has a height equal to z (assumed to be
                                         # one at average building height and zero the other heights)
        self.z0_road = z0_road           # Road roughness [m]
        self.z0_roof = z0_roof           # Roof roughness [m]
        self.SensHt_HVAC = SensHt_HVAC            # Sensible waste heat from building
        self.HVAC_street_frac = HVAC_street_frac  # Fraction of Sensible waste heat from building released into the atmosphere at street level
        self.HVAC_atm_frac = HVAC_atm_frac        # Fraction of sensible waste heat from building released into the atmosphere
        self.windMin = windMin                    # minimum wind speed

        # Define explicit and implicit parts of source and sink terms due to building
        self.srex_vx_h = numpy.zeros(self.nz)  # Term in momentum equation
        self.srex_vy_h = numpy.zeros(self.nz)  # Term in momentum equation
        self.srex_tke_h = numpy.zeros(self.nz) # Term in turbulent kinetic energy equation
        self.srex_th_h = numpy.zeros(self.nz)  # Term in energy equation
        self.srim_vx_v = numpy.zeros(self.nz)  # Term in momentum equation
        self.srim_vy_v = numpy.zeros(self.nz)  # Term in momentum equation
        self.srex_tke_v = numpy.zeros(self.nz) # Term in turbulent kinetic energy equation
        self.srim_th_v = numpy.zeros(self.nz)  # Term in energy equation
        self.srex_th_v = numpy.zeros(self.nz)  # Term in energy equation
        self.srex_qn_h = numpy.zeros(self.nz)  # Term in humidity equation


    def BuildingDrag_UTC(self,thb,qhb,tvb,FractionsGround,FractionsRoof,TWallSun,TWallShade,TGround,TRoof):

        """
        ------
        INPUT:
        thb: Sink/source terms in temperature equation in 1-D model caused by roof and ground
        qhb: Sink/source terms in specific humidity equation in 1-D model
        tvb: Sink/source terms in temperature equation in 1-D model caused by walls
        FractionsGround:
        FractionsRoof:
        TWallSun: Sunlit wall temperature [K]
        TWallShade: Shaded wall temperature [K]
        TGround: Weighted-average ground temperature [K]
        TRoof: Weighted-average roof temperature [K]
        -------
        OUTPUT:
        srex_vx_h: Explicite term in momentum equation [m s^-2]
        srex_vy_h: Explicite term in momentum equation [m s^-2]
        srex_tke_h: Explicite term in turbulent kinetic energy equation [m s^-3]
        srex_th_h: Explicite term in energy equation [K s^-1]
        srim_vx_v: Implicite term in momentum equation [s^-1]
        srim_vy_v: Implicite term in momentum equation [s^-1]
        srex_tke_v: Explicite term in turbulent kinetic energy equation [m^2 s^-3]
        srim_th_v: Implicite term in energy equation [s^-1]
        srex_th_v: Explicite term in energy equation [K s^-1]
        srex_qn_h: Explicite term in humidity equation [kgw kga^-1 s^-1]
        """

        #-----------------------------
        # Fluxes from Ground and Roof
        #-----------------------------
        # Calculate momentum flux from ground
        self.FluxFlatG = self.Flux_Flat(self.z0_road, self.vx[0], self.vy[0],self.th[0], self.th_ref[0],TGround, self.windMin)
        uhb = self.FluxFlatG[0] # momentum flux in x direction [m^2 s^-2]
        vhb = self.FluxFlatG[1] # momentum flux in y direction [m^2 s^-2]
        ehb = self.FluxFlatG[2] # tke flux [m^2 s^-3]

        # Term in momentum equation [m s^-2]
        self.srex_vx_h[0] = (uhb / self.nd) / self.dz * (1 - self.lambdap) / self.vol[0] * self.roadfrac
        # Term in momentum equation [m s^-2]
        self.srex_vy_h[0] = (vhb / self.nd) / self.dz * (1 - self.lambdap) / self.vol[0] * self.roadfrac
        # Term in energy equation [K s^-1]
        self.srex_th_h[0] = (thb.ground_imp/self.nd) / self.dz * (1 - self.lambdap) / self.vol[0] * FractionsGround.fimp + \
                            (thb.ground_bare/self.nd) / self.dz * (1 - self.lambdap) / self.vol[0] * FractionsGround.fbare + \
                            (thb.ground_veg/self.nd) / self.dz * (1 - self.lambdap) / self.vol[0] * FractionsGround.fveg + \
                            self.HVAC_atm_frac*self.HVAC_street_frac*(self.SensHt_HVAC/(self.rho[0]*self.Cp)/self.dz)*self.lambdap/(1-self.lambdap)

        # Term in turbulent kinetic energy equation [m s^-3]
        self.srex_tke_h[0] = (ehb/self.nd)/ self.dz * (1 - self.lambdap) / self.vol[0] * self.roadfrac
        # Term in humidity equation [kgw kga^-1 s^-1]
        self.srex_qn_h[0] = (qhb.ground_imp/self.nd) / self.dz * (1 - self.lambdap) / self.vol[0] * FractionsGround.fimp + \
                            (qhb.ground_bare/self.nd) / self.dz * (1 - self.lambdap) / self.vol[0] * FractionsGround.fbare + \
                            (qhb.ground_veg/self.nd) / self.dz * (1 - self.lambdap) / self.vol[0] * FractionsGround.fveg

        # Calculate fluxes from horizontal surfaces other than the ground
        for i in range(0, self.nz_u + 1):
            # At roof level for simple and non-probabilistic canyon
            if self.ss[i] > 0:
                self.FluxFlatR = self.Flux_Flat(self.z0_roof, self.vx[i],self.vy[i],self.th[i],self.th_ref[i], TRoof, self.windMin)
                uhb = self.FluxFlatR[0]
                vhb = self.FluxFlatR[1]
                ehb = self.FluxFlatR[2]

                thb_Roofs_Imp = thb.roof_imp
                thb_Roofs_Veg = thb.roof_veg
                qhb_Roofs_Imp = qhb.roof_imp
                qhb_Roofs_Veg = qhb.roof_veg

            else:
                uhb = 0
                vhb = 0
                ehb = 0
                thb_Roofs_Imp = 0
                thb_Roofs_Veg = 0
                qhb_Roofs_Imp = 0
                qhb_Roofs_Veg = 0

            # Term in momentum equation [m s^-2]
            self.srex_vx_h[i] += (uhb / self.nd) * (self.ss[i] * self.lambdap / self.vol[i] / self.dz)
            # Term in momentum equation [m s^-2]
            self.srex_vy_h[i] += (vhb / self.nd) * (self.ss[i] * self.lambdap / self.vol[i] / self.dz)
            # Term in energy equation [K s^-1]
            self.srex_th_h[i] += (thb_Roofs_Imp / self.nd) * (self.ss[i] * self.lambdap / self.vol[i] / self.dz) * FractionsRoof.fimp + \
                                 (thb_Roofs_Veg / self.nd) * (self.ss[i] * self.lambdap / self.vol[i] / self.dz) * FractionsRoof.fveg + \
                                 self.HVAC_atm_frac*(1-self.HVAC_street_frac)*(self.SensHt_HVAC /(self.rho[i]*self.Cp)/self.dz)*self.ss[i]*self.lambdap/(1 - self.lambdap)


            # Term in turbulent kinetic energy equation [m s^-3]
            self.srex_tke_h[i] += (ehb/self.nd) * (self.ss[i] * self.lambdap / self.vol[i] / self.dz)
            # Term in humidity equation [kgw kga^-1 s^-1]
            self.srex_qn_h[i] += (qhb_Roofs_Imp / self.nd) * (self.ss[i] * self.lambdap / self.vol[i] / self.dz) * FractionsRoof.fimp + \
                                 (qhb_Roofs_Veg / self.nd) * (self.ss[i] * self.lambdap / self.vol[i] / self.dz) * FractionsRoof.fveg

        # ------------------
        # Fluxes from Walls
        # ------------------
        for i in range(0, self.nz_u):

            # Calculate momentum and tke fluxes from wall
            uva, vva, _uvb_, _vvb_, _tva_, _tvb_, evb = \
                self.Flux_Wall(self.vx[i], self.vy[i], self.th[i], self.Cdrag[i], TWallSun, self.rho[i],self.windMin)

            # Term in momentum equation [s^-1]
            self.srim_vx_v[i] = uva * self.lambdaf * self.pb[i] / max(1e-6, self.hmean) / self.vol[i] / self.nd
            # Term in momentum equation [s^-1]
            self.srim_vy_v[i] = vva * self.lambdaf * self.pb[i] / max(1e-6, self.hmean) / self.vol[i] / self.nd
            # Term in turbulent kinetic energy equation [m^2 s^-3]
            self.srex_tke_v[i] = evb * self.lambdaf * self.pb[i] / max(1e-6, self.hmean) / self.vol[i] / self.nd

            # Option 1: Meili et al, 2020
            '''
            w = 0
            u = numpy.sqrt(self.vx[i]**2+self.vy[i]**2)
            tva_sun = 0
            tva_shade = 0
            RES = self.Cp*self.rho[i]*(11.8+4.2*numpy.sqrt(u**2+w**2))**(-1)
            tvb_sun = (TWallSun-self.th[i])/RES
            tvb_shade = (TWallShade-self.th[i])/RES
            '''
            # Option 2: Martilli et al, 2001
            tva_sun = 0
            tva_shade = 0
            tvb_sun = tvb.wall_sun[i]
            tvb_shade = tvb.wall_shade[i]

            # Implicit Term in energy equation
            self.srim_th_v[i] = tva_sun * self.lambdaf * self.pb[i] / max(1e-6, self.hmean) / self.vol[i] / self.nd + \
                                tva_shade * self.lambdaf * self.pb[i] / max(1e-6, self.hmean) / self.vol[i] / self.nd
            # Explicit Term in energy equation [K s^-1]
            self.srex_th_v[i] = tvb_sun * self.lambdaf * self.pb[i] / max(1e-6, self.hmean) / self.vol[i] / self.nd + \
                                tvb_shade * self.lambdaf * self.pb[i] / max(1e-6, self.hmean) / self.vol[i] / self.nd


    def Flux_Flat(self,z0,vx,vy,th_air,th_ref,pts,windMin):

        Utot = (vx**2+vy**2)**0.5
        zz = self.dz/2

        Utot = max(Utot,windMin)

        # Compute bulk Richardson number using near surface temperatures
        Ri = 2 * 9.81 * zz * (th_air - pts) / ((th_air + pts) * (Utot ** 2))
        # Calculation from Louis, 1979 (eq. 11 and 12)
        b = 9.4
        cm = 7.4
        ch = 5.3
        R = 0.74
        a = self.Ck/math.log(zz/z0)

        if Ri > 0:
            fm = 1/(1+0.5*b*Ri)**2
            fh = fm
        else:
            c = b*cm*a*a*(zz/z0)**0.5
            fm = 1-b*Ri/(1+c*(-Ri)**0.5)
            c = c*ch/cm
            fh = 1-b*Ri/(1+c*(-Ri)**0.5)

        fbuw = -a*a*Utot*Utot*fm
        fbpt = -a*a*Utot*(th_air-pts)*fh/R
        ustar = (-fbuw)**0.5
        tstar = -fbpt/ustar

        # x component momentum flux from horizontal surfaces [m^2 s^-2]
        uhb = -ustar*ustar*vx/Utot
        # y component momentum flux from horizontal surfaces [m^2 s^-2]
        vhb = -ustar*ustar*vy/Utot
        # Heat flux from horizontal surfaces [K m s^-1]
        thb = -ustar*tstar
        # Turbulent flux of TKE from horizontal surfaces [m^2 s^-3]
        ehb = -(9.81/th_ref)*ustar*tstar

        return uhb,vhb,ehb,thb,ustar

    def Flux_Wall(self,vx,vy,th,Cdrag,ptw,rho,windMin):

        vett = (vx**2+vy**2)**0.5
        vett = max(vett, windMin)
        # Implicit term of x component momentum flux from vertical surfaces [m s^-1]
        uva = -Cdrag*vett
        # Implicit term of y component momentum flux from vertical surfaces [m s^-1]
        vva = -Cdrag*vett
        # Explicit term of x component momentum flux from vertical surfaces
        uvb = 0
        # Explicit term of y component momentum flux from vertical surfaces
        vvb = 0

        # Calculation for S_theta_wall in eq. 5.5 (Krayenhoff, PhD thesis)
        # Convective heat transfer coefficient [W K^-1 m^-2]
        hc = 5.678*(1.09+0.23*(vett/0.3048))
        # Using energy balance for a control volume inside the urban unit, the convective heat transfer coefficient should be limited
        # hc must be less than (rho * cp / dt) * [(1-lambdap) * Hmean / (4 * lambdaf * dz)]
        if hc > ((rho*self.Cp/self.dt)*((1-self.lambdap)*self.hmean)/(4*self.lambdaf*self.dz)):
            hc = (rho*self.Cp/self.dt)*((1-self.lambdap)*self.hmean)/(4*self.lambdaf*self.dz)
        # Term in energy equation [K m s^-1]
        tvb = (hc/(rho*self.Cp))*(ptw-th)
        tva = 0

        evb = Cdrag*(abs(vett)**3)

        return uva,vva, uvb, vvb, tva, tvb, evb
