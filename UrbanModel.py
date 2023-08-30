import numpy
from scipy.interpolate import interp1d
from DragLength import Drag_Coef,Length_Scale,Length_Scale_StabilityCorrection
from Resistance_Functions import Ressitance_Calculations
from ColModel import ColumnModelCal
import copy

"""
Urban canopy model
Developed by Mohsen Moradi and Amir A. Aliabadi
Atmospheric Innovations Research (AIR) Laboratory, University of Guelph, Guelph, Canada
Last update: May 2021
"""

class UCM_Def():

    def __init__(self,vx_init,vy_init,tke_init,th_init,qn_init,pr_init,rho_init,th_ref,nz):

        class VerticalProfUrban_Def():
            pass
        self.VerticalProfUrban = VerticalProfUrban_Def()
        self.VerticalProfUrban.vx = [vx_init for i in range(nz)]
        self.VerticalProfUrban.vy = [vy_init for i in range(nz)]
        self.VerticalProfUrban.s = [numpy.sqrt(vx_init**2+vy_init**2) for i in range(nz)]
        self.VerticalProfUrban.tke = [tke_init for i in range(nz)]
        self.VerticalProfUrban.th = [th_init for i in range(nz)]
        self.VerticalProfUrban.qn = [qn_init for i in range(nz)]
        self.VerticalProfUrban.presProf = [pr_init for i in range(nz)]
        self.VerticalProfUrban.th_ref = [th_ref for i in range(nz)]
        self.VerticalProfUrban.rho = [rho_init for i in range(nz)]
        self.VerticalProfUrban.Hflux = [rho_init for i in range(nz)]
        self.VerticalProfUrban.LEflux = [rho_init for i in range(nz)]

        class Vx_eq_Def():
            pass
        self.Vx_eq = Vx_eq_Def()
        self.Vx_eq.srex = [0 for i in range(nz)]
        self.Vx_eq.srex_h = [0 for i in range(nz)]
        self.Vx_eq.srex_dpdx = [0 for i in range(nz)]
        self.Vx_eq.srim = [0 for i in range(nz)]
        self.Vx_eq.srim_v = [0 for i in range(nz)]
        self.Vx_eq.srim_veg = [0 for i in range(nz)]
        self.Vx_eq.diffusion = [0 for i in range(nz+1)]

        class Vy_eq_Def():
            pass
        self.Vy_eq = Vy_eq_Def()
        self.Vy_eq.srex = [0 for i in range(nz)]
        self.Vy_eq.srex_h = [0 for i in range(nz)]
        self.Vy_eq.srex_dpdx = [0 for i in range(nz)]
        self.Vy_eq.srim = [0 for i in range(nz)]
        self.Vy_eq.srim_v = [0 for i in range(nz)]
        self.Vy_eq.srim_veg = [0 for i in range(nz)]
        self.Vy_eq.diffusion = [0 for i in range(nz+1)]

        class thb_Def():
            pass
        self.thb = thb_Def()
        self.thb.ground_imp = 0
        self.thb.ground_bare = 0
        self.thb.ground_veg = 0
        self.thb.tree = 0
        self.thb.roof_imp = 0
        self.thb.roof_veg = 0

        class tvb_Def():
            pass
        self.tvb = tvb_Def()
        self.tvb.wall_sun = 0
        self.tvb.wall_shade = 0

        class qhb_Def():
            pass
        self.qhb = qhb_Def()
        self.qhb.ground_imp = 0
        self.qhb.ground_bare = 0
        self.qhb.ground_veg = 0
        self.qhb.tree = 0
        self.qhb.roof_imp = 0
        self.qhb.roof_veg = 0

    def UCMCal(self,TemperatureR,TemperatureC,FractionsGround,FractionsRoof,ForcTemp,Forcq_Top,ForcWindSpeed,ForcWindDir,
               Geometry_m,ParCalculation,ParVegGround,ParVegRoof,thb_roof,thb_ground,tvb_wall,qhb_roof,qhb_ground,
               Hflux_canyon,LEflux_canyon,Hflux_roof,LEflux_roof,geometry,ParTree,ColParam,BEM,Ustar,Rural_Model_name):
        """
        ------
        INPUT:
        TemperatureR: Temperature of the exterior surfaces of roof [K]
        TemperatureC: Temperature of the exterior surfaces of urban canyon [K]
        FractionsGround: Fractions of ground covered by vegetation, impervious, and bare soil
        FractionsRoof: Fractions of roof covered by vegetation and impervious
        ForcTemp: Temperature boundary conditions at the top of the domain [K]
        Forcq_Top: Specific humidity boundary conditions at the top of the domain [kg kg^-1]
        ForcWindSpeed: Wind speed boundary conditions at the top of the domain [m s^-1]. Note: It will be used if any
                       wind speed at the top is given, otherwise the pressure gradient is used
        ForcWindDir: Wind direction at the top of the domain [deg].
                     Note: If the "Rural_Model_name" is "Forcing_extFile", it is exactly at the top, otherwise it is from epw file.
        Geometry_m: Geometric parameters
        ParCalculation: General calculation parameters
        ParVegGround: Ground vegetation parameters
        ParVegRoof: Roof vegetation parameters
        thb_roof: Sink/source terms in temperature equation in 1-D model caused by roof
        thb_ground: Sink/source terms in temperature equation in 1-D model caused by ground
        tvb_wall: Sink/source terms in temperature equation in 1-D model caused by walls
        qhb_roof: Sink/source terms in humidity equation in 1-D model caused by roof
        qhb_ground: Sink/source terms in humidity equation in 1-D model caused by ground
        geometry: Normalized geometric parameters
        ParTree: Trees parameter
        ColParam: Column model parameters
        BEM: Building energy variables
        -------
        OUTPUT:
        VerticalProfUrban: Vertical profile of temperature, humidity, wind speed, and turbulent kinetic energy
        """

        ForcingVariable = [ForcTemp,ForcWindSpeed,ForcWindDir,Forcq_Top]

        self.thb.ground_imp = copy.copy(thb_ground.ground_imp)    # source term in temperature equation originated from impervious ground [K m s^-1]
        self.thb.ground_bare = copy.copy(thb_ground.ground_bare)  # source term in temperature equation originated from bare ground [K m s^-1]
        self.thb.ground_veg = copy.copy(thb_ground.ground_veg)    # source term in temperature equation originated from vegetated ground [K m s^-1]
        self.thb.tree = copy.copy(thb_ground.tree)                # source term in temperature equation originated from trees [K m s^-1]
        self.thb.roof_imp = copy.copy(thb_roof.roof_imp)          # source term in temperature equation originated from impervious roof [K m s^-1]
        self.thb.roof_veg = copy.copy(thb_roof.roof_veg)          # source term in temperature equation originated from vegetated roof [K m s^-1]

        self.tvb.wall_sun = copy.copy(tvb_wall.Hwsun_z)           # source term in temperature equation originated from sunlit wall [K m s^-1]
        self.tvb.wall_shade = copy.copy(tvb_wall.Hwshade_z)       # source term in temperature equation originated from shaded wall [K m s^-1]

        self.qhb.ground_imp = copy.copy(qhb_ground.ground_imp)    # source term in humidity equation originated from impervious ground [kg kg^-1 m s^-1]
        self.qhb.ground_bare = copy.copy(qhb_ground.ground_bare)  # source term in humidity equation originated from bare ground [kg kg^-1 m s^-1]
        self.qhb.ground_veg = copy.copy(qhb_ground.ground_veg)    # source term in humidity equation originated from vegetated ground [kg kg^-1 m s^-1]
        self.qhb.tree = copy.copy(qhb_ground.tree)                # source term in humidity equation originated from trees [kg kg^-1 m s^-1]
        self.qhb.roof_imp = copy.copy(qhb_roof.roof_imp)          # source term in humidity equation originated from impervious roof [kg kg^-1 m s^-1]
        self.qhb.roof_veg = copy.copy(qhb_roof.roof_veg)          # source term in humidity equation originated from vegetated roof [kg kg^-1 m s^-1]

        # Boolean operator for presence and absence of trees, vegetated ground, bare ground, and impervious ground
        Ctree = int(ParTree.trees == 1)
        Cgveg = int(FractionsGround.fveg > 0)
        Cgbare = int(FractionsGround.fbare > 0)
        Cgimp = int(FractionsGround.fimp > 0)
        Croof = int(FractionsRoof.fimp+FractionsRoof.fveg > 0)

        if FractionsRoof.fveg > 0:
            hc_roof = copy.copy(ParVegRoof.hc)
        else:
            hc_roof = 0

        # "ss(z)" Probability that a building has a height equal to z (In the current version of the model a simple
        # canyon is considered so this probability is one at building average height h mean (nz_u) but zero elsewhere.)
        ss = numpy.zeros(Geometry_m.nz + 1)
        ss[Geometry_m.nz_u] = 1

        # "pb(z)" Probability that a building has a height greater or equal to z (In the current version of the model a simple
        # canyon is considered. So, "pb" is one within the canyon and zero above the canyon.)
        pb = numpy.zeros(Geometry_m.nz + 1)
        for i in range(0, Geometry_m.nz + 1):
            if i <= Geometry_m.nz_u:
                pb[i] = 1
            else:
                pb[i] = 0

        # Volume fraction of air in each urban unit cell
        vol = numpy.zeros(Geometry_m.nz)
        # Fraction of air at the interface between cells
        sf = numpy.zeros(Geometry_m.nz+1)
        for i in range(0, Geometry_m.nz):
            vol[i] = 1 - Geometry_m.lambdap * pb[i]
            # "sf" is calculated from Nazarian's code (https://github.com/nenazarian/MLUCM/blob/master/Column_Model/column_lkPro.f90)
            sf[i] = 1 - Geometry_m.lambdap * ss[i]
        sf[Geometry_m.nz] = 1

        # Interpolate LAD profile
        f_LAD = interp1d(ColParam.h_LAD, ColParam.LAD)

        # Length of an infinite urban unit (length of flat earth corresponding to an infinite urban unit)
        Lurban = Geometry_m.Width_canyon + Geometry_m.Width_roof
        # Calculate LAI per unit area of the canyon [m^2 m^-2]
        LAI_canyon = sum(ColParam.LAD)*Geometry_m.dz
        # Length of projected two side of leaf area in an infinite urban unit
        Lleaf = 2*Geometry_m.Width_canyon*LAI_canyon

        # Coefficient for the destruction of turbulent dissipation rate
        Ceps = 1 / 1.14

        # Coefficient used in the equation of diffusion coefficient
        Ck = 0.4

        # Calculate section drag coefficient (Cdrag) due to buildings
        # Displacement height [m]
        disp = Geometry_m.Height_canyon*(Geometry_m.lambdap)**(0.15)
        # Calculate total urban heat flux per unit flat area of the earth [W m^-2]
        hfx = +(FractionsGround.fimp*(Hflux_canyon.HfluxGroundImp+LEflux_canyon.LEfluxGroundImp)+
                FractionsGround.fbare*(Hflux_canyon.HfluxGroundBare+LEflux_canyon.LEfluxGroundBare)+
                FractionsGround.fveg*(Hflux_canyon.HfluxGroundVeg+LEflux_canyon.LEfluxGroundVeg))*(Geometry_m.Width_canyon/Lurban)
        for i in range(0,Geometry_m.nz_u+1):
            hfx = hfx + (tvb_wall.Hfluxwsun+tvb_wall.Hfluxwshade)*(Geometry_m.dz/Lurban) + \
                  (FractionsRoof.fimp*(Hflux_roof.HfluxRoofImp+LEflux_roof.LEfluxRoofImp)+
                   FractionsRoof.fveg*(Hflux_roof.HfluxRoofVeg+LEflux_roof.LEfluxRoofVeg))*(Geometry_m.Width_roof/Lurban) + \
                  (Hflux_canyon.HfluxTree+LEflux_canyon.LEfluxTree)*(Lleaf/Lurban)

        # Calculate the length scales
        dls , dlk = Length_Scale_StabilityCorrection(Geometry_m.nz,Geometry_m.z,Geometry_m.Height_canyon,Ceps,Ck,hfx,
                                                     self.Vx_eq,self.Vy_eq,self.VerticalProfUrban,Geometry_m.lambdaf,
                                                     disp,Geometry_m.dz,vol)
        # Calculate the drag coefficient
        Cdrag = Drag_Coef(Geometry_m.nz,Geometry_m.lambdaf,pb)

        # Calculate total urban sensible heat flux [W m^-2]
        self.UrbanFlux_H = 0
        self.UrbanFlux_H = (FractionsGround.fimp*Hflux_canyon.HfluxGroundImp + FractionsGround.fbare*Hflux_canyon.HfluxGroundBare +
                            FractionsGround.fveg*Hflux_canyon.HfluxGroundVeg)*(Geometry_m.Width_canyon/Lurban) + \
                           (FractionsRoof.fimp*(Hflux_roof.HfluxRoofImp)+FractionsRoof.fveg*(Hflux_roof.HfluxRoofVeg))*(Geometry_m.Width_roof/Lurban) +\
                           (Hflux_canyon.HfluxTree)*(Lleaf/Lurban)
        for i in range(0,Geometry_m.nz_u+1):
            self.UrbanFlux_H = self.UrbanFlux_H + (tvb_wall.Hfluxwsun+tvb_wall.Hfluxwshade)*(Geometry_m.dz/Lurban)


        # Calculate total urban latent heat flux [W m^-2]
        self.UrbanFlux_LE = 0
        self.UrbanFlux_LE = (FractionsGround.fimp*LEflux_canyon.LEfluxGroundImp + FractionsGround.fbare*LEflux_canyon.LEfluxGroundBare +
                             FractionsGround.fveg*LEflux_canyon.LEfluxGroundVeg)*(Geometry_m.Width_canyon/Lurban) + \
                            (FractionsRoof.fimp*LEflux_roof.LEfluxRoofImp + FractionsRoof.fveg*LEflux_roof.LEfluxRoofVeg)*(Geometry_m.Width_roof/Lurban) + \
                            LEflux_canyon.LEfluxTree*(Lleaf/Lurban)


        # Gas constant [J kg^-1 K^-1]
        R = 287.0
        # Air Specific heat [J kg^-1 K^-1]
        Cp = 1004.0
        Treal = [self.VerticalProfUrban.th[i]*(self.VerticalProfUrban.presProf[i]/self.VerticalProfUrban.presProf[0])**(R/Cp)
                 for i in range(Geometry_m.nz)]
        # Calculate density profile [kg m^-3]
        self.VerticalProfUrban.rho = [self.VerticalProfUrban.presProf[i]/(R*Treal[i]) for i in range(Geometry_m.nz)]

        # Calculate total waste heat from buildings
        SensHt_HVAC = 0
        for i in range(0, len(BEM)):
            SensHt_HVAC = SensHt_HVAC + BEM[i].building.sensWaste

        # Calculate aerodynamic roughness for the ground and roof
        ResistanceCal = Ressitance_Calculations()
        z0_roof, _zoh_, zom_ground_roof, _zoh_ground_, _disp_h_, _zom_H_, _zom_L_, _zoh_H_, _zoh_L_, _d_H_, _d_L_, _zom_other_= \
            ResistanceCal.Urban_roughness(hc_roof, 0, 0, 0, Croof)
        z0_road, _zoh_, zom_ground_ground, _zoh_ground_, _disp_h_, _zom_H_, _zom_L_, _zoh_H_, _zoh_L_, _d_H_, _d_L_, _zom_other_ = \
            ResistanceCal.Urban_roughness(Ctree*Geometry_m.Height_tree, Cgveg*ParVegGround.hc, Cgbare, Cgimp, 0)

        # Call the 1-D model
        vx, vy, tke, th, qn, presProf, th_eq, vx_eq, vy_eq, tke_eq, qn_eq, self.ColumnModel_OtherParam = \
            ColumnModelCal(zom_ground_ground,zom_ground_roof,Ceps,Cdrag,Ck,self.thb,self.qhb,self.tvb,FractionsGround,
                           FractionsRoof,TemperatureC,TemperatureR,ForcingVariable,self.VerticalProfUrban,Geometry_m,geometry,
                           ColParam,f_LAD,dlk,dls,pb,ss,sf,vol,Cp,Ustar,SensHt_HVAC,ParCalculation.dts,Rural_Model_name)

        self.VerticalProfUrban.vx = vx
        self.VerticalProfUrban.vy = vy
        self.VerticalProfUrban.s = [numpy.sqrt(vx[i]**2 + vy[i]**2) for i in range(len(vx))]
        self.VerticalProfUrban.tke = tke
        self.VerticalProfUrban.th = th
        self.VerticalProfUrban.qn = qn
        self.VerticalProfUrban.presProf = presProf
        self.VerticalProfUrban.Hflux = [self.VerticalProfUrban.rho[i]*Cp*self.ColumnModel_OtherParam.wth[i] for i in range(Geometry_m.nz)]
        self.VerticalProfUrban.LEflux = [self.VerticalProfUrban.rho[i]*ParCalculation.Lv*self.ColumnModel_OtherParam.wqn[i] for i in range(Geometry_m.nz)]

        self.T_eq = copy.copy(th_eq)
        self.Vx_eq = copy.copy(vx_eq)
        self.Vy_eq = copy.copy(vy_eq)
        self.TKE_eq = copy.copy(tke_eq)
        self.Qn_eq = copy.copy(qn_eq)