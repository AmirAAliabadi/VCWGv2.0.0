import os
import numpy
import math
import copy

'''
Water Functions:
Developed by Mohsen Moradi and Amir A. Aliabadi
Atmospheric Innovations Research (AIR) Laboratory, University of Guelph, Guelph, Canada
Last update: May 2021
Originally developed by Naika Meili 
'''

class RadiationFunctions(object):

    def TotalLWRabsorbed(self,TemperatureC,geometry,MeteoData,FractionsGround,PropOpticalGround,PropOpticalWall,
                         PropOpticalTree,ParTree,ViewFactor):

        """
        ------
        INPUT:
        TemperatureC: Temperature of the exterior surfaces [K]
        geometry: Normalized geometric parameters
        MeteoData: Forcing variables
        FractionsGround:
        PropOpticalGround: Optical properties of the ground (albedo and emissivity)
        PropOpticalWall: Optical properties of the wall (albedo and emissivity)
        PropOpticalTree: Optical properties of the trees (albedo and emissivity)
        ParTree: Trees parameter
        ViewFactor:
        -------
        OUTPUT:
        LWRin_t: Incoming longwave radiations [W m^-2]
        LWRout_t: Outgoing longwave radiations [W m^-2]
        LWRabs_t: Absorbed longwave radiations [W m^-2]
        LWREB_t: Energy conservation of longwave radiations [W m^-2]
        """

        Tgrimp = TemperatureC[0]
        Tgbare = TemperatureC[1]
        Tgveg = TemperatureC[2]
        Twsun = TemperatureC[3]
        Twshade = TemperatureC[4]
        Ttree = TemperatureC[5]


        if ParTree.trees == 1:
            # Longwave radiation absorbed without trees
            LWRin_nT, LWRout_nT, LWRabs_nT, LWREB_nT = self.LWRabsorbedNoTree(geometry.hcanyon,geometry.wcanyon,MeteoData.LWR,
                                                                              FractionsGround.fveg,FractionsGround.fbare,FractionsGround.fimp,
                                                                              PropOpticalWall.emissivity,PropOpticalGround.eveg,
                                                                              PropOpticalGround.ebare,PropOpticalGround.eimp,Tgrimp,Tgbare,
                                                                              Tgveg,Twsun,Twshade,ViewFactor)
            # Longwave radiation absorbed with trees
            LWRin_T, LWRout_T, LWRabs_T, LWREB_T = self.LWRabsorbedWithTrees(geometry.hcanyon, geometry.wcanyon, geometry.radius_tree,
                                                                             MeteoData.LWR, FractionsGround.fveg, FractionsGround.fbare,
                                                                             FractionsGround.fimp, PropOpticalWall.emissivity,
                                                                             PropOpticalTree.emissivity, PropOpticalGround.eveg,
                                                                             PropOpticalGround.ebare, PropOpticalGround.eimp, Tgrimp,
                                                                             Tgbare, Tgveg, Twsun, Twshade, Ttree, ViewFactor)

            class LWRin_t_Def():
                pass
            LWRin_t = LWRin_t_Def()
            LWRin_t.LWRinGroundImp = ParTree.ftree*LWRin_T.LWRinGroundImp + (1-ParTree.ftree)*LWRin_nT.LWRinGroundImp
            LWRin_t.LWRinGroundBare = ParTree.ftree*LWRin_T.LWRinGroundBare + (1-ParTree.ftree)*LWRin_nT.LWRinGroundBare
            LWRin_t.LWRinGroundVeg = ParTree.ftree*LWRin_T.LWRinGroundVeg + (1-ParTree.ftree)*LWRin_nT.LWRinGroundVeg
            LWRin_t.LWRinTree = ParTree.ftree*LWRin_T.LWRinTree + (1-ParTree.ftree)*LWRin_nT.LWRinTree
            LWRin_t.LWRinWallSun = ParTree.ftree*LWRin_T.LWRinWallSun + (1-ParTree.ftree)*LWRin_nT.LWRinWallSun
            LWRin_t.LWRinWallShade = ParTree.ftree*LWRin_T.LWRinWallShade + (1-ParTree.ftree)*LWRin_nT.LWRinWallShade
            LWRin_t.LWRinTotalGround = ParTree.ftree*LWRin_T.LWRinTotalGround + (1-ParTree.ftree)*LWRin_nT.LWRinTotalGround
            LWRin_t.LWRinTotalCanyon = ParTree.ftree*LWRin_T.LWRinTotalCanyon + (1-ParTree.ftree)*LWRin_nT.LWRinTotalCanyon

            class LWRout_t_Def():
                pass
            LWRout_t = LWRout_t_Def()
            LWRout_t.LWRoutGroundImp = ParTree.ftree*LWRout_T.LWRoutGroundImp + (1-ParTree.ftree)*LWRout_nT.LWRoutGroundImp
            LWRout_t.LWRoutGroundBare = ParTree.ftree*LWRout_T.LWRoutGroundBare + (1-ParTree.ftree)*LWRout_nT.LWRoutGroundBare
            LWRout_t.LWRoutGroundVeg = ParTree.ftree*LWRout_T.LWRoutGroundVeg + (1-ParTree.ftree)*LWRout_nT.LWRoutGroundVeg
            LWRout_t.LWRoutTree = ParTree.ftree*LWRout_T.LWRoutTree + (1-ParTree.ftree)*LWRout_nT.LWRoutTree
            LWRout_t.LWRoutWallSun = ParTree.ftree*LWRout_T.LWRoutWallSun + (1-ParTree.ftree)*LWRout_nT.LWRoutWallSun
            LWRout_t.LWRoutWallShade = ParTree.ftree*LWRout_T.LWRoutWallShade + (1-ParTree.ftree)*LWRout_nT.LWRoutWallShade
            LWRout_t.LWRoutTotalGround = ParTree.ftree*LWRout_T.LWRoutTotalGround + (1-ParTree.ftree)*LWRout_nT.LWRoutTotalGround
            LWRout_t.LWRoutTotalCanyon = ParTree.ftree*LWRout_T.LWRoutTotalCanyon + (1-ParTree.ftree)*LWRout_nT.LWRoutTotalCanyon

            class LWRabs_t_Def():
                pass
            LWRabs_t = LWRabs_t_Def()
            LWRabs_t.LWRabsGroundImp = ParTree.ftree*LWRabs_T.LWRabsGroundImp + (1-ParTree.ftree)*LWRabs_nT.LWRabsGroundImp
            LWRabs_t.LWRabsGroundBare = ParTree.ftree*LWRabs_T.LWRabsGroundBare + (1-ParTree.ftree)*LWRabs_nT.LWRabsGroundBare
            LWRabs_t.LWRabsGroundVeg = ParTree.ftree*LWRabs_T.LWRabsGroundVeg + (1-ParTree.ftree)*LWRabs_nT.LWRabsGroundVeg
            LWRabs_t.LWRabsTree = ParTree.ftree*LWRabs_T.LWRabsTree + (1-ParTree.ftree)*LWRabs_nT.LWRabsTree
            LWRabs_t.LWRabsWallSun = ParTree.ftree*LWRabs_T.LWRabsWallSun + (1-ParTree.ftree)*LWRabs_nT.LWRabsWallSun
            LWRabs_t.LWRabsWallShade = ParTree.ftree*LWRabs_T.LWRabsWallShade + (1-ParTree.ftree)*LWRabs_nT.LWRabsWallShade
            LWRabs_t.LWRabsTotalGround = ParTree.ftree*LWRabs_T.LWRabsTotalGround + (1-ParTree.ftree)*LWRabs_nT.LWRabsTotalGround
            LWRabs_t.LWRabsTotalCanyon = ParTree.ftree*LWRabs_T.LWRabsTotalCanyon + (1-ParTree.ftree)*LWRabs_nT.LWRabsTotalCanyon

            class LWREB_t_Def():
                pass
            LWREB_t = LWREB_t_Def()
            LWREB_t.LWREBGroundImp = ParTree.ftree*LWREB_T.LWREBGroundImp + (1-ParTree.ftree)*LWREB_nT.LWREBGroundImp
            LWREB_t.LWREBGroundBare = ParTree.ftree*LWREB_T.LWREBGroundBare + (1-ParTree.ftree)*LWREB_nT.LWREBGroundBare
            LWREB_t.LWREBGroundVeg = ParTree.ftree*LWREB_T.LWREBGroundVeg + (1-ParTree.ftree)*LWREB_nT.LWREBGroundVeg
            LWREB_t.LWREBTree = ParTree.ftree*LWREB_T.LWREBTree + (1-ParTree.ftree)*LWREB_nT.LWREBTree
            LWREB_t.LWREBWallSun = ParTree.ftree*LWREB_T.LWREBWallSun + (1-ParTree.ftree)*LWREB_nT.LWREBWallSun
            LWREB_t.LWREBWallShade = ParTree.ftree*LWREB_T.LWREBWallShade + (1-ParTree.ftree)*LWREB_nT.LWREBWallShade
            LWREB_t.LWREBTotalGround = ParTree.ftree*LWREB_T.LWREBTotalGround + (1-ParTree.ftree)*LWREB_nT.LWREBTotalGround
            LWREB_t.LWREBTotalCanyon = ParTree.ftree*LWREB_T.LWREBTotalCanyon + (1-ParTree.ftree)*LWREB_nT.LWREBTotalCanyon

            # The absorbed radiation by the tree is not averaged as it is per tree surface
            LWRin_t.LWRinTree = LWRin_T.LWRinTree
            LWRout_t.LWRoutTree = LWRout_T.LWRoutTree
            LWRabs_t.LWRabsTree = LWRabs_T.LWRabsTree
            LWREB_t.LWREBTree = LWREB_T.LWREBTree

        elif ParTree.trees == 0:
            # Longwave radiation absorbed without trees
            LWRin_nT, LWRout_nT, LWRabs_nT, LWREB_nT = self.LWRabsorbedNoTree(geometry.hcanyon,geometry.wcanyon,MeteoData.LWR,
                                                                              FractionsGround.fveg,FractionsGround.fbare,FractionsGround.fimp,
                                                                              PropOpticalWall.emissivity,PropOpticalGround.eveg,
                                                                              PropOpticalGround.ebare,PropOpticalGround.eimp,
                                                                              Tgrimp,Tgbare,Tgveg,Twsun,Twshade,ViewFactor)
            LWRin_t = LWRin_nT
            LWRout_t = LWRout_nT
            LWRabs_t = LWRabs_nT
            LWREB_t = LWREB_nT

        return LWRin_t,LWRout_t,LWRabs_t,LWREB_t

    def TotalSWRabsorbed(self,geometry,FractionsGround,ParTree,PropOpticalGround,PropOpticalWall,PropOpticalTree,
                         ParVegTree,MeteoData,SunPosition,ViewFactor):

        """
        ------
        INPUT:
        geometry: Normalized geometric parameters
        FractionsGround:
        ParTree: Trees parameter
        PropOpticalGround: Optical properties of the ground (albedo and emissivity)
        PropOpticalWall: Optical properties of the wall (albedo and emissivity)
        PropOpticalTree: Optical properties of the trees (albedo and emissivity)
        ParVegTree: Trees parameter
        MeteoData: Forcing variables
        SunPosition: Sun angles
        ViewFactor:
        -------
        OUTPUT:
        SWRin_t: Incoming shortwave radiations [W m^-2]
        SWRout_t: Outgoing shortwave radiations [W m^-2]
        SWRabs_t: Absorbed shortwave radiations [W m^-2]
        SWRabsDir_t: Absorbed direct shortwave radiations [W m^-2]
        SWRabsDiff_t: Absorbed diffuse shortwave radiations [W m^-2]
        SWREB_t: Energy conservation of shortwave radiations [W m^-2]
        """

        if ParTree.trees == 1:
            # Shortwave radiation absorbed without trees
            SWRin_nT, SWRout_nT, SWRabs_nT, SWRabsDir_nT, SWRabsDiff_nT, SWREB_nT = \
                self.SWRabsorbedNoTrees(geometry.hcanyon,geometry.wcanyon,FractionsGround.fveg,FractionsGround.fbare,FractionsGround.fimp,
                                        PropOpticalWall.albedo,PropOpticalGround.aveg,PropOpticalGround.abare,PropOpticalGround.aimp,
                                        MeteoData.SW_dir,MeteoData.SW_diff,SunPosition.theta_Z,SunPosition.theta_n,ViewFactor,ParVegTree)

            # Shortwave radiation absorbed with trees
            SWRin_T, SWRout_T, SWRabs_T, SWRabsDir_T, SWRabsDiff_T, SWREB_T = \
                self.SWRabsorbedWithTrees(geometry.hcanyon,geometry.wcanyon,geometry.htree,geometry.radius_tree,geometry.distance_tree,
                                          FractionsGround.fveg,FractionsGround.fbare,FractionsGround.fimp,PropOpticalWall.albedo,
                                          PropOpticalGround.aveg,PropOpticalGround.abare,PropOpticalGround.aimp,PropOpticalTree.albedo,
                                          ParVegTree.LAI,MeteoData.SW_dir,MeteoData.SW_diff,SunPosition.theta_Z,SunPosition.theta_n,
                                          ViewFactor,ParVegTree)

            class SWRin_t_Def():
                pass
            SWRin_t = SWRin_t_Def()
            SWRin_t.SWRinGroundImp = ParTree.ftree * SWRin_T.SWRinGroundImp + (1 - ParTree.ftree) * SWRin_nT.SWRinGroundImp
            SWRin_t.SWRinGroundBare = ParTree.ftree * SWRin_T.SWRinGroundBare + (1 - ParTree.ftree) * SWRin_nT.SWRinGroundBare
            SWRin_t.SWRinGroundVeg = ParTree.ftree * SWRin_T.SWRinGroundVeg + (1 - ParTree.ftree) * SWRin_nT.SWRinGroundVeg
            SWRin_t.SWRinTree = ParTree.ftree * SWRin_T.SWRinTree + (1 - ParTree.ftree) * SWRin_nT.SWRinTree
            SWRin_t.SWRinWallSun = ParTree.ftree * SWRin_T.SWRinWallSun + (1 - ParTree.ftree) * SWRin_nT.SWRinWallSun
            SWRin_t.SWRinWallShade = ParTree.ftree * SWRin_T.SWRinWallShade + (1 - ParTree.ftree) * SWRin_nT.SWRinWallShade
            SWRin_t.SWRinTotalGround = ParTree.ftree * SWRin_T.SWRinTotalGround + (1 - ParTree.ftree) * SWRin_nT.SWRinTotalGround
            SWRin_t.SWRinTotalCanyon = ParTree.ftree * SWRin_T.SWRinTotalCanyon + (1 - ParTree.ftree) * SWRin_nT.SWRinTotalCanyon

            class SWRout_t_Def():
                pass
            SWRout_t = SWRout_t_Def()
            SWRout_t.SWRoutGroundImp = ParTree.ftree * SWRout_T.SWRoutGroundImp + (1 - ParTree.ftree) * SWRout_nT.SWRoutGroundImp
            SWRout_t.SWRoutGroundBare = ParTree.ftree * SWRout_T.SWRoutGroundBare + (1 - ParTree.ftree) * SWRout_nT.SWRoutGroundBare
            SWRout_t.SWRoutGroundVeg = ParTree.ftree * SWRout_T.SWRoutGroundVeg + (1 - ParTree.ftree) * SWRout_nT.SWRoutGroundVeg
            SWRout_t.SWRoutTree = ParTree.ftree * SWRout_T.SWRoutTree + (1 - ParTree.ftree) * SWRout_nT.SWRoutTree
            SWRout_t.SWRoutWallSun = ParTree.ftree * SWRout_T.SWRoutWallSun + (1 - ParTree.ftree) * SWRout_nT.SWRoutWallSun
            SWRout_t.SWRoutWallShade = ParTree.ftree * SWRout_T.SWRoutWallShade + (1 - ParTree.ftree) * SWRout_nT.SWRoutWallShade
            SWRout_t.SWRoutTotalGround = ParTree.ftree * SWRout_T.SWRoutTotalGround + (1 - ParTree.ftree) * SWRout_nT.SWRoutTotalGround
            SWRout_t.SWRoutTotalCanyon = ParTree.ftree * SWRout_T.SWRoutTotalCanyon + (1 - ParTree.ftree) * SWRout_nT.SWRoutTotalCanyon

            class SWRabs_t_Def():
                pass

            SWRabs_t = SWRabs_t_Def()
            SWRabs_t.SWRabsGroundImp = ParTree.ftree * SWRabs_T.SWRabsGroundImp + (1 - ParTree.ftree) * SWRabs_nT.SWRabsGroundImp
            SWRabs_t.SWRabsGroundBare = ParTree.ftree * SWRabs_T.SWRabsGroundBare + (1 - ParTree.ftree) * SWRabs_nT.SWRabsGroundBare
            SWRabs_t.SWRabsGroundVeg = ParTree.ftree * SWRabs_T.SWRabsGroundVeg + (1 - ParTree.ftree) * SWRabs_nT.SWRabsGroundVeg
            SWRabs_t.SWRabsTree = ParTree.ftree * SWRabs_T.SWRabsTree + (1 - ParTree.ftree) * SWRabs_nT.SWRabsTree
            SWRabs_t.SWRabsWallSun = ParTree.ftree * SWRabs_T.SWRabsWallSun + (1 - ParTree.ftree) * SWRabs_nT.SWRabsWallSun
            SWRabs_t.SWRabsWallShade = ParTree.ftree * SWRabs_T.SWRabsWallShade + (1 - ParTree.ftree) * SWRabs_nT.SWRabsWallShade
            SWRabs_t.SWRabsTotalGround = ParTree.ftree * SWRabs_T.SWRabsTotalGround + (1 - ParTree.ftree) * SWRabs_nT.SWRabsTotalGround
            SWRabs_t.SWRabsTotalCanyon = ParTree.ftree * SWRabs_T.SWRabsTotalCanyon + (1 - ParTree.ftree) * SWRabs_nT.SWRabsTotalCanyon

            class SWRabsDir_t_Def():
                pass

            SWRabsDir_t = SWRabsDir_t_Def()
            SWRabsDir_t.SWRabsGroundImp = ParTree.ftree * SWRabsDir_T.SWRabsGroundImp + (1 - ParTree.ftree) * SWRabsDir_nT.SWRabsGroundImp
            SWRabsDir_t.SWRabsGroundBare = ParTree.ftree * SWRabsDir_T.SWRabsGroundBare + (1 - ParTree.ftree) * SWRabsDir_nT.SWRabsGroundBare
            SWRabsDir_t.SWRabsGroundVeg = ParTree.ftree * SWRabsDir_T.SWRabsGroundVeg + (1 - ParTree.ftree) * SWRabsDir_nT.SWRabsGroundVeg
            SWRabsDir_t.SWRabsTree = ParTree.ftree * SWRabsDir_T.SWRabsTree + (1 - ParTree.ftree) * SWRabsDir_nT.SWRabsTree
            SWRabsDir_t.SWRabsWallSun = ParTree.ftree * SWRabsDir_T.SWRabsWallSun + (1 - ParTree.ftree) * SWRabsDir_nT.SWRabsWallSun
            SWRabsDir_t.SWRabsWallShade = ParTree.ftree * SWRabsDir_T.SWRabsWallShade + (1 - ParTree.ftree) * SWRabsDir_nT.SWRabsWallShade
            SWRabsDir_t.SWRabsTotalGround = ParTree.ftree * SWRabsDir_T.SWRabsTotalGround + (1 - ParTree.ftree) * SWRabsDir_nT.SWRabsTotalGround
            SWRabsDir_t.SWRabsTotalCanyon = ParTree.ftree * SWRabs_T.SWRabsTotalCanyon + (1 - ParTree.ftree) * SWRabsDir_nT.SWRabsTotalCanyon

            class SWRabsDiff_t_Def():
                pass

            SWRabsDiff_t = SWRabsDiff_t_Def()
            SWRabsDiff_t.SWRabsGroundImp = ParTree.ftree * SWRabsDiff_T.SWRabsGroundImp + (1 - ParTree.ftree) * SWRabsDiff_nT.SWRabsGroundImp
            SWRabsDiff_t.SWRabsGroundBare = ParTree.ftree * SWRabsDiff_T.SWRabsGroundBare + (1 - ParTree.ftree) * SWRabsDiff_nT.SWRabsGroundBare
            SWRabsDiff_t.SWRabsGroundVeg = ParTree.ftree * SWRabsDiff_T.SWRabsGroundVeg + (1 - ParTree.ftree) * SWRabsDiff_nT.SWRabsGroundVeg
            SWRabsDiff_t.SWRabsTree = ParTree.ftree * SWRabsDiff_T.SWRabsTree + (1 - ParTree.ftree) * SWRabsDiff_nT.SWRabsTree
            SWRabsDiff_t.SWRabsWallSun = ParTree.ftree * SWRabsDiff_T.SWRabsWallSun + (1 - ParTree.ftree) * SWRabsDiff_nT.SWRabsWallSun
            SWRabsDiff_t.SWRabsWallShade = ParTree.ftree * SWRabsDiff_T.SWRabsWallShade + (1 - ParTree.ftree) * SWRabsDiff_nT.SWRabsWallShade
            SWRabsDiff_t.SWRabsTotalGround = ParTree.ftree * SWRabsDiff_T.SWRabsTotalGround + (1 - ParTree.ftree) * SWRabsDiff_nT.SWRabsTotalGround
            SWRabsDiff_t.SWRabsTotalCanyon = ParTree.ftree * SWRabs_T.SWRabsTotalCanyon + (1 - ParTree.ftree) * SWRabsDiff_nT.SWRabsTotalCanyon

            class SWREB_t_Def():
                pass

            SWREB_t = SWREB_t_Def()
            SWREB_t.SWREBGroundImp = ParTree.ftree * SWREB_T.SWREBGroundImp + (1 - ParTree.ftree) * SWREB_nT.SWREBGroundImp
            SWREB_t.SWREBGroundBare = ParTree.ftree * SWREB_T.SWREBGroundBare + (1 - ParTree.ftree) * SWREB_nT.SWREBGroundBare
            SWREB_t.SWREBGroundVeg = ParTree.ftree * SWREB_T.SWREBGroundVeg + (1 - ParTree.ftree) * SWREB_nT.SWREBGroundVeg
            SWREB_t.SWREBTree = ParTree.ftree * SWREB_T.SWREBTree + (1 - ParTree.ftree) * SWREB_nT.SWREBTree
            SWREB_t.SWREBWallSun = ParTree.ftree * SWREB_T.SWREBWallSun + (1 - ParTree.ftree) * SWREB_nT.SWREBWallSun
            SWREB_t.SWREBWallShade = ParTree.ftree * SWREB_T.SWREBWallShade + (1 - ParTree.ftree) * SWREB_nT.SWREBWallShade
            SWREB_t.SWREBTotalGround = ParTree.ftree * SWREB_T.SWREBTotalGround + (1 - ParTree.ftree) * SWREB_nT.SWREBTotalGround
            SWREB_t.SWREBTotalCanyon = ParTree.ftree * SWREB_T.SWREBTotalCanyon + (1 - ParTree.ftree) * SWREB_nT.SWREBTotalCanyon

            # The absorbed radiation by the tree is not averaged as it is per tree surface
            SWRin_t.SWRinTree = SWRin_T.SWRinTree
            SWRout_t.SWRoutTree = SWRout_T.SWRoutTree
            SWRabs_t.SWRabsTree = SWRabs_T.SWRabsTree
            SWRabsDir_t.SWRabsTree = SWRabsDir_T.SWRabsTree
            SWRabsDiff_t.SWRabsTree = SWRabsDiff_T.SWRabsTree
            SWREB_t.SWREBTree = SWREB_T.SWREBTree

        elif ParTree.trees == 0:
            SWRin_nT, SWRout_nT, SWRabs_nT, SWRabsDir_nT, SWRabsDiff_nT, SWREB_nT = \
                self.SWRabsorbedNoTrees(geometry.hcanyon,geometry.wcanyon,FractionsGround.fveg,FractionsGround.fbare,FractionsGround.fimp,
                                        PropOpticalWall.albedo,PropOpticalGround.aveg,PropOpticalGround.abare,PropOpticalGround.aimp,
                                        MeteoData.SW_dir,MeteoData.SW_diff,SunPosition.theta_Z,SunPosition.theta_n,ViewFactor,ParVegTree)

            SWRin_t = SWRin_nT
            SWRout_t = SWRout_nT
            SWRabs_t = SWRabs_nT
            SWRabsDir_t = SWRabsDir_nT
            SWRabsDiff_t = SWRabsDiff_nT
            SWREB_t = SWREB_nT
        return SWRin_t,SWRout_t,SWRabs_t,SWRabsDir_t,SWRabsDiff_t,SWREB_t

    def LWRabsorbedNoTree(self,h_can,w_can,LWR,fgveg,fgbare,fgimp,ew,egveg,egbare,egimp,Tgimp,Tgbare,Tgveg,Twsun,Twshade,ViewFactor):

        F_gs_nT = ViewFactor.F_gs_nT
        F_gw_nT = ViewFactor.F_gw_nT
        F_ww_nT = ViewFactor.F_ww_nT
        F_wg_nT = ViewFactor.F_wg_nT
        F_ws_nT = ViewFactor.F_ws_nT
        F_sg_nT = ViewFactor.F_sg_nT
        F_sw_nT = ViewFactor.F_sw_nT

        # normalized surface areas
        A_s = w_can
        A_g = w_can
        A_w = h_can
        bolzm = 5.67 * 10 ** (-8)

        SVF = numpy.zeros(3)
        SVF[0] = F_gs_nT + 2 * F_gw_nT
        SVF[1] = F_ww_nT + F_wg_nT + F_ws_nT
        SVF[2] = F_sg_nT + 2 * F_sw_nT

        SVF2 = numpy.zeros(3)
        SVF2[0] = F_gs_nT + 2 * F_ws_nT * h_can
        SVF2[1] = F_sg_nT + 2 * F_wg_nT * h_can
        SVF2[2] = F_ww_nT + F_sw_nT / h_can + F_gw_nT / h_can

        for i in range(0,len(SVF)):
            if SVF[i] < 0.999 or SVF[i] > 1.001:
                print('The view factors do not add up to 1 for a canyon with trees')

        # Solve for infinite reflections equation A*X=C
        if fgimp > 0:
            Cimp = 1
        else:
            Cimp = 0
        if fgbare > 0:
            Cbare = 1
        else:
            Cbare = 0
        if fgveg > 0:
            Cveg = 1
        else:
            Cveg = 0

        #  View factor matrix to solve for infinite reflections equation
        Tij = numpy.array([[1,0,0, -(1-egveg)*F_gw_nT*Cveg, -(1-egveg)*F_gw_nT*Cveg, -(1-egveg)*F_gs_nT*Cveg],
               [0,1,0, -(1-egbare)*F_gw_nT*Cbare, -(1-egbare)*F_gw_nT*Cbare, -(1-egbare)*F_gs_nT*Cbare],
               [0,0,1, -(1-egimp)*F_gw_nT*Cimp, -(1-egimp)*F_gw_nT*Cimp, -(1-egimp)*F_gs_nT*Cimp],
               [-(1-ew)*F_wg_nT*fgveg*Cveg, -(1-ew)*F_wg_nT*fgbare*Cbare, -(1-ew)*F_wg_nT*fgimp*Cimp, 1, -(1-ew)*F_ww_nT, -(1-ew)*F_ws_nT],
               [-(1-ew)*F_wg_nT*fgveg*Cveg, -(1-ew)*F_wg_nT*fgbare*Cbare, -(1-ew)*F_wg_nT*fgimp*Cimp, -(1-ew)*F_ww_nT, 1, -(1-ew)*F_ws_nT],
               [0, 0, 0, 0, 0, 1]])

        # Emitted radiation per surface
        Omega_i = numpy.array([(egveg*bolzm*(Tgveg)**4*Cveg),(egbare * bolzm * (Tgbare) ** 4 * Cbare),(egimp * bolzm * (Tgimp) ** 4 * Cimp),(ew * bolzm * (Twsun) ** 4),(ew * bolzm * (Twshade) ** 4),LWR])

        # Outgoing radiation per surface
        # Outgoing radiation [W/m^2] per m^2 surface area

        B_i = numpy.linalg.solve(Tij,Omega_i)

        if B_i[5] != LWR:
            print('Incoming lonwave radiation and emitted longwave radiation from the sky after the matrix inversion are not equal')

        # Incoming longwave radiation at each surface A_i
        Tij2 = numpy.array([[0, 0, 0, F_gw_nT*Cveg, F_gw_nT*Cveg, F_gs_nT*Cveg],
                [0, 0, 0, F_gw_nT*Cbare, F_gw_nT*Cbare, F_gs_nT*Cbare],
                [0, 0, 0, F_gw_nT*Cimp, F_gw_nT*Cimp, F_gs_nT*Cimp],
                [F_wg_nT*fgveg*Cveg, F_wg_nT*fgbare*Cbare, F_wg_nT*fgimp*Cimp, 0, F_ww_nT, F_ws_nT],
                [F_wg_nT*fgveg*Cveg, F_wg_nT*fgbare*Cbare, F_wg_nT*fgimp*Cimp, F_ww_nT, 0, F_ws_nT],
                [0, 0, 0, 0, 0, 0]])

        A_i = numpy.dot(Tij2,B_i)
        e_i = [egveg,egbare,egimp,ew,ew,0]
        A_i2 = [(B_i[i]-Omega_i[i])/(1-e_i[i]) for i in range(0,len(e_i))]
        Qnet_i2 = A_i - B_i

        # Absorbed longwave radiation
        e_i = [egveg, egbare, egimp, ew, ew, 0]
        Qnet_i = [(e_i[i] * B_i[i] - Omega_i[i]) / (1 - e_i[i]) for i in range(0,len(e_i))]
        for i in range(0,len(e_i)):
            if e_i[i] == 1:
                Qnet_i[i] = A_i[i] - Omega_i[i]

        # Assumption: The sky has a fixed emission of LWR. Hence, Qnet is 0
        Qnet_i[5] = 0

        # Assignment
        LWRout_i = B_i	    # Outgoing radiation [W/m^2] per m^2 surface area
        LWRemit_i = Omega_i # Emitted radiation [W/m^2] per m^2 surface area
        LWRin_i = A_i  	    # Incoming radiation [W/m^2] per m^2 surface area
        LWRnet_i = Qnet_i	# Net absorbed radiation [W/m^2] per m^2 surface area

        # Energy Balance
        LWRin_atm = LWR
        TotalLWRSurface_in = LWRin_i[0]*fgveg*A_g/A_g + LWRin_i[1]*fgbare*A_g/A_g + LWRin_i[2]*fgimp*A_g/A_g +\
                             LWRin_i[3]*A_w/A_g + LWRin_i[4]*A_w/A_g

        TotalLWRSurface_abs	=	LWRnet_i[0]*fgveg*A_g/A_g + LWRnet_i[1]*fgbare*A_g/A_g + LWRnet_i[2]*fgimp*A_g/A_g +\
                                LWRnet_i[3]*A_w/A_g + LWRnet_i[4]*A_w/A_g

        TotalLWRSurface_out	=	LWRout_i[0]*fgveg*A_g/A_s+LWRout_i[1]*fgbare*A_g/A_s+LWRout_i[2]*fgimp*A_g/A_s +\
                                LWRout_i[3]*A_w/A_s+LWRout_i[4]*A_w/A_s



        TotalLWRref_to_atm	=	LWRout_i[0]*F_sg_nT*fgveg + LWRout_i[1]*F_sg_nT*fgbare + LWRout_i[2]*F_sg_nT*fgimp + \
                                LWRout_i[3]*F_sw_nT + LWRout_i[4]*F_sw_nT

        EBSurface = TotalLWRSurface_in - TotalLWRSurface_abs - TotalLWRSurface_out
        EBCanyon = LWRin_atm - TotalLWRSurface_abs - TotalLWRref_to_atm

        # Energy balance
        if abs(EBSurface)>=10**(-1):
            print('EBSurface is not 0',EBSurface)
        if abs(EBCanyon)>=10**(-1):
            print('EBCanyon is not 0',EBCanyon)

        class LWRin_nT_Def():
            pass
        LWRin_nT = LWRin_nT_Def()
        LWRin_nT.LWRinGroundImp = LWRin_i[2] * Cimp
        LWRin_nT.LWRinGroundBare = LWRin_i[1] * Cbare
        LWRin_nT.LWRinGroundVeg = LWRin_i[0] * Cveg
        LWRin_nT.LWRinTree = 0
        LWRin_nT.LWRinWallSun = LWRin_i[3]
        LWRin_nT.LWRinWallShade = LWRin_i[4]
        LWRin_nT.LWRinTotalGround = fgveg * LWRin_i[0] + fgbare * LWRin_i[1] + fgimp * LWRin_i[2]
        LWRin_nT.LWRinTotalCanyon = LWRin_i[0] * fgveg * A_g / A_g + LWRin_i[1] * fgbare * A_g / A_g + \
                                    LWRin_i[2] * fgimp * A_g / A_g + LWRin_i[3] * A_w / A_g + LWRin_i[4] * A_w / A_g

        # Outgoing longwave radiation
        class LWRout_nT_Def():
            pass
        LWRout_nT = LWRout_nT_Def()
        LWRout_nT.LWRoutGroundImp = LWRout_i[2] * Cimp
        LWRout_nT.LWRoutGroundBare = LWRout_i[1] * Cbare
        LWRout_nT.LWRoutGroundVeg = LWRout_i[0] * Cveg
        LWRout_nT.LWRoutTree = 0
        LWRout_nT.LWRoutWallSun = LWRout_i[3]
        LWRout_nT.LWRoutWallShade = LWRout_i[4]
        LWRout_nT.LWRoutTotalGround = fgveg * LWRout_i[0] + fgbare * LWRout_i[1] + fgimp * LWRout_i[2]
        LWRout_nT.LWRoutTotalCanyon = LWRout_i[0] * fgveg * A_g / A_g + LWRout_i[1] * fgbare * A_g / A_g + \
                                      LWRout_i[2] * fgimp * A_g / A_g + LWRout_i[3] * A_w / A_g + LWRout_i[4] * A_w / A_g

        # Absorbed longwave radiation
        class LWRabs_nT_Def():
            pass
        LWRabs_nT = LWRabs_nT_Def()
        LWRabs_nT.LWRabsGroundImp = LWRnet_i[2] * Cimp
        LWRabs_nT.LWRabsGroundBare = LWRnet_i[1] * Cbare
        LWRabs_nT.LWRabsGroundVeg = LWRnet_i[0] * Cveg
        LWRabs_nT.LWRabsTree = 0
        LWRabs_nT.LWRabsWallSun = LWRnet_i[3]
        LWRabs_nT.LWRabsWallShade = LWRnet_i[4]
        LWRabs_nT.LWRabsTotalGround = fgveg * LWRnet_i[0] + fgbare * LWRnet_i[1] + fgimp * LWRnet_i[2]
        LWRabs_nT.LWRabsTotalCanyon = LWRnet_i[0] * fgveg * A_g / A_g + LWRnet_i[1] * fgbare * A_g / A_g + \
                                      LWRnet_i[2] * fgimp * A_g / A_g + LWRnet_i[3] * A_w / A_g + LWRnet_i[4] * A_w / A_g

        # Energy Balance of longwave radiation
        class LWREB_nT_Def():
            pass
        LWREB_nT = LWREB_nT_Def()
        LWREB_nT.LWREBGroundImp = LWRin_nT.LWRinGroundImp - LWRout_nT.LWRoutGroundImp - LWRabs_nT.LWRabsGroundImp
        LWREB_nT.LWREBGroundBare = LWRin_nT.LWRinGroundBare - LWRout_nT.LWRoutGroundBare - LWRabs_nT.LWRabsGroundBare
        LWREB_nT.LWREBGroundVeg = LWRin_nT.LWRinGroundVeg - LWRout_nT.LWRoutGroundVeg - LWRabs_nT.LWRabsGroundVeg
        LWREB_nT.LWREBTree = 0
        LWREB_nT.LWREBWallSun = LWRin_nT.LWRinWallSun - LWRout_nT.LWRoutWallSun - LWRabs_nT.LWRabsWallSun
        LWREB_nT.LWREBWallShade = LWRin_nT.LWRinWallShade - LWRout_nT.LWRoutWallShade - LWRabs_nT.LWRabsWallShade
        LWREB_nT.LWREBTotalGround = LWRin_nT.LWRinTotalGround - LWRout_nT.LWRoutTotalGround - LWRabs_nT.LWRabsTotalGround
        LWREB_nT.LWREBTotalCanyon = LWRin_nT.LWRinTotalCanyon - LWRout_nT.LWRoutTotalCanyon - LWRabs_nT.LWRabsTotalCanyon

        if abs(LWREB_nT.LWREBGroundImp)>=10**(-6):
            print('LWREB_nT.LWREBGroundImp is not 0.')
        if abs(LWREB_nT.LWREBGroundBare)>=10**(-6):
            print('LWREB_nT.LWREBGroundBare is not 0.')
        if abs(LWREB_nT.LWREBGroundVeg	)>=10**(-6):
            print('LWREB_nT.LWREBGroundVeg	 is not 0.')
        if abs(LWREB_nT.LWREBWallSun)>=10**(-6):
            print('LWREB_nT.LWREBWallSun is not 0.')
        if abs(LWREB_nT.LWREBWallShade)>=10**(-6):
            print('LWREB_nT.LWREBWallShade is not 0.')
        if abs(LWREB_nT.LWREBTotalGround)>=10**(-6):
            print('LWREB_nT.LWREBTotalGround is not 0.')
        if abs(LWREB_nT.LWREBTotalCanyon)>=10**(-6):
            print('LWREB_nT.LWREBTotalCanyon is not 0.')

        return LWRin_nT,LWRout_nT,LWRabs_nT,LWREB_nT

    def LWRabsorbedWithTrees(self,h_can,w_can,r_tree,LWR,fgveg,fgbare,fgimp,ew,et,egveg,egbare,egimp,Tgimp,Tgbare,Tgveg,
                             Twsun,Twshade,Ttree,ViewFactor):
        F_gs_T = ViewFactor.F_gs_T
        F_gt_T = ViewFactor.F_gt_T
        F_gw_T = ViewFactor.F_gw_T
        F_ww_T = ViewFactor.F_ww_T
        F_wt_T = ViewFactor.F_wt_T
        F_wg_T = ViewFactor.F_wg_T
        F_ws_T = ViewFactor.F_ws_T
        F_sg_T = ViewFactor.F_sg_T
        F_sw_T = ViewFactor.F_sw_T
        F_st_T = ViewFactor.F_st_T
        F_tg_T = ViewFactor.F_tg_T
        F_tw_T = ViewFactor.F_tw_T
        F_ts_T = ViewFactor.F_ts_T
        F_tt_T = ViewFactor.F_tt_T

        # normalized surface areas
        A_s = w_can
        A_g = w_can
        A_w = h_can
        # There are 2 trees. Hence, the area of tree is twice a circle
        A_t = 2 * 2 * numpy.pi * r_tree
        bolzm = 5.67 * 10 ** (-8)

        # Check if view factors add up to 1
        SVF = numpy.zeros(4)
        SVF[0] = F_gs_T + F_gt_T + 2 * F_gw_T
        SVF[1] = F_ww_T + F_wt_T + F_wg_T + F_ws_T
        SVF[2] = F_sg_T + 2 * F_sw_T + F_st_T
        SVF[3] = F_ts_T + 2 * F_tw_T + F_tt_T + F_tg_T

        SVF2 = numpy.zeros(4)
        SVF2[0] = F_gs_T + 2 * F_ws_T + F_ts_T
        SVF2[1] = F_sg_T + 2 * F_wg_T + F_tg_T
        SVF2[2] = F_ww_T + F_sw_T + F_gw_T + F_tw_T
        SVF2[3] = F_gt_T + 2 * F_wt_T + F_tt_T + F_st_T

        for i in range(0,len(SVF)):
            if SVF[i]<0.999 or SVF[i]>1.001:
                print('The view factors do not add up to 1 for a canyon with trees')

        # Solve for infinite reflections equation A*X=C
        if fgimp > 0:
            Cimp = 1
        else:
            Cimp = 0
        if fgbare > 0:
            Cbare = 1
        else:
            Cbare = 0
        if fgveg > 0:
            Cveg = 1
        else:
            Cveg = 0

        Tij = numpy.array([[1,0,0, -(1-egveg)*F_gw_T*Cveg, -(1-egveg)*F_gw_T*Cveg, -(1-egveg)*F_gt_T*Cveg, -(1-egveg)*F_gs_T*Cveg],
               [0,1,0, -(1-egbare)*F_gw_T*Cbare, -(1-egbare)*F_gw_T*Cbare, -(1-egbare)*F_gt_T*Cbare, -(1-egbare)*F_gs_T*Cbare],
               [0,0,1, -(1-egimp)*F_gw_T*Cimp, -(1-egimp)*F_gw_T*Cimp, -(1-egimp)*F_gt_T*Cimp, -(1-egimp)*F_gs_T*Cimp],
               [-(1-ew)*F_wg_T*fgveg*Cveg, -(1-ew)*F_wg_T*fgbare*Cbare, -(1-ew)*F_wg_T*fgimp*Cimp, 1, -(1-ew)*F_ww_T, -(1-ew)*F_wt_T, -(1-ew)*F_ws_T],
               [-(1-ew)*F_wg_T*fgveg*Cveg, -(1-ew)*F_wg_T*fgbare*Cbare, -(1-ew)*F_wg_T*fgimp*Cimp, -(1-ew)*F_ww_T, 1, -(1-ew)*F_wt_T, -(1-ew)*F_ws_T],
               [-(1-et)*F_tg_T*fgveg*Cveg, -(1-et)*F_tg_T*fgbare*Cbare, -(1-et)*F_tg_T*fgimp*Cimp, -(1-et)*F_tw_T, -(1-et)*F_tw_T, 1-(1-et)*F_tt_T, -(1-et)*F_ts_T],
               [0, 0, 0, 0, 0, 0, 1]])

        Omega_i = numpy.array([(egveg*bolzm*(Tgveg)**4*Cveg),
                   (egbare * bolzm * (Tgbare) ** 4 * Cbare),
                   (egimp * bolzm * (Tgimp) ** 4 * Cimp),
                   (ew * bolzm * (Twsun) ** 4),
                   (ew * bolzm * (Twshade) ** 4),
                   (et * bolzm * (Ttree) ** 4),
                   LWR])

        # Outgoing radiation per surface
        # Outgoing radiation [W/m^2] per m^2 surface area
        B_i = numpy.linalg.solve(Tij,Omega_i)

        if B_i[6] != LWR:
            print('Incoming lonwave radiation and emitted longwave radiation from the sky after the matrix inversion are not equal')

        Tij2 = numpy.array([[0, 0, 0, F_gw_T*Cveg, F_gw_T*Cveg, F_gt_T*Cveg, F_gs_T*Cveg],
                [0, 0, 0, F_gw_T*Cbare, F_gw_T*Cbare, F_gt_T*Cbare, F_gs_T*Cbare],
                [0, 0, 0, F_gw_T*Cimp, F_gw_T*Cimp, F_gt_T*Cimp, F_gs_T*Cimp],
                [F_wg_T*fgveg*Cveg, F_wg_T*fgbare*Cbare, F_wg_T*fgimp*Cimp, 0, F_ww_T, F_wt_T, F_ws_T],
                [F_wg_T*fgveg*Cveg, F_wg_T*fgbare*Cbare, F_wg_T*fgimp*Cimp, F_ww_T, 0, F_wt_T, F_ws_T],
                [F_tg_T*fgveg*Cveg, F_tg_T*fgbare*Cbare, F_tg_T*fgimp*Cimp, F_tw_T, F_tw_T, F_tt_T, F_ts_T],
                [0, 0, 0, 0, 0, 0, 0]])

        A_i = numpy.dot(Tij2, B_i)
        e_i = [egveg, egbare, egimp, ew, ew, et, 0]
        A_i2 = [(B_i[i] - Omega_i[i]) / (1 - e_i[i]) for i in range(0,len(e_i))]
        Qnet_i2 = A_i - B_i

        # Absorbed longwave radiation
        e_i = [egveg, egbare, egimp, ew, ew, et, 0]
        Qnet_i = [(e_i[i] * B_i[i] - Omega_i[i]) / (1 - e_i[i]) for i in range(0,len(e_i))]
        Qnet_i = [(e_i[i] * B_i[i] - Omega_i[i]) / (1 - e_i[i]) for i in range(0, len(e_i))]
        for i in range(0, len(e_i)):
            if e_i[i] == 1:
                Qnet_i[i] = A_i[i] - Omega_i[i]

        # Assumption: The sky has a fixed emission of LWR. Hence, Qnet is 0
        Qnet_i[6] = 0

        # Assignment
        LWRout_i = B_i       # Outgoing radiation [W/m^2] per m^2 surface area
        LWRemit_i = Omega_i  # Emitted radiation [W/m^2] per m^2 surface area
        LWRin_i = A_i        # Incoming radiation [W/m^2] per m^2 surface area
        LWRnet_i = Qnet_i    # Net absorbed radiation [W/m^2] per m^2 surface area

        # Energy balance
        LWRin_atm = LWR
        TotalLWRSurface_in = LWRin_i[0]*fgveg*A_g/A_g + LWRin_i[1]*fgbare*A_g/A_g + LWRin_i[2]*fgimp*A_g/A_g +\
                             LWRin_i[3]*A_w/A_g + LWRin_i[4]*A_w/A_g + LWRin_i[5]*A_t/A_g

        TotalLWRSurface_abs	=	LWRnet_i[0]*fgveg*A_g/A_g + LWRnet_i[1]*fgbare*A_g/A_g + LWRnet_i[2]*fgimp*A_g/A_g +\
                                LWRnet_i[3]*A_w/A_g + LWRnet_i[4]*A_w/A_g + LWRnet_i[5]*A_t/A_g

        TotalLWRSurface_out	=	LWRout_i[0]*fgveg*A_g/A_s+LWRout_i[1]*fgbare*A_g/A_s+LWRout_i[2]*fgimp*A_g/A_s +\
                                LWRout_i[3]*A_w/A_s+LWRout_i[4]*A_w/A_s+LWRout_i[5]*A_t/A_s

        TotalLWRref_to_atm	=	LWRout_i[0]*F_sg_T*fgveg + LWRout_i[1]*F_sg_T*fgbare + LWRout_i[2]*F_sg_T*fgimp + \
                                LWRout_i[3]*F_sw_T + LWRout_i[4]*F_sw_T + LWRout_i[5]*F_st_T

        EBSurface = TotalLWRSurface_in - TotalLWRSurface_abs - TotalLWRSurface_out
        EBCanyon = LWRin_atm - TotalLWRSurface_abs - TotalLWRref_to_atm

        # Energy balance
        if abs(EBSurface)>=10**(-1):
            print('EBSurface is not 0',EBSurface)
        if abs(EBCanyon)>=10**(-1):
            print('EBCanyon is not 0',EBCanyon)


        class LWRin_T_Def():
            pass
        LWRin_T = LWRin_T_Def()
        LWRin_T.LWRinGroundImp = LWRin_i[2] * Cimp
        LWRin_T.LWRinGroundBare = LWRin_i[1] * Cbare
        LWRin_T.LWRinGroundVeg = LWRin_i[0] * Cveg
        LWRin_T.LWRinTree = LWRin_i[5]
        LWRin_T.LWRinWallSun = LWRin_i[3]
        LWRin_T.LWRinWallShade = LWRin_i[4]
        LWRin_T.LWRinTotalGround = fgveg * LWRin_i[0] + fgbare * LWRin_i[1] + fgimp * LWRin_i[2]
        LWRin_T.LWRinTotalCanyon = LWRin_i[0]*fgveg*A_g/A_g + LWRin_i[1]*fgbare*A_g/A_g + \
                                   LWRin_i[2]*fgimp*A_g/A_g + LWRin_i[3]*A_w/A_g + LWRin_i[4]*A_w/A_g + LWRin_i[5]*A_t/A_g

        # Outgoing longwave radiation
        class LWRout_T_Def():
            pass
        LWRout_T = LWRout_T_Def()
        LWRout_T.LWRoutGroundImp = LWRout_i[2] * Cimp
        LWRout_T.LWRoutGroundBare = LWRout_i[1] * Cbare
        LWRout_T.LWRoutGroundVeg = LWRout_i[0] * Cveg
        LWRout_T.LWRoutTree = LWRout_i[5]
        LWRout_T.LWRoutWallSun = LWRout_i[3]
        LWRout_T.LWRoutWallShade = LWRout_i[4]
        LWRout_T.LWRoutTotalGround = fgveg * LWRout_i[0] + fgbare * LWRout_i[1] + fgimp * LWRout_i[2]
        LWRout_T.LWRoutTotalCanyon = LWRout_i[0] * fgveg * A_g / A_g + LWRout_i[1] * fgbare * A_g / A_g + \
                                      LWRout_i[2] * fgimp * A_g / A_g + LWRout_i[3] * A_w / A_g + LWRout_i[4] * A_w / A_g + LWRout_i[5]*A_t/A_g

        # Absorbed longwave radiation
        class LWRabs_T_Def():
            pass
        LWRabs_T = LWRabs_T_Def()
        LWRabs_T.LWRabsGroundImp = LWRnet_i[2] * Cimp
        LWRabs_T.LWRabsGroundBare = LWRnet_i[1] * Cbare
        LWRabs_T.LWRabsGroundVeg = LWRnet_i[0] * Cveg
        LWRabs_T.LWRabsTree = LWRnet_i[5]
        LWRabs_T.LWRabsWallSun = LWRnet_i[3]
        LWRabs_T.LWRabsWallShade = LWRnet_i[4]
        LWRabs_T.LWRabsTotalGround = fgveg * LWRnet_i[0] + fgbare * LWRnet_i[1] + fgimp * LWRnet_i[2]
        LWRabs_T.LWRabsTotalCanyon = LWRnet_i[0] * fgveg * A_g / A_g + LWRnet_i[1] * fgbare * A_g / A_g + \
                                      LWRnet_i[2] * fgimp * A_g / A_g + LWRnet_i[3] * A_w / A_g + LWRnet_i[4] * A_w / A_g + LWRnet_i[5]*A_t/A_g

        # Energy Balance of longwave radiation
        class LWREB_T_Def():
            pass
        LWREB_T = LWREB_T_Def()
        LWREB_T.LWREBGroundImp = LWRin_T.LWRinGroundImp - LWRout_T.LWRoutGroundImp - LWRabs_T.LWRabsGroundImp
        LWREB_T.LWREBGroundBare = LWRin_T.LWRinGroundBare - LWRout_T.LWRoutGroundBare - LWRabs_T.LWRabsGroundBare
        LWREB_T.LWREBGroundVeg = LWRin_T.LWRinGroundVeg - LWRout_T.LWRoutGroundVeg - LWRabs_T.LWRabsGroundVeg
        LWREB_T.LWREBTree = LWRin_T.LWRinTree - LWRout_T.LWRoutTree - LWRabs_T.LWRabsTree
        LWREB_T.LWREBWallSun = LWRin_T.LWRinWallSun - LWRout_T.LWRoutWallSun - LWRabs_T.LWRabsWallSun
        LWREB_T.LWREBWallShade = LWRin_T.LWRinWallShade - LWRout_T.LWRoutWallShade - LWRabs_T.LWRabsWallShade
        LWREB_T.LWREBTotalGround = LWRin_T.LWRinTotalGround - LWRout_T.LWRoutTotalGround - LWRabs_T.LWRabsTotalGround
        LWREB_T.LWREBTotalCanyon = LWRin_T.LWRinTotalCanyon - LWRout_T.LWRoutTotalCanyon - LWRabs_T.LWRabsTotalCanyon

        if abs(LWREB_T.LWREBGroundImp)>=10**(-6):
            print('LWREB_T.LWREBGroundImp is not 0.')
        if abs(LWREB_T.LWREBGroundBare)>=10**(-6):
            print('LWREB_T.LWREBGroundBare is not 0.')
        if abs(LWREB_T.LWREBGroundVeg	)>=10**(-6):
            print('LWREB_T.LWREBGroundVeg	 is not 0.')
        if abs(LWREB_T.LWREBWallSun)>=10**(-6):
            print('LWREB_T.LWREBWallSun is not 0.')
        if abs(LWREB_T.LWREBWallShade)>=10**(-6):
            print('LWREB_T.LWREBWallShade is not 0.')
        if abs(LWREB_T.LWREBTotalGround)>=10**(-6):
            print('LWREB_T.LWREBTotalGround is not 0.')
        if abs(LWREB_T.LWREBTotalCanyon)>=10**(-6):
            print('LWREB_T.LWREBTotalCanyon is not 0.')
        if abs(LWREB_T.LWREBTree)>=10**(-6):
            print('LWREB_T.LWREBTree is not 0.')

        return LWRin_T,LWRout_T,LWRabs_T,LWREB_T


    def SWRabsorbedNoTrees(self,h_can,w_can,fgveg,fgbare,fgimp,aw,agveg,agbare,agimp,SWR_dir,SWR_diff,theta_Z,theta_n,
                           ViewFactor,ParVegTree):
        F_gs_nT = ViewFactor.F_gs_nT
        F_gw_nT = ViewFactor.F_gw_nT
        F_ww_nT = ViewFactor.F_ww_nT
        F_wg_nT = ViewFactor.F_wg_nT
        F_ws_nT = ViewFactor.F_ws_nT
        F_sg_nT = ViewFactor.F_sg_nT
        F_sw_nT = ViewFactor.F_sw_nT

        # normalized surface areas
        A_s = w_can
        A_g = w_can
        A_w = h_can

        SVF = numpy.zeros(3)
        SVF[0] = F_gs_nT + 2 * F_gw_nT
        SVF[1] = F_ww_nT + F_wg_nT + F_ws_nT
        SVF[2] = F_sg_nT + 2 * F_sw_nT

        SVF2 = numpy.zeros(3)
        SVF2[0] = F_gs_nT + 2 * F_ws_nT * h_can
        SVF2[1] = F_sg_nT + 2 * F_wg_nT * h_can
        SVF2[2] = F_ww_nT + F_sw_nT / h_can + F_gw_nT / h_can

        for i in range(0, len(SVF)):
            if SVF[i] < 0.999 or SVF[i] > 1.001:
                print('The view factors do not add up to 1 for a canyon with trees')

        SWRdir_ground, SWRdir_wallsun, SWRdir_wallshade, NOTUSED = self.DirectSWRSurfaces(h_can, w_can, numpy.nan, numpy.nan, numpy.nan,
                                                                                           theta_Z,theta_n, SWR_dir, numpy.nan,
                                                                                          0,ParVegTree)

        # Balance direct shortwave radiation in
        EB_SWRdir = SWR_dir - (SWRdir_ground * A_g / A_g + SWRdir_wallsun * A_w / A_g + SWRdir_wallshade * A_w / A_g)
        EB_SWRdiff = SWR_diff - (F_sg_nT * SWR_diff + F_sw_nT * SWR_diff + F_sw_nT * SWR_diff)

        if abs(EB_SWRdir)>=10**(-1):
            print('EB_SWRdir is not 0',EB_SWRdir)
        if abs(EB_SWRdiff)>=10**(-1):
            print('EB_SWRdiff is not 0',EB_SWRdiff)

        if fgimp > 0:
            Cimp = 1
        else:
            Cimp = 0
        if fgbare > 0:
            Cbare = 1
        else:
            Cbare = 0
        if fgveg > 0:
            Cveg = 1
        else:
            Cveg = 0

        ai = [agveg,agbare,agimp,aw,aw,0]

        # View factor matrix to solve for infinite reflections equation
        Tij = numpy.array([[1,0,0, -agveg*F_gw_nT*Cveg, -agveg*F_gw_nT*Cveg, -agveg*F_gs_nT*Cveg],
               [0,1,0, -agbare*F_gw_nT*Cbare, -agbare*F_gw_nT*Cbare, -agbare*F_gs_nT*Cbare],
               [0,0,1, -agimp*F_gw_nT*Cimp, -agimp*F_gw_nT*Cimp, -agimp*F_gs_nT*Cimp],
               [-aw*F_wg_nT*fgveg*Cveg,-aw*F_wg_nT*fgbare*Cbare,-aw*F_wg_nT*fgimp*Cimp, 1, -aw*F_ww_nT, -aw*F_ws_nT],
               [-aw*F_wg_nT*fgveg*Cveg,-aw*F_wg_nT*fgbare*Cbare,-aw*F_wg_nT*fgimp*Cimp, -aw*F_ww_nT, 1, -aw*F_ws_nT],
               [0, 0, 0, 0, 0, 1]])

        # Incoming shortwave radiation from sky
        Omega_i = numpy.array([agveg*SWRdir_ground*Cveg,
                   agbare * SWRdir_ground * Cbare,
                   agimp * SWRdir_ground * Cimp,
                   aw * SWRdir_wallsun,
                   aw * 0,
                   SWR_diff])

        # Outgoing radiation per surface
        B_i = numpy.linalg.solve(Tij,Omega_i)

        if B_i[5] != SWR_diff:
            print('Incoming lonwave radiation and emitted longwave radiation from the sky after the matrix inversion are not equal')

        # Incoming shortwave radiation at each surface A_i
        Tij2 = numpy.array([[0, 0, 0, F_gw_nT*Cveg, F_gw_nT*Cveg, F_gs_nT*Cveg],
                [0, 0, 0, F_gw_nT*Cbare, F_gw_nT*Cbare, F_gs_nT*Cbare],
                [0, 0, 0, F_gw_nT*Cimp, F_gw_nT*Cimp, F_gs_nT*Cimp],
                [F_wg_nT*fgveg*Cveg, F_wg_nT*fgbare*Cbare, F_wg_nT*fgimp*Cimp, 0, F_ww_nT, F_ws_nT],
                [F_wg_nT*fgveg*Cveg, F_wg_nT*fgbare*Cbare, F_wg_nT*fgimp*Cimp, F_ww_nT, 0, F_ws_nT],
                [0, 0, 0, 0, 0, 0]])

        SWRdir_i = numpy.array([SWRdir_ground*Cveg,
                    SWRdir_ground * Cbare,
                    SWRdir_ground * Cimp,
                    SWRdir_wallsun,
                    0,
                    0])

        A_i1	=	numpy.dot(Tij2,B_i)+SWRdir_i	# Incoming radiation [W/m^2] per m^2 surface area
        A_i			=	B_i/ai		    # Incoming radiation [W/m^2] per m^2 surface area
        for i in range(0,len(ai)):
            if ai[i] == 0:
                A_i[i] = A_i1[i]
        A_i[5]		=	0				# Assumption: The sky has a fixed emission of LWR. Hence, Qnet is 0.
        # Absorbed shortwave radiation at ech surface Qnet_i
        Qnet_i		=	A_i-B_i


        # Assignment
        SWRout_i = B_i	    # Outgoing radiation [W/m^2] per m^2 surface area
        SWRin_i	= A_i	    # Incoming radiation [W/m^2] per m^2 surface area
        SWRnet_i = Qnet_i	# Net absorbed radiation [W/m^2] per m^2 surface area

        # Energy balance
        SWRin_atm = SWR_dir + SWR_diff

        TotalSWRSurface_in = SWRin_i[0] * fgveg * A_g / A_g + SWRin_i[1] * fgbare * A_g / A_g + SWRin_i[2] * fgimp * A_g / A_g +\
                             SWRin_i[3] * A_w / A_g + SWRin_i[4] * A_w / A_g

        TotalSWRSurface_abs = SWRnet_i[0] * fgveg * A_g / A_g + SWRnet_i[1] * fgbare * A_g / A_g + SWRnet_i[2] * fgimp * A_g / A_g +\
                              SWRnet_i[3] * A_w / A_g + SWRnet_i[4] * A_w / A_g

        TotalSWRSurface_out = SWRout_i[0] * fgveg * A_g / A_s + SWRout_i[1] * fgbare * A_g / A_s + SWRout_i[2] * fgimp * A_g / A_s + \
                              SWRout_i[3] * A_w / A_s + SWRout_i[4] * A_w / A_s

        TotalSWRref_to_atm = SWRout_i[0] * F_sg_nT * fgveg + SWRout_i[1] * F_sg_nT * fgbare + SWRout_i[2] * F_sg_nT * fgimp +\
                             SWRout_i[3] * F_sw_nT + SWRout_i[4] * F_sw_nT

        EBSurface = TotalSWRSurface_in - TotalSWRSurface_abs - TotalSWRSurface_out
        EBCanyon = SWRin_atm - TotalSWRSurface_abs - TotalSWRref_to_atm

        if abs(EBSurface)>=10**(-1):
            print('EBSurface is not 0',EBSurface)
        if abs(EBCanyon)>=10**(-1):
            print('EBCanyon is not 0',EBCanyon)

        # Incoming shortwave radiation
        class SWRin_nT_Def():
            pass
        SWRin_nT = SWRin_nT_Def()
        SWRin_nT.SWRinGroundImp = SWRin_i[2] * Cimp
        SWRin_nT.SWRinGroundBare = SWRin_i[1] * Cbare
        SWRin_nT.SWRinGroundVeg = SWRin_i[0] * Cveg
        SWRin_nT.SWRinTree = 0
        SWRin_nT.SWRinWallSun = SWRin_i[3]
        SWRin_nT.SWRinWallShade = SWRin_i[4]
        SWRin_nT.SWRinTotalGround = fgveg * SWRin_i[0] + fgbare * SWRin_i[1] + fgimp * SWRin_i[2]
        SWRin_nT.SWRinTotalCanyon = SWRin_i[0] * fgveg * A_g / A_g + SWRin_i[1] * fgbare * A_g / A_g + SWRin_i[2] * fgimp * A_g / A_g + \
                                    SWRin_i[3] * A_w / A_g + SWRin_i[4] * A_w / A_g

        # Outgoing shortwave radiation
        class SWRout_nT_Def():
            pass
        SWRout_nT = SWRout_nT_Def()
        SWRout_nT.SWRoutGroundImp = SWRout_i[2] * Cimp
        SWRout_nT.SWRoutGroundBare = SWRout_i[1] * Cbare
        SWRout_nT.SWRoutGroundVeg = SWRout_i[0] * Cveg
        SWRout_nT.SWRoutTree = 0
        SWRout_nT.SWRoutWallSun = SWRout_i[3]
        SWRout_nT.SWRoutWallShade = SWRout_i[4]
        SWRout_nT.SWRoutTotalGround = fgveg * SWRout_i[0] + fgbare * SWRout_i[1] + fgimp * SWRout_i[2]
        SWRout_nT.SWRoutTotalCanyon = SWRout_i[0] * fgveg * A_g / A_g + SWRout_i[1] * fgbare * A_g / A_g + SWRout_i[2] * fgimp * A_g / A_g + \
                                      SWRout_i[3] * A_w / A_g + SWRout_i[4] * A_w / A_g

        # Absorbed shortwave radiation
        class SWRabs_nT_Def():
            pass
        SWRabs_nT = SWRabs_nT_Def()
        SWRabs_nT.SWRabsGroundImp = SWRnet_i[2] * Cimp
        SWRabs_nT.SWRabsGroundBare = SWRnet_i[1] * Cbare
        SWRabs_nT.SWRabsGroundVeg = SWRnet_i[0] * Cveg
        SWRabs_nT.SWRabsTree = 0
        SWRabs_nT.SWRabsWallSun = SWRnet_i[3]
        SWRabs_nT.SWRabsWallShade = SWRnet_i[4]
        SWRabs_nT.SWRabsTotalGround = fgveg * SWRnet_i[0] + fgbare * SWRnet_i[1] + fgimp * SWRnet_i[2]
        SWRabs_nT.SWRabsTotalCanyon = SWRnet_i[0] * fgveg * A_g / A_g + SWRnet_i[1] * fgbare * A_g / A_g + SWRnet_i[2] * fgimp * A_g / A_g + \
                                      SWRnet_i[3] * A_w / A_g + SWRnet_i[4] * A_w / A_g

        # Direct absorbed shortwave radiation
        class SWRabsDir_nT_Def():
            pass
        SWRabsDir_nT = SWRabsDir_nT_Def()
        SWRabsDir_nT.SWRabsGroundImp = (1-agimp)*SWRdir_ground*Cimp
        SWRabsDir_nT.SWRabsGroundBare = (1-agbare)*SWRdir_ground*Cbare
        SWRabsDir_nT.SWRabsGroundVeg = (1-agveg)*SWRdir_ground*Cveg
        SWRabsDir_nT.SWRabsTree = 0
        SWRabsDir_nT.SWRabsWallSun = (1-aw)*SWRdir_wallsun
        SWRabsDir_nT.SWRabsWallShade = (1-aw)*SWRdir_wallshade
        SWRabsDir_nT.SWRabsTotalGround = fgveg*(1-agveg)*SWRdir_ground+fgbare*(1-agbare)*SWRdir_ground+fgimp*(1-agimp)*SWRdir_ground
        SWRabsDir_nT.SWRabsTotalCanyon = fgveg*(1-agveg)*SWRdir_ground*A_g/A_g+fgbare*(1-agbare)*SWRdir_ground*A_g/A_g+\
                                         fgimp*(1-agimp)*SWRdir_ground*A_g/A_g + (1-aw)*SWRdir_wallsun*A_w/A_g + (1-aw)*SWRdir_wallshade*A_w/A_g

        # Diffuse absorbed shortwave radiation
        class SWRabsDiff_nT_Def():
            pass
        SWRabsDiff_nT = SWRabsDiff_nT_Def()
        SWRabsDiff_nT.SWRabsGroundImp = (SWRabs_nT.SWRabsGroundImp-SWRabsDir_nT.SWRabsGroundImp)*Cimp
        SWRabsDiff_nT.SWRabsGroundBare = (SWRabs_nT.SWRabsGroundBare-SWRabsDir_nT.SWRabsGroundBare)*Cbare
        SWRabsDiff_nT.SWRabsGroundVeg = (SWRabs_nT.SWRabsGroundVeg-SWRabsDir_nT.SWRabsGroundVeg)*Cveg
        SWRabsDiff_nT.SWRabsTree = 0
        SWRabsDiff_nT.SWRabsWallSun = SWRabs_nT.SWRabsWallSun-SWRabsDir_nT.SWRabsWallSun
        SWRabsDiff_nT.SWRabsWallShade = SWRabs_nT.SWRabsWallShade-SWRabsDir_nT.SWRabsWallShade
        SWRabsDiff_nT.SWRabsTotalGround = SWRabs_nT.SWRabsTotalGround-SWRabsDir_nT.SWRabsTotalGround
        SWRabsDiff_nT.SWRabsTotalCanyon = SWRabs_nT.SWRabsTotalCanyon-SWRabsDir_nT.SWRabsTotalCanyon

        # Energy Balance of shortwave radiation
        class SWREB_nT_Def():
            pass
        SWREB_nT = SWREB_nT_Def()
        SWREB_nT.SWREBGroundImp = SWRin_nT.SWRinGroundImp - SWRout_nT.SWRoutGroundImp - SWRabs_nT.SWRabsGroundImp
        SWREB_nT.SWREBGroundBare = SWRin_nT.SWRinGroundBare - SWRout_nT.SWRoutGroundBare - SWRabs_nT.SWRabsGroundBare
        SWREB_nT.SWREBGroundVeg = SWRin_nT.SWRinGroundVeg - SWRout_nT.SWRoutGroundVeg - SWRabs_nT.SWRabsGroundVeg
        SWREB_nT.SWREBTree = 0
        SWREB_nT.SWREBWallSun = SWRin_nT.SWRinWallSun-SWRout_nT.SWRoutWallSun - SWRabs_nT.SWRabsWallSun
        SWREB_nT.SWREBWallShade = SWRin_nT.SWRinWallShade-SWRout_nT.SWRoutWallShade - SWRabs_nT.SWRabsWallShade
        SWREB_nT.SWREBTotalGround = SWRin_nT.SWRinTotalGround-SWRout_nT.SWRoutTotalGround - SWRabs_nT.SWRabsTotalGround
        SWREB_nT.SWREBTotalCanyon = SWRin_nT.SWRinTotalCanyon-SWRout_nT.SWRoutTotalCanyon - SWRabs_nT.SWRabsTotalCanyon

        if abs(SWREB_nT.SWREBGroundImp) >= 10**(-6):
            print('SWREB_nT.SWREBGroundImp is not 0')
        elif abs(SWREB_nT.SWREBGroundBare) >= 10**(-6):
            print('SWREB_nT.SWREBGroundBare is not 0')
        elif abs(SWREB_nT.SWREBGroundVeg) >= 10**(-6):
            print('SWREB_nT.SWREBGroundVeg	 is not 0')
        elif abs(SWREB_nT.SWREBWallSun) >= 10**(-6):
            print('SWREB_nT.SWREBWallSun is not 0')
        elif abs(SWREB_nT.SWREBWallShade) >= 10**(-6):
            print('SWREB_nT.SWREBWallShade is not 0')
        elif abs(SWREB_nT.SWREBTotalGround) >= 10**(-6):
            print('SWREB_nT.SWREBTotalGround is not ')
        elif abs(SWREB_nT.SWREBTotalCanyon) >= 10**(-6):
            print('SWREB_nT.SWREBTotalCanyon is not 0')

        return SWRin_nT,SWRout_nT,SWRabs_nT,SWRabsDir_nT,SWRabsDiff_nT,SWREB_nT

    def SWRabsorbedWithTrees(self,h_can,w_can,h_tree,r_tree,d_tree,fgveg,fgbare,fgimp,aw,agveg,agbare,agimp,at,LAIt,
                             SWR_dir,SWR_diff,theta_Z,theta_n,ViewFactor,ParVegTree):

        F_gs_T = ViewFactor.F_gs_T
        F_gt_T = ViewFactor.F_gt_T
        F_gw_T = ViewFactor.F_gw_T
        F_ww_T = ViewFactor.F_ww_T
        F_wt_T = ViewFactor.F_wt_T
        F_wg_T = ViewFactor.F_wg_T
        F_ws_T = ViewFactor.F_ws_T
        F_sg_T = ViewFactor.F_sg_T
        F_sw_T = ViewFactor.F_sw_T
        F_st_T = ViewFactor.F_st_T
        F_tg_T = ViewFactor.F_tg_T
        F_tw_T = ViewFactor.F_tw_T
        F_ts_T = ViewFactor.F_ts_T
        F_tt_T = ViewFactor.F_tt_T

        # normalized surface areas
        A_s = w_can
        A_g = w_can
        A_w = h_can
        # There are 2 trees. Hence, the area of tree is twice a circle
        A_t = 2 * 2 * numpy.pi * r_tree

        # load shortwave radiation
        SWRdir_ground, SWRdir_wallsun, SWRdir_wallshade, SWRdir_tree = self.DirectSWRSurfaces(h_can, w_can, d_tree, h_tree,
                                                                                              r_tree, theta_Z, theta_n,
                                                                                              SWR_dir, LAIt, 1,ParVegTree)

        # Balance direct shortwave radiation in
        EB_SWRdir = SWR_dir-(SWRdir_ground*A_g/A_g+SWRdir_wallsun*A_w/A_g+SWRdir_wallshade*A_w/A_g+SWRdir_tree*A_t/A_g)
        EB_SWRdiff = SWR_diff-(F_sg_T*SWR_diff+F_sw_T*SWR_diff+F_sw_T*SWR_diff+F_st_T*SWR_diff)

        if abs(EB_SWRdir) >= 10 ** (-1):
            print('EB_SWRdir is not 0',EB_SWRdir)
        if abs(EB_SWRdiff) >= 10 ** (-1):
            print('EB_SWRdiff is not 0',EB_SWRdiff)

        SVF = numpy.zeros(4)
        SVF[0] = F_gs_T+F_gt_T+2*F_gw_T
        SVF[1] = F_ww_T+F_wt_T+F_wg_T+F_ws_T
        SVF[2] = F_sg_T+2*F_sw_T+F_st_T
        SVF[3] = F_ts_T+2*F_tw_T+F_tt_T+F_tg_T

        SVF2 = numpy.zeros(4)
        SVF2[0] = F_gs_T+2*F_ws_T*A_w+F_ts_T*A_t
        SVF2[1] = F_sg_T+2*F_wg_T*A_w+F_tg_T*A_t
        SVF2[2] = F_ww_T+F_sw_T*A_g/A_w+F_gw_T*A_g/A_w+F_tw_T*A_t/A_w
        SVF2[3] = F_gt_T*A_g/A_t+2*F_wt_T*A_w/A_t+F_tt_T

        for i in range(0, len(SVF)):
            if SVF[i] < 0.999 or SVF[i] > 1.001:
                print('The view factors do not add up to 1 for a canyon with trees')

        if fgimp > 0:
            Cimp = 1
        else:
            Cimp = 0
        if fgbare > 0:
            Cbare = 1
        else:
            Cbare = 0
        if fgveg > 0:
            Cveg = 1
        else:
            Cveg = 0

        ai = [agveg, agbare, agimp, aw, aw, at, 0]

        # View factor matrix to solve for infinite reflections equation
        Tij = numpy.array([[1,0,0, -agveg*F_gw_T*Cveg, -agveg*F_gw_T*Cveg, -agveg*F_gt_T*Cveg, -agveg*F_gs_T*Cveg],
               [0,1,0, -agbare*F_gw_T*Cbare, -agbare*F_gw_T*Cbare, -agbare*F_gt_T*Cbare, -agbare*F_gs_T*Cbare],
               [0,0,1, -agimp*F_gw_T*Cimp, -agimp*F_gw_T*Cimp, -agimp*F_gt_T*Cimp, -agimp*F_gs_T*Cimp],
               [-aw*F_wg_T*fgveg*Cveg,-aw*F_wg_T*fgbare*Cbare,-aw*F_wg_T*fgimp*Cimp, 1, -aw*F_ww_T, -aw*F_wt_T, -aw*F_ws_T],
               [-aw*F_wg_T*fgveg*Cveg,-aw*F_wg_T*fgbare*Cbare,-aw*F_wg_T*fgimp*Cimp, -aw*F_ww_T, 1, -aw*F_wt_T, -aw*F_ws_T],
               [-at*F_tg_T*fgveg*Cveg,-at*F_tg_T*fgbare*Cbare,-at*F_tg_T*fgimp*Cimp, -at*F_tw_T, -at*F_tw_T, 1-at*F_tt_T, -at*F_ts_T],
               [0, 0, 0, 0, 0, 0, 1]])

        # Incoming shortwave radiation from sky
        Omega_i = numpy.array([agveg * SWRdir_ground * Cveg,
                   agbare * SWRdir_ground * Cbare,
                   agimp * SWRdir_ground * Cimp,
                   aw * SWRdir_wallsun,
                   aw * 0,
                   at * SWRdir_tree,
                   SWR_diff])

        # Outgoing radiation per surface
        B_i = numpy.linalg.solve(Tij,Omega_i)

        if B_i[6] != SWR_diff:
            print('Incoming lonwave radiation and emitted longwave radiation from the sky after the matrix inversion are not equal')

        # Incoming shortwave radiation at each surface A_i
        Tij2 = numpy.array([[0, 0, 0, F_gw_T*Cveg, F_gw_T*Cveg, F_gt_T*Cveg, F_gs_T*Cveg],
                [0, 0, 0, F_gw_T*Cbare, F_gw_T*Cbare, F_gt_T*Cbare, F_gs_T*Cbare],
                [0, 0, 0, F_gw_T*Cimp, F_gw_T*Cimp, F_gt_T*Cimp, F_gs_T*Cimp],
                [F_wg_T*fgveg*Cveg, F_wg_T*fgbare*Cbare, F_wg_T*fgimp*Cimp, 0, F_ww_T, F_wt_T, F_ws_T],
                [F_wg_T*fgveg*Cveg, F_wg_T*fgbare*Cbare, F_wg_T*fgimp*Cimp, F_ww_T, 0, F_wt_T, F_ws_T],
                [F_tg_T*fgveg*Cveg, F_tg_T*fgbare*Cbare, F_tg_T*fgimp*Cimp, F_tw_T, F_tw_T, F_tt_T, F_ts_T],
                [0, 0, 0, 0, 0, 0, 0]])

        SWRdir_i = numpy.array([SWRdir_ground * Cveg,
                    SWRdir_ground * Cbare,
                    SWRdir_ground * Cimp,
                    SWRdir_wallsun,
                    0,
                    SWRdir_tree,
                    0])

        A_i1 = numpy.dot(Tij2, B_i) + SWRdir_i  # Incoming radiation [W/m^2] per m^2 surface area
        A_i = B_i / ai  # Incoming radiation [W/m^2] per m^2 surface area
        for i in range(0, len(ai)):
            if ai[i] == 0:
                A_i[i] = A_i1[i]
        A_i[6] = 0  # Assumption: The sky has a fixed emission of LWR. Hence, Qnet is 0.
        # Absorbed shortwave radiation at ech surface Qnet_i
        Qnet_i = A_i - B_i

        # Assignment
        SWRout_i = B_i  # Outgoing radiation [W/m^2] per m^2 surface area
        SWRin_i = A_i  # Incoming radiation [W/m^2] per m^2 surface area
        SWRnet_i = Qnet_i  # Net absorbed radiation [W/m^2] per m^2 surface area

        # Energy balance
        SWRin_atm = SWR_dir + SWR_diff

        TotalSWRSurface_in = SWRin_i[0] * fgveg * A_g / A_g + SWRin_i[1] * fgbare * A_g / A_g + SWRin_i[2] * fgimp * A_g / A_g +\
                             SWRin_i[3] * A_w / A_g + SWRin_i[4] * A_w / A_g + SWRin_i[5]*A_t/A_g

        TotalSWRSurface_abs = SWRnet_i[0] * fgveg * A_g / A_g + SWRnet_i[1] * fgbare * A_g / A_g + SWRnet_i[2] * fgimp * A_g / A_g +\
                              SWRnet_i[3] * A_w / A_g + SWRnet_i[4] * A_w / A_g + SWRnet_i[5]*A_t/A_g

        TotalSWRSurface_out = SWRout_i[0] * fgveg * A_g / A_s + SWRout_i[1] * fgbare * A_g / A_s + SWRout_i[2] * fgimp * A_g / A_s + \
                              SWRout_i[3] * A_w / A_s + SWRout_i[4] * A_w / A_s + SWRout_i[5]*A_t/A_s

        TotalSWRref_to_atm = SWRout_i[0] * F_sg_T * fgveg + SWRout_i[1] * F_sg_T * fgbare + SWRout_i[2] * F_sg_T * fgimp +\
                             SWRout_i[3] * F_sw_T + SWRout_i[4] * F_sw_T + SWRout_i[5]*F_st_T

        TotalSWRref_to_atm2 = SWRout_i[0] * F_gs_T * fgveg + SWRout_i[1] * F_gs_T * fgbare + SWRout_i[2] * F_gs_T * fgimp +\
                             SWRout_i[3] * F_ws_T + SWRout_i[4] * F_ws_T + SWRout_i[5]*F_ts_T


        EBSurface = TotalSWRSurface_in - TotalSWRSurface_abs - TotalSWRSurface_out
        EBCanyon = SWRin_atm - TotalSWRSurface_abs - TotalSWRref_to_atm

        if abs(EBSurface)>=10**(-1):
            print('EBSurface is not 0',EBSurface)
        if abs(EBCanyon)>=10**(-1):
            print('EBCanyon is not 0',EBCanyon)


        # Incoming shortwave radiation
        class SWRin_T_Def():
            pass
        SWRin_T = SWRin_T_Def()
        SWRin_T.SWRinGroundImp = SWRin_i[2] * Cimp
        SWRin_T.SWRinGroundBare = SWRin_i[1] * Cbare
        SWRin_T.SWRinGroundVeg = SWRin_i[0] * Cveg
        SWRin_T.SWRinTree = SWRin_i[5]
        SWRin_T.SWRinWallSun = SWRin_i[3]
        SWRin_T.SWRinWallShade = SWRin_i[4]
        SWRin_T.SWRinTotalGround = fgveg * SWRin_i[0] + fgbare * SWRin_i[1] + fgimp * SWRin_i[2]
        SWRin_T.SWRinTotalCanyon = SWRin_i[0] * fgveg * A_g / A_g + SWRin_i[1] * fgbare * A_g / A_g + SWRin_i[2] * fgimp * A_g / A_g + \
                                    SWRin_i[3] * A_w / A_g + SWRin_i[4] * A_w / A_g  +SWRin_i[5]*A_t/A_g

        # Outgoing shortwave radiation
        class SWRout_T_Def():
            pass
        SWRout_T = SWRout_T_Def()
        SWRout_T.SWRoutGroundImp = SWRout_i[2] * Cimp
        SWRout_T.SWRoutGroundBare = SWRout_i[1] * Cbare
        SWRout_T.SWRoutGroundVeg = SWRout_i[0] * Cveg
        SWRout_T.SWRoutTree = SWRout_i[5]
        SWRout_T.SWRoutWallSun = SWRout_i[3]
        SWRout_T.SWRoutWallShade = SWRout_i[4]
        SWRout_T.SWRoutTotalGround = fgveg * SWRout_i[0] + fgbare * SWRout_i[1] + fgimp * SWRout_i[2]
        SWRout_T.SWRoutTotalCanyon = SWRout_i[0] * fgveg * A_g / A_g + SWRout_i[1] * fgbare * A_g / A_g + SWRout_i[2] * fgimp * A_g / A_g + \
                                    SWRout_i[3] * A_w / A_g + SWRout_i[4] * A_w / A_g + SWRout_i[5]*A_t/A_g

        # Absorbed shortwave radiation
        class SWRabs_T_Def():
            pass
        SWRabs_T = SWRabs_T_Def()
        SWRabs_T.SWRabsGroundImp = SWRnet_i[2] * Cimp
        SWRabs_T.SWRabsGroundBare = SWRnet_i[1] * Cbare
        SWRabs_T.SWRabsGroundVeg = SWRnet_i[0] * Cveg
        SWRabs_T.SWRabsTree = SWRnet_i[5]
        SWRabs_T.SWRabsWallSun = SWRnet_i[3]
        SWRabs_T.SWRabsWallShade = SWRnet_i[4]
        SWRabs_T.SWRabsTotalGround = fgveg * SWRnet_i[0] + fgbare * SWRnet_i[1] + fgimp * SWRnet_i[2]
        SWRabs_T.SWRabsTotalCanyon = SWRnet_i[0] * fgveg * A_g / A_g + SWRnet_i[1] * fgbare * A_g / A_g + SWRnet_i[2] * fgimp * A_g / A_g + \
                                    SWRnet_i[3] * A_w / A_g + SWRnet_i[4] * A_w / A_g + SWRnet_i[5]*A_t/A_g

        # Direct absorbed shortwave radiation
        class SWRabsDir_T_Def():
            pass
        SWRabsDir_T = SWRabsDir_T_Def()
        SWRabsDir_T.SWRabsGroundImp = (1-agimp)*SWRdir_ground*Cimp
        SWRabsDir_T.SWRabsGroundBare = (1-agbare)*SWRdir_ground*Cbare
        SWRabsDir_T.SWRabsGroundVeg = (1-agveg)*SWRdir_ground*Cveg
        SWRabsDir_T.SWRabsTree = (1-at)*SWRdir_tree
        SWRabsDir_T.SWRabsWallSun = (1-aw)*SWRdir_wallsun
        SWRabsDir_T.SWRabsWallShade = (1-aw)*SWRdir_wallshade
        SWRabsDir_T.SWRabsTotalGround = fgveg*(1-agveg)*SWRdir_ground+fgbare*(1-agbare)*SWRdir_ground+fgimp*(1-agimp)*SWRdir_ground
        SWRabsDir_T.SWRabsTotalCanyon = fgveg*(1-agveg)*SWRdir_ground*A_g/A_g+fgbare*(1-agbare)*SWRdir_ground*A_g/A_g+\
                                        fgimp*(1-agimp)*SWRdir_ground*A_g/A_g + (1-aw)*SWRdir_wallsun*A_w/A_g + (1-aw)*SWRdir_wallshade*A_w/A_g + \
                                        (1-at)*SWRdir_tree*A_t/A_g


        # Diffuse absorbed shortwave radiation
        class SWRabsDiff_T_Def():
            pass
        SWRabsDiff_T = SWRabsDiff_T_Def()
        SWRabsDiff_T.SWRabsGroundImp = (SWRabs_T.SWRabsGroundImp-SWRabsDir_T.SWRabsGroundImp)*Cimp
        SWRabsDiff_T.SWRabsGroundBare = (SWRabs_T.SWRabsGroundBare-SWRabsDir_T.SWRabsGroundBare)*Cbare
        SWRabsDiff_T.SWRabsGroundVeg = (SWRabs_T.SWRabsGroundVeg-SWRabsDir_T.SWRabsGroundVeg)*Cveg
        SWRabsDiff_T.SWRabsTree = SWRabs_T.SWRabsTree-SWRabsDir_T.SWRabsTree
        SWRabsDiff_T.SWRabsWallSun = SWRabs_T.SWRabsWallSun-SWRabsDir_T.SWRabsWallSun
        SWRabsDiff_T.SWRabsWallShade = SWRabs_T.SWRabsWallShade-SWRabsDir_T.SWRabsWallShade
        SWRabsDiff_T.SWRabsTotalGround = SWRabs_T.SWRabsTotalGround-SWRabsDir_T.SWRabsTotalGround
        SWRabsDiff_T.SWRabsTotalCanyon = SWRabs_T.SWRabsTotalCanyon-SWRabsDir_T.SWRabsTotalCanyon

        # Energy Balance of shortwave radiation
        class SWREB_T_Def():
            pass
        SWREB_T = SWREB_T_Def()
        SWREB_T.SWREBGroundImp = SWRin_T.SWRinGroundImp - SWRout_T.SWRoutGroundImp - SWRabs_T.SWRabsGroundImp
        SWREB_T.SWREBGroundBare = SWRin_T.SWRinGroundBare - SWRout_T.SWRoutGroundBare - SWRabs_T.SWRabsGroundBare
        SWREB_T.SWREBGroundVeg = SWRin_T.SWRinGroundVeg - SWRout_T.SWRoutGroundVeg - SWRabs_T.SWRabsGroundVeg
        SWREB_T.SWREBTree = SWRin_T.SWRinTree - SWRout_T.SWRoutTree - SWRabs_T.SWRabsTree
        SWREB_T.SWREBWallSun = SWRin_T.SWRinWallSun-SWRout_T.SWRoutWallSun - SWRabs_T.SWRabsWallSun
        SWREB_T.SWREBWallShade = SWRin_T.SWRinWallShade-SWRout_T.SWRoutWallShade - SWRabs_T.SWRabsWallShade
        SWREB_T.SWREBTotalGround = SWRin_T.SWRinTotalGround-SWRout_T.SWRoutTotalGround - SWRabs_T.SWRabsTotalGround
        SWREB_T.SWREBTotalCanyon = SWRin_T.SWRinTotalCanyon-SWRout_T.SWRoutTotalCanyon - SWRabs_T.SWRabsTotalCanyon

        if abs(SWREB_T.SWREBGroundImp) >= 10**(-6):
            print('SWREB_T.SWREBGroundImp is not 0')
        elif abs(SWREB_T.SWREBGroundBare) >= 10**(-6):
            print('SWREB_T.SWREBGroundBare is not 0')
        elif abs(SWREB_T.SWREBGroundVeg) >= 10**(-6):
            print('SWREB_T.SWREBGroundVeg	 is not 0')
        elif abs(SWREB_T.SWREBWallSun) >= 10**(-6):
            print('SWREB_T.SWREBWallSun is not 0')
        elif abs(SWREB_T.SWREBWallShade) >= 10**(-6):
            print('SWREB_T.SWREBWallShade is not 0')
        elif abs(SWREB_T.SWREBTotalGround) >= 10**(-6):
            print('SWREB_T.SWREBTotalGround is not ')
        elif abs(SWREB_T.SWREBTotalCanyon) >= 10**(-6):
            print('SWREB_T.SWREBTotalCanyon is not 0')
        elif abs(SWREB_T.SWREBTree) >= 10**(-6):
            print('SWREB_T.SWREBTree is not 0')

        return SWRin_T,SWRout_T,SWRabs_T,SWRabsDir_T,SWRabsDiff_T,SWREB_T


    def DirectSWRSurfaces(self,h_can,w_can,d_tree,h_tree,r_tree,theta_Z,theta_n,SWR_dir,LAIt,trees,ParVegTree):

        Kopt_T = ParVegTree.Kopt

        # Calculation of SWRdir_t
        if trees == 0:
            tau = 0
            SWRdir_t = 0
        else:
            SWR_tree1, SWR_tree2 = self.DirectSWRTrees(h_can,d_tree,h_tree,r_tree,theta_Z,theta_n,SWR_dir)
            # Calculate how much shortwave radiation passes through the trees
            tau = numpy.exp(-Kopt_T * LAIt)
            # averaging over the two trees
            SWRdir_t = (1 - tau) * (SWR_tree1 + SWR_tree2) / 2

        # Calculation of SWRdir_g, SWRdir_wsun, SWRdir_wshd
        if trees == 0:
            X_Shadow, X_Tree, n_Shadow, n_Tree = self.ShadowLengthNoTree(h_can, w_can, theta_Z, theta_n)
        else:
            X_Shadow, X_Tree, n_Shadow, n_Tree = self.ShadowLengthWithTrees(h_can, w_can, d_tree, h_tree, r_tree, theta_Z, theta_n)

        Xsi = math.tan(theta_Z) * abs(math.sin(theta_n))
        #print(math.tan(theta_Z))
        #print(math.sin(theta_n))
        SWRdir_g = SWR_dir * (1 - X_Shadow + tau * X_Tree)
        SWRdir_wsun = SWR_dir * Xsi * (1 - n_Shadow + tau * n_Tree)
        SWRdir_wshd = 0

        A_g = w_can
        A_w = h_can
        A_t = 2 * numpy.pi * r_tree

        check_total = A_g / A_g * SWRdir_g + A_w / A_g * SWRdir_wsun + A_w / A_g * SWRdir_wshd + 2 * A_t / A_g * SWRdir_t

        if abs(check_total - SWR_dir) > 1e-10:
            delta = check_total - SWR_dir
            if delta < 0:
                # the energy excess or deficit is distributed to the trees
                SWRdir_t = SWRdir_t - delta * A_g / (2 * A_t)
            else:
                # the energy excess or deficit is distributed to the trees
                SWRdir_t = SWRdir_t - delta * A_g / (2 * A_t)

        return SWRdir_g,SWRdir_wsun,SWRdir_wshd,SWRdir_t


    def DirectSWRTrees(self,h_can,d_tree,h_tree,r_tree,theta_Z,theta_n,SWR_dir):

        Xsi = math.tan(theta_Z) * abs(math.sin(theta_n))

        tan_theta1 = ((1-d_tree)*(h_can-h_tree)+r_tree*numpy.sqrt((1-d_tree)**2+(h_can-h_tree)**2-r_tree**2))/((h_can-h_tree)**2-r_tree**2)
        tan_theta2 = ((1-d_tree)*(h_can-h_tree)-r_tree*numpy.sqrt((1-d_tree)**2+(h_can-h_tree)**2-r_tree**2))/((h_can-h_tree)**2-r_tree**2)
        tan_theta3 = (d_tree*(h_can-h_tree)+r_tree*numpy.sqrt(d_tree**2+(h_can-h_tree)**2-r_tree**2))/((h_can-h_tree)**2-r_tree**2)
        tan_theta4 = (d_tree*(h_can-h_tree)-r_tree*numpy.sqrt(d_tree**2+(h_can-h_tree)**2-r_tree**2))/((h_can-h_tree)**2-r_tree**2)

        if Xsi >= tan_theta1:
            # Tree 1 is completely shaded
            SWR_tree1 = 0
        elif Xsi < tan_theta1 and Xsi >= tan_theta2:
            # Tree 1 is partially sunlit
            SWR_tree1 = SWR_dir*(r_tree*numpy.sqrt(1+Xsi**2) + (1 - d_tree) - (h_can - h_tree) * Xsi) / (2 * numpy.pi * r_tree)
        elif Xsi < tan_theta2:
            # tree 1 is completely sunlit
            SWR_tree1 = SWR_dir * (2 * r_tree * numpy.sqrt(1 + Xsi**2)) / (2 * numpy.pi * r_tree)
        else:
            # Account for weird angles at night (angles = NaN)
            SWR_tree1 = 0

        if Xsi >= tan_theta3:
            # Tree 2 is completely shaded
            SWR_tree2 = 0
        elif Xsi < tan_theta3 and Xsi >= tan_theta4:
            # Tree 2 is partially sunlit
            SWR_tree2 = SWR_dir * (r_tree * numpy.sqrt(1 + Xsi**2) + d_tree - (h_can - h_tree) * Xsi) / (2 * numpy.pi * r_tree)
        elif Xsi<tan_theta4:
            # tree 1 is completely sunlit
            SWR_tree2 = SWR_dir * (2 * r_tree * numpy.sqrt(1 + Xsi**2)) / (2 * numpy.pi * r_tree)
        else:
            # Account for weird angles at night (angles = NaN)
            SWR_tree2 = 0

        return SWR_tree1,SWR_tree2

    def ShadowLengthNoTree(self,h_can,w_can,theta_Z,theta_n):

        Xsi = math.tan(theta_Z) * abs(math.sin(theta_n))

        # Shadow by the Wall
        # shadow cast on the ground by the wall
        X_Shadow = h_can * Xsi
        # shadow cast on the opposite wall by the wall
        n_Shadow = h_can - w_can / Xsi

        if abs(X_Shadow) < w_can:
            n_Shadow = 0
        else:
            X_Shadow = w_can

        if n_Shadow < h_can:
            n_Shadow = n_Shadow
        else:
            n_Shadow = h_can

        # NOTE : the origin (0,0) is the lower left corner of the canyon
        X_Shadow = X_Shadow
        X_tree = 0
        n_Shadow = n_Shadow / h_can
        n_tree = 0

        return X_Shadow,X_tree,n_Shadow,n_tree


    def ShadowLengthWithTrees(self,h_can,w_can,d_tree,h_tree,r_tree,theta_Z,theta_n):

        Xsi = math.tan(theta_Z) * abs(math.sin(theta_n))

        # Shadow by the Wall
        X_wall = h_can * Xsi
        n_wall = h_can - w_can / Xsi

        if abs(X_wall) < w_can:
            n_wall = 0
        else:
            X_wall = w_can

        x0 = max(0., w_can - h_can * Xsi)
        y0 = max(0., h_can - w_can / Xsi)
        secXsi = numpy.sqrt(1 + Xsi ** 2)
        cosecXsi = numpy.sqrt(1 + 1 / Xsi ** 2)

        # Shadow by the Tree 1
        x1 = max(0, d_tree - h_tree * Xsi - r_tree * secXsi)
        y1 = max(0, h_tree - (w_can - d_tree) / Xsi - r_tree * cosecXsi)
        x2 = max(0, d_tree - h_tree * Xsi + r_tree * secXsi)
        y2 = max(0, h_tree - (w_can - d_tree) / Xsi + r_tree * cosecXsi)

        X_Tree1 = x2 - x1

        # Shadow by the Tree 2
        x3 = max(0, w_can - d_tree - h_tree * Xsi - r_tree * secXsi)
        y3 = max(0, h_tree - d_tree / Xsi - r_tree * cosecXsi)
        x4 = max(0, w_can - d_tree - h_tree * Xsi + r_tree * secXsi)
        y4 = max(0, h_tree - d_tree / Xsi + r_tree * cosecXsi)

        X_Tree2 = x4 - x3
        n_Tree1 = y4 - y3
        n_Tree2 = y2 - y1

        # Total shadow length on the ground by wall, Tree 1, and Tree 2
        delta = max(0, x2 - x0)

        if x0 < x4:
            X_Shadow = w_can - min(x0,x3) + X_Tree1 - delta
            if x0 < x3:
                X_Tree = X_Tree1 - delta
            else:
                X_Tree = X_Tree1 + x0 - x3
        elif x0 >= x4:
            X_Shadow = X_wall + X_Tree1 + X_Tree2
            X_Tree = X_Tree1 + X_Tree2

        # Total shadow length on the wall by wall, Tree 1, and Tree 2
        lowest_shaded = max(y0,y2)
        if y3 > lowest_shaded:
            n_Shadow = n_Tree1 + lowest_shaded
            if y2 > n_wall:
                n_Tree = n_Tree1 + y2 - n_wall
            elif y2 <= n_wall:
                n_Tree = n_Tree1
            elif y1 > n_wall:
                n_Tree = n_Tree1 + n_Tree2
        else:
            n_Shadow = max([y0,y1,y2,y3,y4])
            if (y4 > n_wall):
                n_Tree = y4 - n_wall
            else:
                n_Tree = 0.

        n_Tree = n_Tree / h_can
        n_Shadow = n_Shadow / h_can

        return X_Shadow,X_Tree,n_Shadow,n_Tree


    def VFUrbanCanyon(self,ViewFactorCal_Param,Gemeotry_m,geometry,Person,ParTree,ViewFactor_file_text):


        if ViewFactorCal_Param.OPTION_RAY == 1:

            ViewFactor_file = numpy.loadtxt(ViewFactor_file_text)

            F_gs_nT = ViewFactor_file[0]
            F_gw_nT = ViewFactor_file[1]
            F_ww_nT = ViewFactor_file[2]
            F_wg_nT = ViewFactor_file[3]
            F_ws_nT = ViewFactor_file[4]
            F_sg_nT = ViewFactor_file[5]
            F_sw_nT = ViewFactor_file[6]
            F_gs_T = ViewFactor_file[7]
            F_gt_T = ViewFactor_file[8]
            F_gw_T = ViewFactor_file[9]
            F_ww_T = ViewFactor_file[10]
            F_wt_T = ViewFactor_file[11]
            F_wg_T = ViewFactor_file[12]
            F_ws_T = ViewFactor_file[13]
            F_sg_T = ViewFactor_file[14]
            F_sw_T = ViewFactor_file[15]
            F_st_T = ViewFactor_file[16]
            F_tg_T = ViewFactor_file[17]
            F_tw_T = ViewFactor_file[18]
            F_ts_T = ViewFactor_file[19]
            F_tt_T = ViewFactor_file[20]

            Sum_g = F_gs_T + F_gt_T + 2 * F_gw_T
            Sum_w = F_ww_T + F_wt_T + F_wg_T + F_ws_T
            Sum_s = F_sg_T + 2 * F_sw_T + F_st_T
            Sum_t = F_ts_T + 2 * F_tw_T + F_tt_T + F_tg_T

            F_pg = ViewFactor_file[21]
            F_ps = ViewFactor_file[22]
            F_pt = ViewFactor_file[23]
            F_pw = ViewFactor_file[24]

        else:

            if ParTree.trees == 0:
                geometry.radius_tree = 0
                geometry.htree = -Gemeotry_m.Height_canyon / 10000
                geometry.distance_tree = 0

            # compute view factors with monte carlo ray tracing
            F_gs_T, F_gt_T, F_gw_T, F_ww_T, F_wt_T, F_wg_T, F_ws_T, F_ts_T, F_tw_T, F_tt_T, F_tg_T, F_sg_T, F_sw_T, F_st_T,\
            F_pg, F_ps, F_pw, F_pt, VFRayTracingRaw_T, VFRayTracing_T = \
                self.VFRayTracingReciprocity(Gemeotry_m.Height_canyon, Gemeotry_m.Width_canyon, geometry.radius_tree,geometry.htree,
                                             geometry.distance_tree, Person, int(ViewFactorCal_Param.MCSampleSize), int(ViewFactorCal_Param.NRays))

        # calculate view factors with analytical solutions
        F_gs_nT, F_gt_nT, F_gw_nT, F_ww_nT, F_wt_nT, F_wg_nT, F_ws_nT, F_ts_nT, F_tw_nT, F_tt_nT, F_tg_nT, F_sg_nT, \
        F_sw_nT,F_st_nT, ViewFactor_nT = self.VFAnalytical(Gemeotry_m.Height_canyon, Gemeotry_m.Width_canyon)


        VF_values = [F_gs_nT,F_gw_nT,F_ww_nT,F_wg_nT,F_ws_nT,F_sg_nT,F_sw_nT,F_gs_T,F_gt_T,F_gw_T,F_ww_T,F_wt_T,F_wg_T,F_ws_T,
                     F_sg_T,F_sw_T,F_st_T,F_tg_T,F_tw_T,F_ts_T,F_tt_T, F_pg, F_ps, F_pt, F_pw]
        outputFile_VF = open(ViewFactor_file_text, "w")
        outputFile_VF.write("#### \t Vertical City Weather Generator (VCWG)  \t #### \n")
        outputFile_VF.write("# View Factors \n")
        outputFile_VF.write("# F_gs_nT	F_gw_nT	F_ww_nT	F_wg_nT	F_ws_nT	F_sg_nT	F_sw_nT	F_gs_T	F_gt_T	F_gw_T	F_ww_T	F_wt_T	F_wg_T	F_ws_T	F_sg_T	F_sw_T	F_st_T	F_tg_T	F_tw_T	F_ts_T	F_tt_T F_pg, F_ps, F_pt, F_pw \n")
        for i in range(25):
            outputFile_VF.write("%f " % (VF_values[i]))
        outputFile_VF.close()


        class ViewFactor_Def():
            pass
        ViewFactor = ViewFactor_Def()
        ViewFactor.F_gs_nT = F_gs_nT
        ViewFactor.F_gw_nT = F_gw_nT
        ViewFactor.F_ww_nT = F_ww_nT
        ViewFactor.F_wg_nT = F_wg_nT
        ViewFactor.F_ws_nT = F_ws_nT
        ViewFactor.F_sg_nT = F_sg_nT
        ViewFactor.F_sw_nT = F_sw_nT
        ViewFactor.F_gs_T = F_gs_T
        ViewFactor.F_gt_T = F_gt_T
        ViewFactor.F_gw_T = F_gw_T
        ViewFactor.F_ww_T = F_ww_T
        ViewFactor.F_wt_T = F_wt_T
        ViewFactor.F_wg_T = F_wg_T
        ViewFactor.F_ws_T = F_ws_T
        ViewFactor.F_sg_T = F_sg_T
        ViewFactor.F_sw_T = F_sw_T
        ViewFactor.F_st_T = F_st_T
        ViewFactor.F_tg_T = F_tg_T
        ViewFactor.F_tw_T = F_tw_T
        ViewFactor.F_ts_T = F_ts_T
        ViewFactor.F_tt_T = F_tt_T

        class ViewFactorPoint_Def():
            pass
        ViewFactorPoint = ViewFactorPoint_Def()
        ViewFactorPoint.F_pg = F_pg
        ViewFactorPoint.F_ps = F_ps
        ViewFactorPoint.F_pt = F_pt
        ViewFactorPoint.F_pw = F_pw

        return  ViewFactor,ViewFactorPoint

    def VFRayTracingReciprocity(self,H, W, a, ht, d, Person, MCSampleSize, NRays):

        _F_gs_T_, _F_gt_T_, _F_gw_T_, _F_ww_T_, _F_wt_T_, _F_wg_T_, _F_ws_T_, _F_ts_T_, _F_tw_T_, _F_tt_T_, _F_tg_T_, \
        _F_sg_T_, _F_sw_T_, _F_st_T_, F_pg, F_ps, F_pw, F_pt, VFRayTracingRaw_T = \
            self.VFRayTracing(H, W, a, ht, d, Person, MCSampleSize, NRays)

        h = H / W
        w = W / W
        ratio = h / w

        Sum = numpy.zeros(4)
        Sum2 = numpy.zeros(4)

        if a == 0:

            # The view factor taken from the ray tracing is F_gs_T
            F_gs_T = VFRayTracingRaw_T.F_gs_T
            # factor 0.5 because there are 2 walls that are seen by the ground
            F_gw_T = 0.5 * (1 - F_gs_T)
            F_gt_T = 0

            F_sg_T = F_gs_T * w / w
            F_sw_T = F_gw_T * w / w
            F_st_T = 0

            F_wg_T = F_gw_T * w / h
            F_ws_T = F_sw_T * w / h
            F_ww_T = 1 - F_wg_T - F_ws_T
            F_wt_T = 0

            F_tg_T = 0
            F_ts_T = 0
            F_tw_T = 0
            F_tt_T = 0

            Sum[0] = F_gs_T + 2 * F_gw_T
            Sum[1] = F_ww_T + F_wg_T + F_ws_T
            Sum[2] = F_sg_T + 2 * F_sw_T
            Sum[3] = 0

            Sum2[0] = F_sg_T * w / w + 2 * F_wg_T * h / w
            Sum2[1] = F_ww_T * h / h + F_gw_T * w / h + F_sw_T * w / h
            Sum2[2] = F_gs_T * w / w + 2 * F_ws_T * h / w
            Sum2[3] = 0

        else:
            # The view factors taken from the ray tracing are F_st_T, F_gs_T,F_gt_T, F_wt_T
            Atree = 2 * 2 * numpy.pi * a

            F_gs_T = VFRayTracingRaw_T.F_gs_T
            F_gt_T = VFRayTracingRaw_T.F_gt_T
            # factor 0.5 because there are 2 walls that are seen by the ground
            F_gw_T = 0.5 * (1 - F_gs_T - F_gt_T)

            F_sg_T = F_gs_T * w / w
            F_st_T = VFRayTracingRaw_T.F_st_T
            # factor 0.5 because there are 2 walls that are seen by the ground
            F_sw_T = 0.5 * (1 - F_sg_T - F_st_T)

            F_wg_T = F_gw_T * w / h
            F_ws_T = F_sw_T * w / h
            F_wt_T = VFRayTracingRaw_T.F_wt_T
            F_ww_T = 1 - F_wg_T - F_ws_T - F_wt_T

            F_ts_T = F_st_T * w / Atree
            F_tw_T = F_wt_T * h / Atree
            F_tg_T = F_gt_T * w / Atree
            F_tt_T = 1 - F_ts_T - 2 * F_tw_T - F_tg_T

            Sum[0] = F_gs_T + 2 * F_gw_T + F_gt_T
            Sum[1] = F_ww_T + F_wg_T + F_ws_T + F_wt_T
            Sum[2] = F_sg_T + 2 * F_sw_T + F_st_T
            Sum[3] = F_tg_T + 2 * F_tw_T + F_ts_T + F_tt_T

            Sum2[0] = F_sg_T * w / w + 2 * F_wg_T * h / w + F_tg_T * Atree / w
            Sum2[1] = F_ww_T * h / h + F_gw_T * w / h + F_sw_T * w / h + F_tw_T * Atree / h
            Sum2[2] = F_gs_T * w / w + 2 * F_ws_T * h / w + F_ts_T * Atree / w
            Sum2[3] = F_gt_T * w / Atree + 2 * F_wt_T * h / Atree + F_st_T * w / Atree + F_tt_T * Atree / Atree


        # Check sum
        if a > 0 and any(Sum < 0.9999) or any(Sum > 1.0001): ### Check any
            print('The view factor do not add up to 1. Please check the ray tracing algorithm.')
        elif a == 0 and any(Sum[0:3] < 0.9999) or any(Sum[0:3] > 1.0001):
            print('The view factor do not add up to 1. Please check the ray tracing algorithm.')

        # Assign view factors to struct
        class VFRayTracing_T_Def():
            pass
        VFRayTracing_T = VFRayTracing_T_Def()
        VFRayTracing_T.F_gs_T = F_gs_T
        VFRayTracing_T.F_gt_T = F_gt_T
        VFRayTracing_T.F_gw_T = F_gw_T
        VFRayTracing_T.F_ww_T = F_ww_T
        VFRayTracing_T.F_wt_T = F_wt_T
        VFRayTracing_T.F_wg_T = F_wg_T
        VFRayTracing_T.F_ws_T = F_ws_T
        VFRayTracing_T.F_sg_T = F_sg_T
        VFRayTracing_T.F_sw_T = F_sw_T
        VFRayTracing_T.F_st_T = F_st_T
        VFRayTracing_T.F_tg_T = F_tg_T
        VFRayTracing_T.F_tw_T = F_tw_T
        VFRayTracing_T.F_ts_T = F_ts_T
        VFRayTracing_T.F_tt_T = F_tt_T
        VFRayTracing_T.F_pg = F_pg
        VFRayTracing_T.F_ps = F_ps
        VFRayTracing_T.F_pw = F_pw
        VFRayTracing_T.F_pt = F_pt

        return F_gs_T,F_gt_T,F_gw_T,F_ww_T,F_wt_T,F_wg_T,F_ws_T,F_ts_T,F_tw_T,F_tt_T,F_tg_T,F_sg_T,F_sw_T,F_st_T,F_pg,\
               F_ps,F_pw,F_pt, VFRayTracingRaw_T,VFRayTracing_T

    def VFAnalytical(self,H,W):

        # Sky view factors without trees: Harman et al. 2004
        ratio = H / W

        F_gs_nT = numpy.sqrt(1 + (ratio) ** 2) - ratio
        F_gt_nT = 0
        # factor 0.5 because there are 2 walls that are seen by the ground
        F_gw_nT = 0.5 * (1 - F_gs_nT)

        F_ww_nT = numpy.sqrt(1 + (1 / ratio) ** 2) - 1 / ratio
        F_wt_nT = 0
        F_wg_nT = 0.5 * (1 - F_ww_nT)
        F_ws_nT = 0.5 * (1 - F_ww_nT)

        F_ts_nT = 0
        F_tw_nT = 0
        F_tt_nT = 0
        F_tg_nT = 0

        F_sg_nT = F_gs_nT
        F_sw_nT = ratio * F_ws_nT
        F_st_nT = 0

        # Check for unity of the sum of the view factors
        h = H / W
        w = W / W

        Sum_g = F_gs_nT + F_gt_nT + F_gw_nT * 2
        Sum_w = F_ww_nT + F_wt_nT + F_wg_nT + F_ws_nT
        Sum_t = F_ts_nT + 2 * F_tw_nT + F_tt_nT + F_tg_nT
        Sum_s = F_sg_nT + 2 * F_sw_nT + F_st_nT

        Sum_g2 = F_wg_nT * h / w * 2 + F_sg_nT * w / w
        Sum_w2 = F_gw_nT * w / h + F_ww_nT * h / h + F_sw_nT * w / h
        Sum_t2 = 0
        Sum_s2 = F_gs_nT * w / w + 2 * F_ws_nT * h / w

        class ViewFactor_nT_Def():
            pass
        ViewFactor_nT = ViewFactor_nT_Def()
        ViewFactor_nT.F_gs_nT = F_gs_nT
        ViewFactor_nT.F_gw_nT = F_gw_nT
        ViewFactor_nT.F_ww_nT = F_ww_nT
        ViewFactor_nT.F_wg_nT = F_wg_nT
        ViewFactor_nT.F_ws_nT = F_ws_nT
        ViewFactor_nT.F_sg_nT = F_sg_nT
        ViewFactor_nT.F_sw_nT = F_sw_nT

        return F_gs_nT,F_gt_nT,F_gw_nT,F_ww_nT,F_wt_nT,F_wg_nT,F_ws_nT,F_ts_nT,F_tw_nT,F_tt_nT,F_tg_nT,F_sg_nT,F_sw_nT,\
               F_st_nT,ViewFactor_nT

    def VFRayTracing(self,H,W,a,ht,d,Person,MCSampleSize,NRays):

        # Emitting surface
        # 1 = from wall 1
        # 2 = from wall 2
        # 3 = from ground
        # 4 = from tree 1
        # 5 = from tree 2
        # 6 = from sky
        # 7 = from point p

        View_factor = numpy.zeros((7,6))
        for option_surface in range(7):
            VG, VW1, VW2, VS, VT1, VT2 = self.View_Factors_Geometry(H, W, a, ht, d, Person,option_surface, MCSampleSize,NRays)
            # towards wall 1
            View_factor[option_surface, 0] = VW1
            # towards wall 2
            View_factor[option_surface, 1] = VW2
            # towards ground
            View_factor[option_surface, 2] = VG
            # towards tree 1
            View_factor[option_surface, 3] = VT1
            # towards tree 2
            View_factor[option_surface, 4] = VT2
            # towards sky
            View_factor[option_surface, 5] = VS


        # Elimination of self-view factor and rescaling
        # Wall 1 to wall 1 self view factor elimination
        View_factor[0,:] = [View_factor[0,i]/ (1 - View_factor[0,0]) for i in range(6)]
        View_factor[0,0] = 0
        # Wall 2 to wall 2 self view factor elimination
        View_factor[1,:] = [View_factor[1,i] / (1 - View_factor[1,1]) for i in range(6)]
        View_factor[1,1] = 0
        # Ground to Ground self view factor elimination
        View_factor[2, :] = [View_factor[2, i] / (1 - View_factor[2, 2]) for i in range(6)]
        View_factor[2, 2] = 0
        # Tree 1 to tree 1 self view factor elimination
        View_factor[3, :] = [View_factor[3, i] / (1 - View_factor[3, 3]) for i in range(6)]
        View_factor[3, 3] = 0
        # Tree 2 to tree 2 self view factor elimination
        View_factor[4, :] = [View_factor[4, i] / (1 - View_factor[4, 4]) for i in range(6)]
        View_factor[4, 4] = 0
        # Sky to Sky self view factor elimination
        View_factor[5, :] = [View_factor[5, i] / (1 - View_factor[5, 5]) for i in range(6)]
        View_factor[5, 5] = 0

        # View factor assignment
        F_gs_T = View_factor[2,5]
        F_gt_T = View_factor[2,3] + View_factor[2,4]
        F_gw_T = (View_factor[2,0] + View_factor[2,1]) / 2

        F_ww_T = View_factor[0,1]
        F_wt_T = View_factor[0,3] + View_factor[0,4]
        F_wg_T = View_factor[0,2]
        F_ws_T = View_factor[0,5]

        if a > 0:
            F_ts_T = View_factor[3,5]
            F_tw_T = (View_factor[3,0] + View_factor[3,1]) / 2
            F_tt_T = View_factor[3,4]
            F_tg_T = View_factor[3,2]
        else:
            F_ts_T = 0
            F_tw_T = 0
            F_tt_T = 0
            F_tg_T = 0

        F_sg_T = View_factor[5,2]
        F_sw_T = (View_factor[5,0] + View_factor[5,1]) / 2
        F_st_T = View_factor[5,3] + View_factor[5,4]

        F_pg = View_factor[6,2]
        F_ps = View_factor[6,5]
        F_pw = (View_factor[6,0] + View_factor[6,1]) / 2
        F_pt = View_factor[6,3] + View_factor[6,4]

        # Check sum
        Sum = numpy.zeros(5)
        Sum[0] = F_gs_T + F_gt_T + F_gw_T * 2
        Sum[1] = F_ww_T + F_wt_T + F_wg_T + F_ws_T
        Sum[2] = F_sg_T + 2 * F_sw_T + F_st_T
        Sum[3] = F_pg + F_ps + 2 * F_pw + F_pt
        Sum[4] = F_ts_T + 2 * F_tw_T + F_tt_T + F_tg_T

        if a == 0:
            F_ts_T = 0
            F_tw_T = 0
            F_tt_T = 0
            F_tg_T = 0

        if a > 0 and any(Sum < 0.9999) or any(Sum > 1.0001):
            print('The view factor do not add up to 1. Please check the ray tracing algorithm.')
        elif a == 0 and any(Sum[0:4] < 0.9999) and any(Sum[0:4] > 1.0001):
            print('The view factor do not add up to 1. Please check the ray tracing algorithm.')

        class VFRayTracingRaw_T_Def():
            pass
        VFRayTracingRaw_T = VFRayTracingRaw_T_Def()
        VFRayTracingRaw_T.F_gs_T = F_gs_T
        VFRayTracingRaw_T.F_gt_T = F_gt_T
        VFRayTracingRaw_T.F_gw_T = F_gw_T
        VFRayTracingRaw_T.F_ww_T = F_ww_T
        VFRayTracingRaw_T.F_wt_T = F_wt_T
        VFRayTracingRaw_T.F_wg_T = F_wg_T
        VFRayTracingRaw_T.F_ws_T = F_ws_T
        VFRayTracingRaw_T.F_sg_T = F_sg_T
        VFRayTracingRaw_T.F_sw_T = F_sw_T
        VFRayTracingRaw_T.F_st_T = F_st_T
        VFRayTracingRaw_T.F_tg_T = F_tg_T
        VFRayTracingRaw_T.F_tw_T = F_tw_T
        VFRayTracingRaw_T.F_ts_T = F_ts_T
        VFRayTracingRaw_T.F_tt_T = F_tt_T
        VFRayTracingRaw_T.F_pg = F_pg
        VFRayTracingRaw_T.F_ps = F_ps
        VFRayTracingRaw_T.F_pw = F_pw
        VFRayTracingRaw_T.F_pt = F_pt

        return F_gs_T,F_gt_T,F_gw_T,F_ww_T,F_wt_T,F_wg_T,F_ws_T,F_ts_T,F_tw_T,F_tt_T,F_tg_T,F_sg_T,F_sw_T,F_st_T,F_pg,\
               F_ps,F_pw,F_pt,VFRayTracingRaw_T

    def View_Factors_Geometry(self,H,W,a,ht,d,Person,OPTION_SURFACE,MCSampleSize,NRays):

        # Geometry specification
        h = H / W
        w = W / W

        pz = Person.PositionPz / W
        px = Person.PositionPx / W

        # Roof
        x1a = [0,1]
        z1a = [h,h]

        x1b = [1 + w,2 + w]
        z1b = [h,h]

        # Ground
        x2 = [1,1 + w]
        z2 = [0,0]

        # Wall 1
        x3 = [1,1]
        z3 = [h,0]

        # Wall 2
        x4 = [1 + w,1 + w]
        z4 = [0,h]

        # Sky
        x5 = [1,1 + w]
        z5 = [h,h]

        # Tree 1
        xc = 1 + d * w
        yc = ht * w
        r = a * w
        ang = numpy.arange(start=0,stop=2*numpy.pi,step=0.02)
        xt = [r * math.cos(ang[i]) for i in range(len(ang))]
        yt = [r * math.sin(ang[i]) for i in range(len(ang))]
        if r == 0:
            xc = 0
            yc = 0

        # Tree 2
        xc2 = 1 + w - d * w
        ang = numpy.arange(start=0,stop=2*numpy.pi,step=0.02)
        if r == 0:
            xc2 = 0
            yc = 0

        # Person
        xcp6 = 1 + px
        ycp6 = pz
        rp6 = 1 / 1000
        xp6 = [rp6 * math.cos(ang[i]) for i in range(len(ang))]
        yp6 = [rp6 * math.sin(ang[i]) for i in range(len(ang))]

        # Monte Carlo Parameters
        RandSZ = numpy.random.uniform(0,1,MCSampleSize)
        # Uniformly distributed "random" values in the interval [0,1]
        DeltaRays = numpy.arange(0,1+1/(NRays/2),1/(NRays/2))

        # polar angle (zenith)
        AnlgeDist = [math.asin(DeltaRays[i]) for i in range(len(DeltaRays))]
        # convert it to altitude/elevation angle in first quadrant
        RayAngleQ1 = [numpy.pi/2 - AnlgeDist[i] for i in range(len(AnlgeDist))][::-1]
        # Angle in second quadrant
        RayAngleQ2 = [numpy.pi / 2 + AnlgeDist[i] for i in range(len(AnlgeDist))]
        # for a horizontal planar surface
        RayAngle = RayAngleQ1[0:-1] + RayAngleQ2

        # Ray Angle is defined as the altitude angle on a horizontal surface
        # starting on the "right side" (first quadrat of coordinate system). It can
        # be used directly for the ground surface but needs to be shifted by +pi/2
        # and -pi/2 for the wall surfaces as the walls are vertical in our
        # coordinate system. It also needs to be shifted according to the
        # orientation of the tangent of the emitting point on the tree circle.

        # The emitting point needs to be slightly moved away from the surface.
        # Otherwise, it will be counted as crossing itself. stc defines how much a point is moved away from the surface
        # How far is the starting point away from the surface.
        stc = 10**(-10)

        # Vector definition
        if OPTION_SURFACE == 0:
            print('OPTION_SURFACE: ', OPTION_SURFACE)
            # View Factor from Wall-1
            # Randomly distributed emitting points
            YSv = [h * RandSZ[i] for i in range(len(RandSZ))]
            XSv = [(1+stc)*1 for i in range(len(YSv))]

            # Uniformly distributed emitting points
            RayAngle_array = numpy.array([RayAngle])
            dthe = numpy.ones((len(XSv), 1)) @ (RayAngle_array - numpy.pi / 2)

        elif OPTION_SURFACE == 1:
            print('OPTION_SURFACE: ', OPTION_SURFACE)
            # View Factor from Wall-2
            # Randomly distributed emitting points
            YSv = [h * RandSZ[i] for i in range(len(RandSZ))]
            XSv = [(1 + w - stc) * 1 for i in range(len(YSv))]

            # Uniformly distributed emitting points
            RayAngle_array = numpy.array([RayAngle])
            dthe = numpy.ones((len(XSv), 1)) @ (RayAngle_array + numpy.pi / 2)

        elif OPTION_SURFACE == 2:
            print('OPTION_SURFACE: ', OPTION_SURFACE)
            # View Factor from ground
            # Randomly distributed emitting points
            XSv = [1+w*RandSZ[i] for i in range(len(RandSZ))]
            YSv = [stc*1 for i in range(len(XSv))]

            # Uniformly distributed emitting points
            RayAngle_array = numpy.array([RayAngle])
            dthe = numpy.ones((len(XSv), 1)) @ (RayAngle_array)

        elif OPTION_SURFACE == 3:
            print('OPTION_SURFACE: ', OPTION_SURFACE)
            # View from Tree-1
            # Randomly distributed emitting points
            ang = [2*numpy.pi*RandSZ[i] for i in range(len(RandSZ))]

            # Uniformly distributed emitting points
            xt = [(r + stc) * math.cos(ang[i]) for i in range(len(ang))]
            yt = [(r + stc) * math.sin(ang[i]) for i in range(len(ang))]
            XSv = [xc + xt[i] for i in range(len(xt))]
            YSv = [yc + yt[i] for i in range(len(yt))]
            if r == 0:
                XSv[:] = [0 for k in range(len(XSv))]
                YSv[:] = [0 for k in range(len(YSv))]

            RayAngle_array = numpy.array([RayAngle])
            ang_array = numpy.array([ang])
            dthe = numpy.ones((len(XSv),1)) @ (RayAngle_array-numpy.pi/2) + ang_array.T

        elif OPTION_SURFACE == 4:
            print('OPTION_SURFACE: ', OPTION_SURFACE)
            # View from Tree-2
            # Randomly distributed emitting points
            ang = [2*numpy.pi*RandSZ[i] for i in range(len(RandSZ))]

            # Uniformly distributed emitting points
            xt = [(r + stc) * math.cos(ang[i]) for i in range(len(ang))]
            yt = [(r + stc) * math.sin(ang[i]) for i in range(len(ang))]
            XSv = [xc2 + xt[i] for i in range(len(xt))]
            YSv = [yc + yt[i] for i in range(len(yt))]
            if r == 0:
                XSv[:] = [0 for k in range(len(XSv))]
                YSv[:] = [0 for k in range(len(YSv))]

            RayAngle_array = numpy.array([RayAngle])
            ang_array = numpy.array([ang])
            dthe = numpy.ones((len(XSv), 1)) @ (RayAngle_array - numpy.pi / 2) + ang_array.T

        elif OPTION_SURFACE == 5:
            print('OPTION_SURFACE: ', OPTION_SURFACE)
            # View Factor from sky
            # Randomly distributed emitting points
            XSv = [(1+w*RandSZ[i]) for i in range(len(RandSZ))]
            YSv = [(h-stc)*1 for i in range(len(XSv))]

            # Uniformly distributed emitting points
            RayAngle_array = numpy.array([RayAngle])
            dthe = numpy.ones((len(XSv), 1)) @ (RayAngle_array + numpy.pi)

        elif OPTION_SURFACE == 6:
            print('OPTION_SURFACE: ', OPTION_SURFACE)
            # View from point for MRT
            ang = [2*numpy.pi*RandSZ[i] for i in range(len(RandSZ))]

            # Uniformly distributed emitting points
            xp6 = [(rp6+stc)*math.cos(ang[i]) for i in range(len(ang))]
            yp6 = [(rp6 + stc)*math.sin(ang[i]) for i in range(len(ang))]
            XSv = [xcp6 + xp6[i] for i in range(len(xp6))]
            YSv = [ycp6 + yp6[i] for i in range(len(yp6))]

            RayAngle_array = numpy.array([RayAngle])
            ang_array = numpy.array([ang])
            dthe = numpy.ones((len(XSv), 1)) @ (RayAngle_array - numpy.pi / 2) + ang_array.T

        # Parameters of the search
        # maximum ray length, maximum search distance
        dmax = numpy.sqrt(h ** 2 + w ** 2) + numpy.sqrt(h ** 2 + w ** 2) / 100
        # Search step size for tree detection
        sz = w / 1000
        # plots graph
        GRAPH = 0

        VG, VW1, VW2, VS, VT1, VT2  = self.ViewFactorsComputation(XSv, YSv, dmax, sz, dthe, GRAPH, x2, z2, x3,z3, x4, z4,
                                                                  xc, yc, r, xc2, x5, z5)

        return VG,VW1,VW2,VS,VT1,VT2

    def ViewFactorsComputation(self,XSv,YSv,dmax,sz,dthe,GRAPH,x2,z2,x3,z3,x4,z4,xc,yc,r,xc2,x5,z5):

        # pass of search
        spass = numpy.sqrt(2)*sz
        # search distance [m]
        SD = numpy.arange(start=spass,step=spass,stop=dmax)

        np = len(XSv)
        VGv = numpy.zeros(np)
        VW1v = numpy.zeros(np)
        VW2v = numpy.zeros(np)
        VSv = numpy.zeros(np)
        VT1v = numpy.zeros(np)
        VT2v = numpy.zeros(np)

        # For the number of emitting points
        for ii in range(np):
            # search angle [angular degree]
            Z = copy.copy(dthe[ii,:])

            XS = XSv[ii]
            YS = YSv[ii]
            VG = 0
            VW1 = 0
            VW2 = 0
            VS = 0
            VT1 = 0
            VT2 = 0

            # For the number of rays emitted from each emitting point
            for k in range(len(Z)):

                ## Check
                xp = [SD[ip]*numpy.cos(Z[k]) for ip in range(len(SD))]
                yp = [SD[ip]*numpy.sin(Z[k]) for ip in range(len(SD))]

                # Ground
                l1 = [x2[0],z2[0],x2[1],z2[1]]
                l2 = [XS,YS,XS+xp[-1],YS+yp[-1]]
                out = self.lineSegmentIntersect(l1, l2)
                Sfn = out.intAdjacencyMatrix
                xI = out.intMatrixX
                yI = out.intMatrixY
                if Sfn == 1:
                    # distance
                    D2 = numpy.sqrt(abs(xI-XS)**2 + abs(yI-YS)**2)
                else:
                    D2 = numpy.NaN

                # Wall-1
                l1 = [x3[0],z3[0],x3[1],z3[1]]
                l2 = [XS,YS,XS+xp[-1],YS+yp[-1]]
                out = self.lineSegmentIntersect(l1, l2)
                Sfn = out.intAdjacencyMatrix
                xI = out.intMatrixX
                yI = out.intMatrixY
                if Sfn == 1:
                    # distance
                    D3 = numpy.sqrt(abs(xI-XS)**2 + abs(yI-YS)**2)
                else:
                    D3 = numpy.NaN

                # Wall-2
                l1 = [x4[0], z4[0], x4[1], z4[1]]
                l2 = [XS, YS, XS + xp[-1], YS + yp[-1]]
                out = self.lineSegmentIntersect(l1, l2)
                Sfn = out.intAdjacencyMatrix
                xI = out.intMatrixX
                yI = out.intMatrixY
                if Sfn == 1:
                    # distance
                    D4 = numpy.sqrt(abs(xI - XS) ** 2 + abs(yI - YS) ** 2)
                else:
                    D4 = numpy.NaN

                # Sky
                l1 = [x5[0], z5[0], x5[1], z5[1]]
                l2 = [XS, YS, XS + xp[-1], YS + yp[-1]]
                out = self.lineSegmentIntersect(l1, l2)
                Sfn = out.intAdjacencyMatrix
                xI = out.intMatrixX
                yI = out.intMatrixY
                if Sfn == 1:
                    # distance
                    D5 = numpy.sqrt(abs(xI - XS) ** 2 + abs(yI - YS) ** 2)
                else:
                    D5 = numpy.NaN

                # Tree 1
                # Inside tree
                IC = numpy.zeros(len(xp))
                for i in range(len(xp)):
                    if (XS + xp[i] - xc)**2 + (YS + yp[i] - yc)**2 <= r**2:
                        IC[i] = 1
                    else:
                        IC[i] = 0
                if sum(IC) > 1:

                    sdi = min(numpy.where(IC == 1)[0])
                    DT1 = numpy.sqrt(abs(XS+xp[sdi] - XS) ** 2 + abs(YS+yp[sdi] - YS) ** 2)
                else:
                    DT1 = numpy.NaN

                # Tree 2
                # Inside tree
                IC = numpy.zeros(len(xp))
                for i in range(len(xp)):
                    if (XS + xp[i] - xc2) ** 2 + (YS + yp[i] - yc) ** 2 <= r ** 2:
                        IC[i] = 1
                    else:
                        IC[i] = 0
                if sum(IC) > 1:

                    sdi = min(numpy.where(IC == 1)[0])
                    DT2 = numpy.sqrt(abs(XS + xp[sdi] - XS) ** 2 + abs(YS + yp[sdi] - YS) ** 2)
                else:
                    DT2 = numpy.NaN

                # Assign a count for the surface that the ray is passing through
                # Ground  Wall 1 Wall 2  Tree 1  Tree 2 Sky
                md = min((imin for imin in [D2,D3,D4,DT1,DT2,D5] if not math.isnan(imin)))
                pmin = numpy.where([D2,D3,D4,DT1,DT2,D5] == md)[0]

                if numpy.isnan(md):
                    pass
                else:
                    if pmin == 0:
                        VG = VG + 1
                    elif pmin == 1:
                        VW1 = VW1 + 1
                    elif pmin == 2:
                        VW2 = VW2 + 1
                    elif pmin == 3:
                        VT1 = VT1 + 1
                    elif pmin == 4:
                        VT2 = VT2 + 1
                    else:
                        VS = VS + 1

            # Calculates the view factors for each emitting point
            VG = VG / len(Z)
            VW1 = VW1 / len(Z)
            VW2 = VW2 / len(Z)
            VS = VS / len(Z)
            VT1 = VT1 / len(Z)
            VT2 = VT2 / len(Z)

            # This should be 1
            Sum_view = sum([VG,VW1,VW2,VS,VT1,VT2])

            VGv[ii] = VG / Sum_view
            VW1v[ii] = VW1 / Sum_view
            VW2v[ii] = VW2 / Sum_view
            VSv[ii] = VS / Sum_view
            VT1v[ii] = VT1 / Sum_view
            VT2v[ii] = VT2 / Sum_view

        # Calcualtes the mean view factor of all the emitting points together
        VG = numpy.mean(VGv)
        VW1 = numpy.mean(VW1v)
        VW2 = numpy.mean(VW2v)
        VS = numpy.mean(VSv)
        VT1 = numpy.mean(VT1v)
        VT2 = numpy.mean(VT2v)

        return VG,VW1,VW2,VS,VT1,VT2

    def lineSegmentIntersect(self,XY1,XY2):

        n_rows_1 = 1 ## check
        n_cols_1 = len(XY1)
        n_rows_2 = 1 ## check
        n_cols_2 = len(XY2)

        if n_cols_1 != 4 or n_cols_2 != 4:
            print('Arguments must be a Nx4 matrices.')

        # Prepare matrices for vectorized computation of line intersection points.
        ## check
        X1 = XY1[0]
        X2 = XY1[2]
        Y1 = XY1[1]
        Y2 = XY1[3]
        ## check
        X3 = XY2[0]
        X4 = XY2[2]
        Y3 = XY2[1]
        Y4 = XY2[3]

        X4_X3 = (X4 - X3)
        Y1_Y3 = (Y1 - Y3)
        Y4_Y3 = (Y4 - Y3)
        X1_X3 = (X1 - X3)
        X2_X1 = (X2 - X1)
        Y2_Y1 = (Y2 - Y1)

        numerator_a = X4_X3 * Y1_Y3 - Y4_Y3 * X1_X3
        numerator_b = X2_X1 * Y1_Y3 - Y2_Y1 * X1_X3
        denominator = Y4_Y3 * X2_X1 - X4_X3 * Y2_Y1

        u_a = numerator_a / denominator
        u_b = numerator_b / denominator

        # Find the adjacency matrix A of intersecting lines.
        INT_X = X1 + X2_X1 * u_a
        INT_Y = Y1 + Y2_Y1 * u_a
        if (u_a >= 0) and (u_a <= 1) and (u_b >= 0) and (u_b <= 1):
            INT_B = 1
        else:
            INT_B = 0
        if denominator == 0:
            PAR_B = 1
        else:
            PAR_B = 0
        if (numerator_a == 0 and numerator_b == 0 and PAR_B == 1):
            COINC_B = 1
        else:
            COINC_B = 0

        # Arrange output.
        class out_Def():
            pass
        out = out_Def()
        out.intAdjacencyMatrix = INT_B
        out.intMatrixX = INT_X * INT_B
        out.intMatrixY = INT_Y * INT_B
        out.intNormalizedDistance1To2 = u_a
        out.intNormalizedDistance2To1 = u_b
        out.parAdjacencyMatrix = PAR_B
        out.coincAdjacencyMatrix = COINC_B

        return out

    def TotalSWR_LWR_Rrual(self,RSMParam,Text,MeteoData,SunPosition,simTime):

        # Stefan Boltzmann constant [W m^-2 K^-4]
        SIGMA = 5.67e-08
        # Calculate outgoing and net longwave radiation in the rural area  [W m^-2]
        # Outgoing longwave radiation [W m^-2]
        L_rural_emt = RSMParam.e_rural * SIGMA * Text**4
        # Incoming longwave radiation [W m^-2]
        L_rural_in = RSMParam.e_rural * MeteoData.LWR
        # Net longwave radiation at the rural surface [W m^-2]
        rural_infra = L_rural_in - L_rural_emt

        # Calculate incoming and net shortwave in the rural
        SDir_rural = max(math.cos(SunPosition.theta_Z) * MeteoData.SW_dir, 0.0)
        # Winter: no vegetation
        if simTime.month < RSMParam.vegStart or simTime.month > RSMParam.vegEnd:
            S_rural = (1 - RSMParam.a_rural) * (SDir_rural + MeteoData.SW_diff)
            S_rural_out = RSMParam.a_rural * (SDir_rural + MeteoData.SW_diff)
        # Summer: effect of vegetation is considered
        else:
            S_rural = ((1-RSMParam.rurVegCover)*(1-RSMParam.a_rural)+RSMParam.rurVegCover*(1-RSMParam.aveg_rural))*(SDir_rural+MeteoData.SW_diff)
            S_rural_out = ((1-RSMParam.rurVegCover)*RSMParam.a_rural+RSMParam.rurVegCover*RSMParam.aveg_rural)*(SDir_rural+MeteoData.SW_diff)


        rural_solRec = SDir_rural + MeteoData.SW_diff
        rural_solAbs = S_rural

        class SWR_Rural_Def():
            pass
        SWR_Rural = SWR_Rural_Def()
        SWR_Rural.SWRinRural = rural_solRec
        SWR_Rural.SWRoutRural = S_rural_out
        SWR_Rural.SWRabsRural = rural_solAbs
        SWR_Rural.Dir = SDir_rural
        SWR_Rural.Diff = MeteoData.SW_diff
        class LWR_Rural_Def():
            pass
        LWR_Rural = LWR_Rural_Def()
        LWR_Rural.LWRinRural = L_rural_in
        LWR_Rural.LWRoutRural = L_rural_emt
        LWR_Rural.LWRabsRural = rural_infra

        return SWR_Rural,LWR_Rural



























































