art_make(
         LIB_LIBRARIES 
	 #larcorealg_Geometry
         larreco_RecoAlg_ClusterRecoUtil
         lardataobj_RecoBase
         ${ROOT_BASIC_LIB_LIST}
         ${ROOT_SPECTRUM}
	
	 LIBRARY_NAME
	 uboonecode_uboone_EventWeightReader

         MODULE_LIBRARIES
                        larcore_Geometry_Geometry_service
                        larsim_Simulation nutools_ParticleNavigation lardataobj_Simulation
                        lardata_Utilities
                        larevt_Filters
                        lardataobj_RawData
                        lardataobj_RecoBase
                        lardata_RecoObjects larreco_Deprecated  
                        larreco_RecoAlg
                        larpandora_LArPandoraInterface
                        nusimdata_SimulationBase
                        larreco_ClusterFinder
                        larreco_Calorimetry
                        larcoreobj_SummaryData
                        lardata_ArtDataHelper
                        nutools_RandomUtils_NuRandomService_service
			uboonecode_uboone_EventWeightReader

                        ${LARRECO_LIB}
                        ${LARDATA_LIB}
                        ${ART_FRAMEWORK_CORE}
                        ${ART_FRAMEWORK_PRINCIPAL}
                        ${ART_FRAMEWORK_SERVICES_REGISTRY}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
                        #art_Persistency_Common canvas
                        #art_Persistency_Provenance canvas
                        #art_Utilities canvas
                        ${MF_MESSAGELOGGER}
                        ${MF_UTILITIES}
                        ${FHICLCPP}
                        #cetlib cetlib_except
                        ${ROOT_GEOM}
                        ${ROOT_XMLIO}
                        ${ROOT_GDML}
                        ${ROOT_BASIC_LIB_LIST}
        )

install_headers()
install_fhicl()
install_source()
