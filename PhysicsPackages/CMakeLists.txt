
# Add Source Files to the global list
add_sources( ./physicsPackage_inter.f90
             ./physicsPackageFactory_func.f90
             ./eigenPhysicsPackage_class.f90
	     ./fixedSourcePhysicsPackage_class.f90
             ./vizPhysicsPackage_class.f90
             ./rayVolPhysicsPackage_class.f90
	     ./randomRayPhysicsPackage_class_org.f90
	     ./dynamicRRPhysicsPackage_class_TI.f90
	     ./dynamicRRPhysicsPackage_class_TCR.f90
	     ./transient_class.f90 # TODO
             )
             #./dynamPhysicsPackage_class.f90)
