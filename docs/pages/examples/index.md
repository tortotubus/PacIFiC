# Examples

## Preview Geometry

### Unit Cube

Prepare a .xml file for use with Grains3D using the *CompositeRigidBody* setting. 

**cube.xml**
```
<Grains3D Type="CompositeRigidBody">
  <Construction>
    <LinkedCell MX="1." MY="1." MZ="1."/>
    <CompositeParticles>
      <CompositeParticle Number="1" Density="1000.">
        <Volume Value="0."/>
        <MomentOfInertiaTensor>
          <Ixx Value="0."/>
          <Ixy Value="0."/>
          <Ixz Value="0."/>
          <Iyy Value="0."/>
          <Iyz Value="0."/>
          <Izz Value="0."/>
        </MomentOfInertiaTensor>
        <AngularPosition Type="Identity"></AngularPosition>
        <Material>matP</Material>
        <ElementaryParticles>
          <ElementaryParticle>
            <Convex CrustThickness="2e-05">
              <Box LX="1e-4" LY="1e-4" LZ="1e-4" />
            </Convex>
            <RelativePosition X="0" Y="0" Z="0"/>
            <AngularPosition Type="Identity"></AngularPosition>
          </ElementaryParticle>
        </ElementaryParticles>
      </CompositeParticle>
    </CompositeParticles>
    <CompositeParticleFeatures GeomTypeID="0">
      <OutputFile Name="composite"/>
      <Grid NonDimDelta="1.e-3" Nmax="32"/>
    </CompositeParticleFeatures>
  </Construction>
</Grains3D>
```

To run, type `grains3d cube.xml`

### Dendrite Snowflake

Prepare a .xml file for use with Grains3D using the "CompositeRigidBody" setting. 
 
 **dendrite.xml**

 ```
 <Grains3D Type="CompositeRigidBody">
  <Construction>
    <LinkedCell MX="1." MY="1." MZ="1."/>
    <CompositeParticles>
      <CompositeParticle Number="1" Density="2410." SpecificShape="Dendrite">
        <Geometry ArmWidth="2e-3" Depth="1e-3" ArmLength="4e-3" CrustThickness="3.914e-05"/>
        <AngularPosition Type="Identity"></AngularPosition>
        <Material>matP</Material>
      </CompositeParticle>   
    </CompositeParticles>
    <CompositeParticleFeatures GeomTypeID="0">
      <OutputFile Name="composite"/>
      <Grid NonDimDelta="1.e-3" Nmax="64"/>
    </CompositeParticleFeatures>
  </Construction>
</Grains3D>
 ```

 To run, type `grains3d dendrite.xml`.

 ## Run Granular Simulations

 ### Falling Snowflakes

 The following script simulates 40 falling dendrite snowflakes using the *Standard* setting in Grains3d.

**sim.xml**
 ```
 <Grains3D Type="Standard">
  <Construction>
    <LinkedCell MX="0.124" MY="0.124" MZ="0.615" CellSizeFactor="1."/>
    <Origin OX="0." OY="0." OZ="0."/>
    <Periodicity PX="0" PY="0" PZ="0"/>  
    <DomainDecomposition NX="1" NY="1" NZ="1"/>
    <Verbosity Level="1"/> 


    <ContactForceModels>
      <ContactForceModel>
	<Material materialA="matP" materialB="matP"/>
        <HookeMemory>
      	  <kn>2.e+7</kn>	  
	  <en>0.85</en>
	  <muc>0.5</muc>
	  <kt>2.e+7</kt>
	  <etarpf>0.</etarpf>
	  <mur>0.01</mur> 
	</HookeMemory>	  
      </ContactForceModel>   

      <ContactForceModel>
	<Material materialA="matP" materialB="matS"/>
        <HookeMemory>
      	  <kn>2.e+7</kn>	  
	  <en>0.85</en>
	  <muc>0.</muc>
	  <kt>2.e+7</kt>
	  <etarpf>0.</etarpf>
	  <mur>0.</mur> 
	</HookeMemory>	  
      </ContactForceModel>                      
    </ContactForceModels>    

    <CompositeParticles>                        
      <CompositeParticle Number="40" Density="2410." SpecificShape="Dendrite">
        <Geometry ArmWidth="2e-3" Depth="1e-3" ArmLength="4e-3" CrustThickness="3.914e-05"/>
        <AngularPosition Type="Identity"></AngularPosition>
        <Material>matP</Material>
      </CompositeParticle>                           
    </CompositeParticles>     

    <Obstacles>
      <CylindricalShell name="UpperCyl">
          <Geometry Height="0.372" Radius="0.062" Width="0.005" N="32" CrustThickness="3.914e-05" ElementaryObstacle="Box"/>
          <Centre X="0.062" Y="0.062" Z="0.426"/>
        <AngularPosition Type="Identity"></AngularPosition>
        <Material>matP</Material>
      </CylindricalShell>
      <Obstacle name="LowerBottom">
          <Convex CrustThickness="3.914e-05">
            <Rectangle LX="0.149" LY="0.149"/>
          </Convex>
          <Centre X="0.062" Y="0.062" Z="0.24"/>
          <AngularPosition Type="Identity"></AngularPosition>
          <Material>matP</Material>
      </Obstacle>       
      <Obstacle name="UpperBottom">
          <Convex CrustThickness="3.914e-05">
            <Rectangle LX="0.149" LY="0.149"/>
          </Convex>
          <Centre X="0.062" Y="0.062" Z="0.615"/>
          <AngularPosition Type="Identity"></AngularPosition>
          <Material>matP</Material>
      </Obstacle>                       
    </Obstacles>   

    <CollisionDetectionAlgorithm>
	<CollisionDetection Type="GJK" Tolerance="1e-08" Acceleration="OFF" History="OFF"/>
	<BoundingVolume Type="OBB"/>
    </CollisionDetectionAlgorithm>  

  </Construction>


  <Forces>
    <Gravity GX="0." GY="0." GZ="-9.807"/>
  </Forces>


  <Simulation>
    <TimeInterval Start="0." End="1.2"/> 
    <TimeStep dt="1.e-6"/>
    <RestartFile Name="Restart/insert" WritingMode="Hybrid" SingleFile="True" Clones="None"/>
    <TimeSave Start="0." End="10." Every="0.01"/>         
    <CollisionalRelativeVelocity value="1.911"/> 
    <MovingObstacles LinkUpdateEvery="1" GeometricallyMove="True"/>    

    <ParticleInsertion>
      <Mode Type="OverTime"/>
      <Order Type="Random"/>
      <InitialAngularPosition Type="Random"/>
      <RandomGeneratorSeed Type="Random"/> 
      <Frequency TryEvery="1"/> 
      <ForceInsertion Value="False"/>    

      <ParticlePosition>
        <Windows>
          <Window Type="Cylinder" Name="IZone">
            <BottomCentre X="0.062" Y="0.062" Z="0.3"/>
            <Cylinder Radius="0.049" Height="0.0" Direction="Z"/>      
          </Window>            
        </Windows>
      </ParticlePosition>  

      <InitialVelocity Mode="Zero">
      </InitialVelocity>
    </ParticleInsertion>

    <PostProcessing>   
      <Writers>
        <Paraview RootName="insert" Directory="Post" WritingMode="Binary" PerType="True" ForceNetwork="True" MPIIO="True"/>
        <RawData Name="Post/insert" WritingMode="Binary"/>
      </Writers>              
    </PostProcessing>
  </Simulation>      
</Grains3D>
 ```

 To run, type `grains3d sim.xml`.