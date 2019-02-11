#Lengths are in nm

#Defined sx and sy to be the simulation size without the PML(opposite of how David is defining it)
sx=150
sy=155.64

material =

geometry = [mp.Block(mp.Vector3(29.16, 97.32, 0),
                     center=mp.Vector3(-60.42, 0),
                     material=material),
            mp.Block(mp.Vector3(29.16, 97.32, 0),
                     center=mp.Vector3(60.42, 0),
                     material=material),
            mp.Block(mp.Vector3(20, 9.16, 0),
                     center=mp.Vector3(-65, -53.24),
                     material=material),
            mp.Block(mp.Vector3(20, 9.16, 0),
                     center=mp.Vector3(65, -53.24),
                     material=material),
            mp.Cylinder(radius=9.16,
                        height=mp.inf,
                        axis=mp.Vector3(0, 0, 1),
                        center=mp.Vector3(-55, -48.66),
                        material=material),
            mp.Cylinder(radius=9.16,
                        height=mp.inf,
                        axis=mp.Vector3(0, 0, 1),
                        center=mp.Vector3(55, -48.66),
                        material=material),
            mp.Prism(vertices=[mp.Vector3(0,0,0), mp.Vector3(8.98, 0, 0), mp.Vector3(0, -97.32, 0)],
                     height=mp.inf,
                     axis=mp.Vector3(0, 0, 1),
                     center=mp.Vector3(-42.85, 16.22, 0),
                     material=material),
            mp.Prism(vertices=[mp.Vector3(0,0,0), mp.Vector3(8.98, 0, 0), mp.Vector3(8.98, -97.32, 0)],
                     height=mp.inf,
                     axis=mp.Vector3(0, 0, 1),
                     center=mp.Vector3(42.85, 16.22, 0),
                     material=material),
            mp.Block(mp.Vector3(150, 29.16, 0),
                     center=mp.Vector3(0, 48.66+(29.16)/2),
                     material=material),
            mp.Ellipsoid(center=mp.Vector3(0, 48.66, 0),
                         size=mp.Vector3(73.72, 26.04, mp.inf),
                         material=mp.Medium(epsilon=1))]
